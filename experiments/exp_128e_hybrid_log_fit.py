"""Sprint 128e: Test logarithmic growth for hybrid chi_F at q>=5.

Hybrid alpha extrapolates to 0 for q>=5 (Sprint 128c).
Test whether chi_F ~ A*(ln N)^beta fits better than power law.

If logarithmic wins: flag as potentially novel finding.
If it loses: close the hybrid thread.
"""
import numpy as np
import json, time, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from lmfit import Model
from db_utils import query, record
from fss_utils import fit_power_law, pairwise_exponents

results = {
    'experiment': '128e_hybrid_log_fit',
    'sprint': 128,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'fits': {},
    'verdict': None,
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results',
                           'sprint_128e_hybrid_log_fit.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)


def power_law(N, A, alpha):
    return A * N**alpha

def log_growth(N, A, beta):
    return A * np.log(N)**beta

def log_linear(N, A, B):
    """chi_F = A * ln(N) + B (simplest sub-power-law)"""
    return A * np.log(N) + B


print("=" * 80)
print("Sprint 128e: Hybrid chi_F -- power law vs logarithmic growth")
print("=" * 80)

q_values = [3, 4, 5, 6, 7]
log_wins = 0
power_wins = 0
total_tested = 0

for q_val in q_values:
    rows = query(quantity='chi_F_exact', model='hybrid', q=q_val)
    if not rows:
        continue

    # Build data (most recent sprint per n)
    data = {}
    for r in rows:
        n, val, sprint = r[4], r[6], r[1]
        if n is not None and val is not None:
            if n not in data or sprint > data[n][1]:
                data[n] = (val, sprint)

    sizes = np.array(sorted(data.keys()), dtype=float)
    chi_F = np.array([data[int(n)][0] for n in sizes])

    if len(sizes) < 4:
        print(f"\n  q={q_val}: only {len(sizes)} sizes, skipping")
        continue

    total_tested += 1
    print(f"\n{'='*70}")
    print(f"  HYBRID q={q_val}: {len(sizes)} sizes")
    print(f"{'='*70}")
    for i, n in enumerate(sizes):
        print(f"    n={int(n):2d}: chi_F = {chi_F[i]:.6f}")

    # Pairwise exponents
    pairs = pairwise_exponents(sizes, chi_F)
    print(f"  Pairwise alpha:")
    for p in pairs:
        print(f"    ({p['n1']},{p['n2']}): {p['alpha']:.4f}")

    # Fit 1: Power law
    m1 = Model(power_law)
    p1 = m1.make_params(A=1.0, alpha=1.0)
    p1['A'].min = 0.001
    r1 = m1.fit(chi_F, p1, N=sizes)
    ss_res_1 = np.sum(r1.residual**2)
    ss_tot = np.sum((chi_F - chi_F.mean())**2)
    r2_1 = 1 - ss_res_1 / ss_tot

    # Fit 2: Logarithmic power
    m2 = Model(log_growth)
    p2 = m2.make_params(A=1.0, beta=2.0)
    p2['A'].min = 0.001
    p2['beta'].min = 0.1
    r2_fit = m2.fit(chi_F, p2, N=sizes)
    ss_res_2 = np.sum(r2_fit.residual**2)
    r2_2 = 1 - ss_res_2 / ss_tot

    # Fit 3: Simple log-linear
    m3 = Model(log_linear)
    p3 = m3.make_params(A=3.0, B=-1.0)
    r3 = m3.fit(chi_F, p3, N=sizes)
    ss_res_3 = np.sum(r3.residual**2)
    r2_3 = 1 - ss_res_3 / ss_tot

    print(f"\n  [1] Power law: A={r1.params['A'].value:.4f}, alpha={r1.params['alpha'].value:.4f}")
    print(f"      R2={r2_1:.8f}, AIC={r1.aic:.2f}")

    print(f"  [2] Log power: A={r2_fit.params['A'].value:.4f}, beta={r2_fit.params['beta'].value:.4f}")
    print(f"      R2={r2_2:.8f}, AIC={r2_fit.aic:.2f}")

    print(f"  [3] Log-linear: A={r3.params['A'].value:.4f}, B={r3.params['B'].value:.4f}")
    print(f"      R2={r2_3:.8f}, AIC={r3.aic:.2f}")

    # Determine winner
    all_fits = [
        ('power_law', r1.aic, r2_1),
        ('log_power', r2_fit.aic, r2_2),
        ('log_linear', r3.aic, r2_3),
    ]
    all_fits.sort(key=lambda x: x[1])
    winner = all_fits[0][0]
    daic = all_fits[1][1] - all_fits[0][1]

    print(f"\n  WINNER: {winner} (dAIC to next = {daic:.2f})")

    if 'log' in winner:
        log_wins += 1
    else:
        power_wins += 1

    results['fits'][f'q{q_val}'] = {
        'sizes': [int(s) for s in sizes],
        'chi_F': [float(v) for v in chi_F],
        'power_law': {
            'alpha': float(r1.params['alpha'].value),
            'alpha_err': float(r1.params['alpha'].stderr) if r1.params['alpha'].stderr else None,
            'A': float(r1.params['A'].value),
            'r_squared': float(r2_1),
            'aic': float(r1.aic),
        },
        'log_power': {
            'beta': float(r2_fit.params['beta'].value),
            'beta_err': float(r2_fit.params['beta'].stderr) if r2_fit.params['beta'].stderr else None,
            'A': float(r2_fit.params['A'].value),
            'r_squared': float(r2_2),
            'aic': float(r2_fit.aic),
        },
        'log_linear': {
            'A': float(r3.params['A'].value),
            'B': float(r3.params['B'].value),
            'r_squared': float(r2_3),
            'aic': float(r3.aic),
        },
        'winner': winner,
        'daic': float(daic),
        'pairwise': [{'n1': p['n1'], 'n2': p['n2'], 'alpha': p['alpha']} for p in pairs],
    }
    save()


# Verdict
print(f"\n{'='*80}")
print("VERDICT")
print(f"{'='*80}")
print(f"\n  Log wins: {log_wins}/{total_tested}")
print(f"  Power wins: {power_wins}/{total_tested}")

if log_wins > power_wins:
    verdict = "LOG_WINS"
    print(f"\n  Logarithmic growth fits better for hybrid q>=5.")
    print(f"  chi_F ~ A*(ln N)^beta is the better model.")
    print(f"  This confirms hybrid transitions are NOT walking (sub-power-law).")
    print(f"  Decision: one more sprint to harden, then close hybrid thread.")
elif power_wins >= total_tested - 1:
    verdict = "POWER_WINS"
    print(f"\n  Power law still wins for most q values.")
    print(f"  Hybrid transitions have power-law chi_F with small exponents.")
    print(f"  Decision: close hybrid thread, redirect to S_q q=4 and q=8.")
else:
    verdict = "MIXED"
    print(f"\n  Mixed results -- no clear winner.")
    print(f"  Decision: close hybrid thread (insufficient signal).")

results['verdict'] = verdict
results['log_wins'] = log_wins
results['power_wins'] = power_wins

save()
print(f"\nResults saved.")
