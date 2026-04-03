"""Sprint 109c: Power-law vs logarithmic correction fits for alpha(q).

Compare two models for pairwise alpha drift:
  Model A (log):   alpha_pair(N) = alpha_inf + c/ln(N)
  Model B (power): alpha_pair(N) = alpha_inf + c/N^omega

For q=3 (continuous): expect power-law with omega ~ 4/5
For q=4 (BKT): expect logarithmic
For q=5 (walking): expect no correction

Test: does power-law q=3 extrapolation recover exact nu=5/6 (alpha=1.4)?
"""
import numpy as np
import json, time, os
from scipy.optimize import curve_fit
from db_utils import record

results = {
    'experiment': '109c_correction_analysis',
    'sprint': 109,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_109c_correction_analysis.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

# Load all data from 109a, 109b, and existing Sprint 108 data
def load_data():
    """Load chi_F data for q=2,3,4,5."""
    data = {}

    # q=3: from 109a
    with open('results/sprint_109a_q3_gpu_extend.json') as f:
        d = json.load(f)
    q3_n = []; q3_chi = []; q3_gap = []; q3_me = []
    for key in sorted(d['data'].keys(), key=lambda x: int(x)):
        dd = d['data'][key]
        q3_n.append(dd['n']); q3_chi.append(dd['chi_F'])
        q3_gap.append(dd['gap_multiplet']); q3_me.append(dd['me_sq'])
    data[3] = {'n': np.array(q3_n, dtype=float), 'chi': np.array(q3_chi),
               'gap': np.array(q3_gap), 'me': np.array(q3_me)}

    # q=2 and q=4: from 109b
    with open('results/sprint_109b_q4_q2_extend.json') as f:
        d = json.load(f)
    for q_val in [2, 4]:
        ns = []; chis = []; gaps = []; mes = []
        for key in sorted(d['data'].keys(), key=lambda k: d['data'][k]['n']):
            dd = d['data'][key]
            if dd['q'] == q_val:
                ns.append(dd['n']); chis.append(dd['chi_F'])
                gaps.append(dd['gap_multiplet']); mes.append(dd['me_sq'])
        data[q_val] = {'n': np.array(ns, dtype=float), 'chi': np.array(chis),
                       'gap': np.array(gaps), 'me': np.array(mes)}

    # q=5: from Sprint 108b
    with open('results/sprint_108b_log_corrections.json') as f:
        d = json.load(f)
    d5 = d['data'].get('q=5', {})
    if d5:
        data[5] = {'n': np.array(d5['Ns'], dtype=float), 'chi': np.array(d5['chi_Fs']),
                   'gap': np.array(d5['gap_ms']), 'me': np.array(d5['me_sqs'])}

    return data

def pairwise_exponents(Ns, values):
    """Compute pairwise log-log slopes."""
    log_N = np.log(Ns)
    log_v = np.log(values)
    pw = []
    N_mid = []
    for i in range(len(Ns)-1):
        dln = log_N[i+1] - log_N[i]
        pw.append((log_v[i+1] - log_v[i]) / dln)
        N_mid.append(0.5 * (Ns[i] + Ns[i+1]))
    return np.array(pw), np.array(N_mid)

def model_log(N, a_inf, c_log):
    return a_inf + c_log / np.log(N)

def model_power(N, a_inf, c_pow, omega):
    return a_inf + c_pow / N**omega

print("=" * 70)
print("Sprint 109c: Power-law vs logarithmic correction analysis")
print("=" * 70)

data = load_data()

exact_alpha = {2: 1.0, 3: 1.4, 4: 2.0}  # alpha = 2/nu - 1; nu: q=2→1, q=3→5/6, q=4→2/3
exact_nu = {2: 1.0, 3: 5/6, 4: 2/3}

for q in [2, 3, 4, 5]:
    d = data[q]
    Ns = d['n']
    chis = d['chi']

    pw_alpha, N_mid = pairwise_exponents(Ns, chis)

    print(f"\n{'=' * 60}")
    print(f"q={q}: {len(Ns)} sizes, N={int(Ns[0])}-{int(Ns[-1])}")
    print(f"{'=' * 60}")
    print(f"  Pairwise alpha: {['%.4f' % a for a in pw_alpha]}")

    # Skip first 2 pairs (small-N artifacts)
    skip = 2
    pw_fit = pw_alpha[skip:]
    N_fit = N_mid[skip:]

    if len(pw_fit) < 3:
        print(f"  Not enough pairs after skipping first {skip}")
        results['data'][f'q={q}'] = {
            'Ns': [int(n) for n in Ns],
            'pairwise_alpha': [float(a) for a in pw_alpha],
            'skip_note': 'not enough data for fit',
        }
        continue

    # Model A: alpha_pair = a_inf + c/ln(N)
    try:
        popt_log, pcov_log = curve_fit(model_log, N_fit, pw_fit, p0=[pw_fit[-1], 1.0])
        pred_log = model_log(N_fit, *popt_log)
        ss_res_log = np.sum((pw_fit - pred_log)**2)
        ss_tot = np.sum((pw_fit - np.mean(pw_fit))**2)
        R2_log = 1 - ss_res_log / ss_tot if ss_tot > 1e-15 else 0
        aic_log = len(pw_fit) * np.log(ss_res_log / len(pw_fit) + 1e-30) + 2 * 2
        a_inf_log = popt_log[0]
        c_log = popt_log[1]
    except:
        R2_log = -1; a_inf_log = np.nan; c_log = np.nan; aic_log = 1e10

    # Model B: alpha_pair = a_inf + c/N^omega (3 params)
    try:
        popt_pow, pcov_pow = curve_fit(model_power, N_fit, pw_fit,
                                        p0=[pw_fit[-1], 1.0, 0.8],
                                        bounds=([0, -10, 0.1], [5, 10, 3.0]),
                                        maxfev=5000)
        pred_pow = model_power(N_fit, *popt_pow)
        ss_res_pow = np.sum((pw_fit - pred_pow)**2)
        R2_pow = 1 - ss_res_pow / ss_tot if ss_tot > 1e-15 else 0
        aic_pow = len(pw_fit) * np.log(ss_res_pow / len(pw_fit) + 1e-30) + 2 * 3
        a_inf_pow = popt_pow[0]
        c_pow = popt_pow[1]
        omega_fit = popt_pow[2]
    except Exception as e:
        print(f"  Power-law fit failed: {e}")
        R2_pow = -1; a_inf_pow = np.nan; c_pow = np.nan; omega_fit = np.nan; aic_pow = 1e10

    # Model C: alpha_pair = a_inf + c/N^omega with omega FIXED at known value
    omega_known = {2: 2.0, 3: 0.8, 4: None, 5: None}
    if omega_known.get(q) is not None:
        om = omega_known[q]
        try:
            def model_fixed(N, a_inf, c_fix):
                return a_inf + c_fix / N**om
            popt_fix, pcov_fix = curve_fit(model_fixed, N_fit, pw_fit, p0=[pw_fit[-1], 1.0])
            pred_fix = model_fixed(N_fit, *popt_fix)
            ss_res_fix = np.sum((pw_fit - pred_fix)**2)
            R2_fix = 1 - ss_res_fix / ss_tot if ss_tot > 1e-15 else 0
            aic_fix = len(pw_fit) * np.log(ss_res_fix / len(pw_fit) + 1e-30) + 2 * 2
            a_inf_fix = popt_fix[0]
        except:
            R2_fix = -1; a_inf_fix = np.nan; aic_fix = 1e10
    else:
        R2_fix = -1; a_inf_fix = np.nan; aic_fix = 1e10

    print(f"\n  Model A (log):   alpha_inf={a_inf_log:.4f}, c={c_log:.3f}, R²={R2_log:.6f}, AIC={aic_log:.2f}")
    print(f"  Model B (power): alpha_inf={a_inf_pow:.4f}, c={c_pow:.3f}, omega={omega_fit:.3f}, R²={R2_pow:.6f}, AIC={aic_pow:.2f}")
    if omega_known.get(q) is not None:
        print(f"  Model C (fixed omega={omega_known[q]:.1f}): alpha_inf={a_inf_fix:.4f}, R²={R2_fix:.6f}, AIC={aic_fix:.2f}")

    if q in exact_alpha:
        print(f"\n  Exact alpha = {exact_alpha[q]:.4f} (nu = {exact_nu[q]:.4f})")
        print(f"  Log  extrapolation: {a_inf_log:.4f} (error {abs(a_inf_log - exact_alpha[q]):.4f}, {abs(a_inf_log - exact_alpha[q])/exact_alpha[q]*100:.1f}%)")
        print(f"  Power extrapolation: {a_inf_pow:.4f} (error {abs(a_inf_pow - exact_alpha[q]):.4f}, {abs(a_inf_pow - exact_alpha[q])/exact_alpha[q]*100:.1f}%)")
        if omega_known.get(q) is not None:
            print(f"  Fixed-omega extrap: {a_inf_fix:.4f} (error {abs(a_inf_fix - exact_alpha[q]):.4f}, {abs(a_inf_fix - exact_alpha[q])/exact_alpha[q]*100:.1f}%)")

        nu_log = 2.0 / (a_inf_log + 1) if a_inf_log > -1 else np.nan
        nu_pow = 2.0 / (a_inf_pow + 1) if a_inf_pow > -1 else np.nan
        print(f"\n  nu from log:   {nu_log:.4f} (exact {exact_nu[q]:.4f})")
        print(f"  nu from power: {nu_pow:.4f} (exact {exact_nu[q]:.4f})")

    # Also fit z_m and beta_me pairwise
    pw_z, _ = pairwise_exponents(Ns, 1.0/d['gap'])  # gap ~ N^{-z} => 1/gap ~ N^z
    pw_b, _ = pairwise_exponents(Ns, d['me'])

    print(f"\n  Pairwise z_m:     {['%.4f' % z for z in pw_z]}")
    print(f"  Pairwise beta_me: {['%.4f' % b for b in pw_b]}")

    results['data'][f'q={q}'] = {
        'Ns': [int(n) for n in Ns],
        'chi_F': [float(c) for c in chis],
        'pairwise_alpha': [float(a) for a in pw_alpha],
        'pairwise_z_m': [float(z) for z in pw_z],
        'pairwise_beta_me': [float(b) for b in pw_b],
        'log_fit': {'alpha_inf': float(a_inf_log), 'c': float(c_log), 'R2': float(R2_log), 'AIC': float(aic_log)},
        'power_fit': {'alpha_inf': float(a_inf_pow), 'c': float(c_pow), 'omega': float(omega_fit), 'R2': float(R2_pow), 'AIC': float(aic_pow)},
    }
    if omega_known.get(q) is not None:
        results['data'][f'q={q}']['fixed_omega_fit'] = {
            'omega': omega_known[q], 'alpha_inf': float(a_inf_fix), 'R2': float(R2_fix), 'AIC': float(aic_fix)
        }

save()

# Summary table
print(f"\n{'=' * 70}")
print("SUMMARY: Correction model comparison")
print(f"{'=' * 70}")
print(f"{'q':>3} {'exact_a':>8} {'a_log':>8} {'a_pow':>8} {'R2_log':>8} {'R2_pow':>8} {'AIC_log':>8} {'AIC_pow':>8} {'omega':>6} {'winner':>8}")
print("-" * 75)
for q in [2, 3, 4, 5]:
    d = results['data'].get(f'q={q}', {})
    if not d:
        continue
    ea = exact_alpha.get(q, '?')
    dl = d.get('log_fit', {})
    dp = d.get('power_fit', {})
    winner = "POWER" if dp.get('AIC', 1e10) < dl.get('AIC', 1e10) else "LOG"
    if q == 5:
        winner = "NONE"
    print(f"{q:>3} {ea if isinstance(ea, str) else f'{ea:.3f}':>8} "
          f"{dl.get('alpha_inf', np.nan):>8.3f} {dp.get('alpha_inf', np.nan):>8.3f} "
          f"{dl.get('R2', np.nan):>8.4f} {dp.get('R2', np.nan):>8.4f} "
          f"{dl.get('AIC', np.nan):>8.1f} {dp.get('AIC', np.nan):>8.1f} "
          f"{dp.get('omega', np.nan):>6.2f} {winner:>8}")

# Record key results
for q in [2, 3, 4]:
    d = results['data'].get(f'q={q}', {})
    dp = d.get('power_fit', {})
    dl = d.get('log_fit', {})
    if dp.get('alpha_inf') is not None:
        record(sprint=109, model='sq_potts', q=q, n=0,
               quantity='alpha_inf_power', value=dp['alpha_inf'],
               method='power_law_fit', notes=f"omega={dp.get('omega', '?')}")
    if dl.get('alpha_inf') is not None:
        record(sprint=109, model='sq_potts', q=q, n=0,
               quantity='alpha_inf_log', value=dl['alpha_inf'],
               method='log_fit')

save()
print("\nResults saved.")
