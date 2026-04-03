#!/usr/bin/env python3
"""Sprint 089b: Level redistribution analysis across q=2,3,5,7.

Key hypothesis: walking-specific effects live in HOW weight redistributes
between ground state, (q-1) multiplet, and tail — not in tail growth rate
(which is universal b≈2.0).

Compute d(%S)/d(ln n) and d(w)/d(ln n) for each level.
Compare redistribution trajectories across q.
"""
import numpy as np
import json, time
from scipy.optimize import curve_fit

results = {
    'experiment': '089b_level_redistribution',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
    'derivatives': {},
    'discriminators': {},
}

def save():
    with open("results/sprint_089b_level_redistribution.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# Load all DMRG entanglement spectrum data
files = {
    2: "results/sprint_088b_entspec_dmrg_q2.json",
    3: "results/sprint_088a_entspec_dmrg_q3.json",
    5: "results/sprint_087a_entspec_dmrg_q5.json",
    7: ["results/sprint_087b_entspec_dmrg_q7.json",
        "results/sprint_089a_entspec_dmrg_q7_ext.json"],
}

data = {}
for q_val, fnames in files.items():
    if isinstance(fnames, str):
        fnames = [fnames]
    all_entries = []
    for fname in fnames:
        with open(fname) as f:
            d = json.load(f)
        all_entries.extend(d['data'])
    # Deduplicate by n (keep later entry if duplicate)
    seen = {}
    for entry in all_entries:
        seen[entry['n']] = entry
    entries = sorted(seen.values(), key=lambda x: x['n'])

    data[q_val] = {
        'ns': np.array([e['n'] for e in entries]),
        'w_max': np.array([e['lam_max'] for e in entries]),
        'w_mult': np.array([e['w_mult'] for e in entries]),
        'w_tail': np.array([e['w_tail'] for e in entries]),
        'S_lev0': np.array([e['S_lev0_frac'] for e in entries]),
        'S_lev1': np.array([e['S_lev1_frac'] for e in entries]),
        'S_tail': np.array([e['S_tail_frac'] for e in entries]),
        'S_total': np.array([e['S_total'] for e in entries]),
        'ent_gap': np.array([e['ent_gap'] for e in entries]),
    }

print("Sprint 089b: Level redistribution analysis")
print("=" * 70, flush=True)

# Complex CFT Re(c)
def Re_c(q):
    if q <= 4:
        return {2: 0.5, 3: 0.8, 4: 1.0}[q]
    alpha = np.arccosh(np.sqrt(q) / 2)
    return 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)

# 1. Weight fractions vs n for all q
print("\n--- Weight fractions (probability) ---")
for q_val in [2, 3, 5, 7]:
    d = data[q_val]
    print(f"\nq={q_val} (Re(c)={Re_c(q_val):.3f}):")
    print(f"{'n':>4} {'w_max':>8} {'w_mult':>8} {'w_tail':>10} {'S_l0':>8} {'S_l1':>8} {'S_tail':>8} {'S_tot':>8}")
    for i in range(len(d['ns'])):
        print(f"{int(d['ns'][i]):4d} {d['w_max'][i]:8.5f} {d['w_mult'][i]:8.5f} {d['w_tail'][i]:10.6f} "
              f"{d['S_lev0'][i]:8.4f} {d['S_lev1'][i]:8.4f} {d['S_tail'][i]:8.4f} {d['S_total'][i]:8.4f}")

# 2. Compute derivatives d(%S)/d(ln n) using consecutive pairs
print("\n\n--- Derivatives d(%S)/d(ln n) from consecutive pairs ---")
print(f"{'q':>3} {'n_pair':>10} {'d(S_l0)':>10} {'d(S_l1)':>10} {'d(S_tail)':>10} {'d(w_max)':>10} {'d(w_mult)':>10} {'d(w_tail)':>10}")

for q_val in [2, 3, 5, 7]:
    d = data[q_val]
    ns = d['ns']
    ln_n = np.log(ns)
    derivs = []
    for i in range(len(ns) - 1):
        dln = ln_n[i+1] - ln_n[i]
        dS0 = (d['S_lev0'][i+1] - d['S_lev0'][i]) / dln
        dS1 = (d['S_lev1'][i+1] - d['S_lev1'][i]) / dln
        dSt = (d['S_tail'][i+1] - d['S_tail'][i]) / dln
        dwm = (d['w_max'][i+1] - d['w_max'][i]) / dln
        dwmult = (d['w_mult'][i+1] - d['w_mult'][i]) / dln
        dwtail = (d['w_tail'][i+1] - d['w_tail'][i]) / dln
        pair = f"{int(ns[i])}-{int(ns[i+1])}"
        print(f"{q_val:3d} {pair:>10} {dS0:10.4f} {dS1:10.4f} {dSt:10.4f} {dwm:10.5f} {dwmult:10.5f} {dwtail:10.6f}")
        derivs.append({
            'n_pair': pair, 'n_lo': int(ns[i]), 'n_hi': int(ns[i+1]),
            'dS_lev0': float(dS0), 'dS_lev1': float(dS1), 'dS_tail': float(dSt),
            'dw_max': float(dwm), 'dw_mult': float(dwmult), 'dw_tail': float(dwtail),
        })
    results['derivatives'][str(q_val)] = derivs

# 3. Average derivatives for n≥10 (asymptotic regime)
print("\n\n--- Asymptotic derivatives (n≥10 pairs) ---")
print(f"{'q':>3} {'Re(c)':>8} {'<d(S_l0)>':>10} {'<d(S_l1)>':>10} {'<d(S_tail)>':>10} {'walking':>12}")

asymp = {}
for q_val in [2, 3, 5, 7]:
    derivs = results['derivatives'][str(q_val)]
    # Filter for pairs where n_lo >= 10
    big = [d for d in derivs if d['n_lo'] >= 10]
    if len(big) == 0:
        # Fall back to largest available pairs
        big = derivs[-2:] if len(derivs) >= 2 else derivs[-1:]

    avg_dS0 = np.mean([d['dS_lev0'] for d in big])
    avg_dS1 = np.mean([d['dS_lev1'] for d in big])
    avg_dSt = np.mean([d['dS_tail'] for d in big])

    walk_type = {2: 'real', 3: 'real', 5: 'walking', 7: 'broken'}[q_val]
    print(f"{q_val:3d} {Re_c(q_val):8.3f} {avg_dS0:10.4f} {avg_dS1:10.4f} {avg_dSt:10.4f} {walk_type:>12}")

    asymp[q_val] = {
        'avg_dS_lev0': float(avg_dS0),
        'avg_dS_lev1': float(avg_dS1),
        'avg_dS_tail': float(avg_dSt),
        'n_pairs': len(big),
        'walking_type': walk_type,
    }

results['asymptotic_derivatives'] = {str(k): v for k, v in asymp.items()}

# 4. Key discriminators
print("\n\n--- Walking discriminators ---")

# 4a. Ratio d(S_lev1)/d(S_tail) — how fast does lev1 feed tail?
print("\nRatio d(S_lev1)/d(S_tail) (how fast multiplet feeds tail):")
for q_val in [2, 3, 5, 7]:
    a = asymp[q_val]
    if abs(a['avg_dS_tail']) > 1e-6:
        ratio = a['avg_dS_lev1'] / a['avg_dS_tail']
    else:
        ratio = float('inf')
    print(f"  q={q_val}: {ratio:.3f} ({a['walking_type']})")
    asymp[q_val]['lev1_tail_ratio'] = float(ratio)

# 4b. %S(lev0) saturation value — extract from largest n
print("\n%S(lev0) saturation value (largest available n):")
for q_val in [2, 3, 5, 7]:
    d = data[q_val]
    S_l0_last = d['S_lev0'][-1]
    n_last = int(d['ns'][-1])
    print(f"  q={q_val}: %S(lev0) = {S_l0_last:.4f} at n={n_last} ({asymp[q_val]['walking_type']})")
    asymp[q_val]['S_lev0_asymp'] = float(S_l0_last)

# 4c. %S(lev1) at largest n
print("\n%S(lev1) at largest available n:")
for q_val in [2, 3, 5, 7]:
    d = data[q_val]
    S_l1_last = d['S_lev1'][-1]
    n_last = int(d['ns'][-1])
    print(f"  q={q_val}: %S(lev1) = {S_l1_last:.4f} at n={n_last} ({asymp[q_val]['walking_type']})")
    asymp[q_val]['S_lev1_asymp'] = float(S_l1_last)

# 4d. NEW: Multiplet entropy TRANSFER RATE = d(%S_lev1)/d(ln n)
# This should be the walking-specific quantity
print("\nMultiplet entropy transfer rate d(%S_lev1)/d(ln n) at large n:")
for q_val in [2, 3, 5, 7]:
    dS1 = asymp[q_val]['avg_dS_lev1']
    print(f"  q={q_val}: {dS1:.4f} ({asymp[q_val]['walking_type']})")

# 4e. Lev0+Lev1 combined (conformal sector) vs tail
print("\nConformal sector %S(lev0+lev1) at largest n:")
for q_val in [2, 3, 5, 7]:
    d = data[q_val]
    conf = d['S_lev0'][-1] + d['S_lev1'][-1]
    n_last = int(d['ns'][-1])
    print(f"  q={q_val}: {conf:.4f} (tail: {d['S_tail'][-1]:.4f}) at n={n_last}")

# 5. Fit %S(lev0) to saturation model: S_l0(n) = S_l0_inf + A/n^alpha
print("\n\n--- %S(lev0) saturation fits: S_l0 = S_inf + A*n^(-alpha) ---")
def sat_model(x, S_inf, A, alpha):
    return S_inf + A * x**(-alpha)

for q_val in [2, 3, 5, 7]:
    d = data[q_val]
    ns = d['ns']
    S0 = d['S_lev0']
    if len(ns) >= 4:
        try:
            popt, pcov = curve_fit(sat_model, ns, S0, p0=[S0[-1], 1.0, 1.0],
                                   bounds=([0, -10, 0], [1, 10, 5]), maxfev=10000)
            S_inf, A, alpha = popt
            predicted = sat_model(ns, *popt)
            ss_res = np.sum((S0 - predicted)**2)
            ss_tot = np.sum((S0 - np.mean(S0))**2)
            R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            print(f"  q={q_val}: S_l0_inf = {S_inf:.4f}, A = {A:.3f}, alpha = {alpha:.2f}, R² = {R2:.4f}")
            asymp[q_val]['S_lev0_inf'] = float(S_inf)
            asymp[q_val]['S_lev0_sat_alpha'] = float(alpha)
        except Exception as e:
            print(f"  q={q_val}: fit failed — {e}")
    else:
        print(f"  q={q_val}: insufficient data points ({len(ns)})")

# 6. Fit %S(lev1) trajectory: S_l1(n) = S_l1_inf + B*n^(-beta)
print("\n--- %S(lev1) saturation fits: S_l1 = S_inf + B*n^(-beta) ---")
for q_val in [2, 3, 5, 7]:
    d = data[q_val]
    ns = d['ns']
    S1 = d['S_lev1']
    if len(ns) >= 4:
        try:
            popt, pcov = curve_fit(sat_model, ns, S1, p0=[S1[-1], 1.0, 1.0],
                                   bounds=([0, -10, 0], [1, 10, 5]), maxfev=10000)
            S_inf, B, beta = popt
            predicted = sat_model(ns, *popt)
            ss_res = np.sum((S1 - predicted)**2)
            ss_tot = np.sum((S1 - np.mean(S1))**2)
            R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            print(f"  q={q_val}: S_l1_inf = {S_inf:.4f}, B = {B:.3f}, beta = {beta:.2f}, R² = {R2:.4f}")
            asymp[q_val]['S_lev1_inf'] = float(S_inf)
            asymp[q_val]['S_lev1_sat_beta'] = float(beta)
        except Exception as e:
            print(f"  q={q_val}: fit failed — {e}")
    else:
        print(f"  q={q_val}: insufficient data points ({len(ns)})")

results['asymptotic_derivatives'] = {str(k): v for k, v in asymp.items()}

# 7. Summary: which quantity best discriminates walking from non-walking?
print("\n\n" + "=" * 70)
print("SUMMARY: Walking discriminators")
print("=" * 70)

# Compile candidates
print(f"\n{'Discriminator':<35} {'q=2':>8} {'q=3':>8} {'q=5':>8} {'q=7':>8} {'spread':>8} {'monotonic':>10}")

candidates = {}

# %S(lev0) at largest n
vals = [asymp[q]['S_lev0_asymp'] for q in [2,3,5,7]]
spread = max(vals) - min(vals)
mono = all(vals[i] >= vals[i+1] for i in range(len(vals)-1))
print(f"{'%S(lev0) at n_max':<35} {vals[0]:8.4f} {vals[1]:8.4f} {vals[2]:8.4f} {vals[3]:8.4f} {spread:8.4f} {'YES' if mono else 'NO':>10}")
candidates['S_lev0_asymp'] = {'vals': vals, 'spread': spread, 'monotonic': mono}

# %S(lev1) at largest n
vals = [asymp[q]['S_lev1_asymp'] for q in [2,3,5,7]]
spread = max(vals) - min(vals)
mono = all(vals[i] <= vals[i+1] for i in range(len(vals)-1))
print(f"{'%S(lev1) at n_max':<35} {vals[0]:8.4f} {vals[1]:8.4f} {vals[2]:8.4f} {vals[3]:8.4f} {spread:8.4f} {'YES' if mono else 'NO':>10}")
candidates['S_lev1_asymp'] = {'vals': vals, 'spread': spread, 'monotonic': mono}

# d(%S_lev1)/d(ln n)
vals = [asymp[q]['avg_dS_lev1'] for q in [2,3,5,7]]
spread = max(vals) - min(vals)
mono = all(vals[i] >= vals[i+1] for i in range(len(vals)-1))
print(f"{'d(%S_lev1)/d(ln n)':<35} {vals[0]:8.4f} {vals[1]:8.4f} {vals[2]:8.4f} {vals[3]:8.4f} {spread:8.4f} {'YES' if mono else 'NO':>10}")
candidates['dS_lev1'] = {'vals': vals, 'spread': spread, 'monotonic': mono}

# d(%S_tail)/d(ln n)
vals = [asymp[q]['avg_dS_tail'] for q in [2,3,5,7]]
spread = max(vals) - min(vals)
mono = all(vals[i] <= vals[i+1] for i in range(len(vals)-1))
print(f"{'d(%S_tail)/d(ln n)':<35} {vals[0]:8.4f} {vals[1]:8.4f} {vals[2]:8.4f} {vals[3]:8.4f} {spread:8.4f} {'YES' if mono else 'NO':>10}")
candidates['dS_tail'] = {'vals': vals, 'spread': spread, 'monotonic': mono}

# lev1/tail ratio
vals = [asymp[q].get('lev1_tail_ratio', 0) for q in [2,3,5,7]]
spread = max(vals) - min(vals) if all(np.isfinite(v) for v in vals) else float('inf')
mono = all(vals[i] >= vals[i+1] for i in range(len(vals)-1)) if all(np.isfinite(v) for v in vals) else False
print(f"{'lev1/tail transfer ratio':<35} {vals[0]:8.3f} {vals[1]:8.3f} {vals[2]:8.3f} {vals[3]:8.3f} {'inf' if not np.isfinite(spread) else f'{spread:8.3f}':>8} {'YES' if mono else 'NO':>10}")

# S_lev0_inf (saturation value)
vals = [asymp[q].get('S_lev0_inf', asymp[q]['S_lev0_asymp']) for q in [2,3,5,7]]
spread = max(vals) - min(vals)
mono = all(vals[i] >= vals[i+1] for i in range(len(vals)-1))
print(f"{'S_lev0 saturation value':<35} {vals[0]:8.4f} {vals[1]:8.4f} {vals[2]:8.4f} {vals[3]:8.4f} {spread:8.4f} {'YES' if mono else 'NO':>10}")

results['discriminators'] = candidates
save()
print(f"\nSaved to results/sprint_089b_level_redistribution.json")

from db_utils import record
for q_val in [2, 3, 5, 7]:
    a = asymp[q_val]
    record(sprint=89, model='sq_potts', q=q_val, n=0,
           quantity='dS_lev1_dlnn', value=a['avg_dS_lev1'],
           method='dmrg_spectrum_derivative',
           notes=f'multiplet entropy transfer rate, asymptotic')
    record(sprint=89, model='sq_potts', q=q_val, n=0,
           quantity='S_lev0_asymp', value=a['S_lev0_asymp'],
           method='dmrg_spectrum',
           notes=f'ground state entropy fraction at largest n')
    if 'S_lev0_inf' in a:
        record(sprint=89, model='sq_potts', q=q_val, n=0,
               quantity='S_lev0_inf', value=a['S_lev0_inf'],
               method='saturation_fit',
               notes=f'extrapolated S_lev0 saturation value')
print("Recorded to DB.")
