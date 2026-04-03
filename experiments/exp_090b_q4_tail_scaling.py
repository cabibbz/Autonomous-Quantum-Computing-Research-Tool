#!/usr/bin/env python3
"""Sprint 090b: q=4 tail weight power-law analysis.

Test whether q=4 tail weight follows the universal b≈2.0 power law.
Compare to q=2 (b=1.98), q=3 (b=2.01), q=5 (b=2.01), q=7 (b=2.07).
"""
import numpy as np
import json, time
from scipy.optimize import curve_fit

results = {
    'experiment': '090b_q4_tail_scaling',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

def save():
    with open("results/sprint_090b_q4_tail_scaling.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# Load all data including new q=4
files = {
    2: "results/sprint_088b_entspec_dmrg_q2.json",
    3: "results/sprint_088a_entspec_dmrg_q3.json",
    4: "results/sprint_090a_entspec_dmrg_q4.json",
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
    seen = {}
    for entry in all_entries:
        seen[entry['n']] = entry
    entries = sorted(seen.values(), key=lambda x: x['n'])
    data[q_val] = entries

def log_fit(x, a, b):
    return a + b * x

def Re_c(q):
    if q <= 4:
        return {2: 0.5, 3: 0.8, 4: 1.0}[q]
    alpha = np.arccosh(np.sqrt(q) / 2)
    return 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)

print("Sprint 090b: q=4 tail weight scaling analysis")
print("=" * 70, flush=True)

# 1. Power-law fit for ALL q including q=4 (n≥8)
print("\n--- Tail weight power law: w_tail ~ n^b (n≥8) ---")
print(f"{'q':>3} {'n_range':>10} {'npts':>5} {'b':>8} {'b_err':>8} {'R2':>8} {'type':>10}")

all_b = {}
for q_val in [2, 3, 4, 5, 7]:
    entries = [e for e in data[q_val] if e['n'] >= 8]
    ns = np.array([e['n'] for e in entries])
    wt = np.array([e['w_tail'] for e in entries])
    mask = wt > 0
    ln_n = np.log(ns[mask])
    ln_wt = np.log(wt[mask])

    popt, pcov = curve_fit(log_fit, ln_n, ln_wt)
    b = popt[1]
    b_err = np.sqrt(pcov[1, 1])
    predicted = log_fit(ln_n, *popt)
    ss_res = np.sum((ln_wt - predicted)**2)
    ss_tot = np.sum((ln_wt - np.mean(ln_wt))**2)
    R2 = 1 - ss_res / ss_tot

    cft_type = {2: 'real', 3: 'real', 4: 'boundary', 5: 'walking', 7: 'broken'}[q_val]
    n_range = f"{int(ns[0])}-{int(ns[-1])}"
    print(f"{q_val:3d} {n_range:>10} {len(ln_n):5d} {b:8.3f} {b_err:8.3f} {R2:8.4f} {cft_type:>10}")
    all_b[q_val] = {'b': float(b), 'b_err': float(b_err), 'R2': float(R2),
                    'n_range': n_range, 'npts': int(len(ln_n))}

bs = [all_b[q]['b'] for q in [2, 3, 4, 5, 7]]
print(f"\nMean exponent (all q): {np.mean(bs):.3f} ± {np.std(bs):.3f}")
print(f"Max deviation from 2.0: {max(abs(b - 2.0) for b in bs):.3f}")

results['power_law'] = all_b
results['power_law']['mean_b'] = float(np.mean(bs))
results['power_law']['std_b'] = float(np.std(bs))

# 2. Entanglement gap scaling
print("\n\n--- Entanglement gap: Δξ vs ln(n) ---")
for q_val in [2, 3, 4, 5, 7]:
    entries = [e for e in data[q_val] if e['n'] >= 8]
    ns = np.array([e['n'] for e in entries])
    gaps = np.array([e['ent_gap'] for e in entries])
    popt, pcov = curve_fit(log_fit, np.log(ns), gaps)
    print(f"  q={q_val}: Δξ = {popt[0]:.3f} + {popt[1]:.3f}·ln(n)  (slope={popt[1]:.3f})")

# 3. %S(tail) comparison at matched sizes
print("\n\n--- %S(tail) at matched sizes ---")
print(f"{'q':>3} {'n=8':>8} {'n=12':>8} {'n=16':>8} {'n=20':>8} {'n=24':>8}")
for q_val in [2, 3, 4, 5, 7]:
    row = f"{q_val:3d}"
    for n_target in [8, 12, 16, 20, 24]:
        match = [e for e in data[q_val] if e['n'] == n_target]
        if match:
            row += f" {match[0]['S_tail_frac']:8.4f}"
        else:
            row += "       —"
    print(row)

# 4. c_eff convergence for q=4
print("\n\n--- c_eff convergence for q=4 ---")
dd = data[4]
c_eff_pairs = []
for i in range(1, len(dd)):
    n1, S1 = dd[i-1]['n'], dd[i-1]['S_total']
    n2, S2 = dd[i]['n'], dd[i]['S_total']
    c_eff = 6 * (S2 - S1) / (np.log(n2) - np.log(n1))
    ratio = c_eff / Re_c(4)
    c_eff_pairs.append({'n1': n1, 'n2': n2, 'c_eff': float(c_eff), 'ratio': float(ratio)})
    print(f"  ({n1},{n2}): c_eff = {c_eff:.4f}, c_eff/Re(c) = {ratio:.4f}")

# Drift rate dc/d(ln n) for q=4
if len(c_eff_pairs) >= 2:
    c_vals = [p['c_eff'] for p in c_eff_pairs]
    n_mids = [np.sqrt(p['n1'] * p['n2']) for p in c_eff_pairs]
    dc_dlnn = (c_vals[-1] - c_vals[0]) / (np.log(n_mids[-1]) - np.log(n_mids[0]))
    print(f"\n  Drift rate dc/d(ln n) ≈ {dc_dlnn:.4f}")
    print(f"  Compare: q=5 (+0.014), q=6 (-0.048), q=7 (-0.091)")
    results['c_eff'] = {
        'pairs': c_eff_pairs,
        'drift_rate': float(dc_dlnn),
    }

# 5. M/[(q-1)/q] at n=16 for complete comparison
print("\n\n--- M/[(q-1)/q] at n=16 (all q) ---")
print(f"{'q':>3} {'M':>8} {'(q-1)/q':>8} {'ratio':>8} {'type':>10}")
M_at_16 = {}
for q_val in [2, 3, 4, 5, 7]:
    match = [e for e in data[q_val] if e['n'] == 16]
    if match:
        e = match[0]
        S0 = e['S_lev0_frac']
        S1 = e['S_lev1_frac']
        M = S1 / (S0 + S1)
        qm1_q = (q_val - 1) / q_val
        ratio = M / qm1_q
        cft_type = {2: 'real', 3: 'real', 4: 'boundary', 5: 'walking', 7: 'broken'}[q_val]
        print(f"{q_val:3d} {M:8.4f} {qm1_q:8.4f} {ratio:8.4f} {cft_type:>10}")
        M_at_16[q_val] = {'M': float(M), 'ratio': float(ratio)}

results['M_at_16'] = {str(k): v for k, v in M_at_16.items()}
save()

# 6. Interpolation: where exactly does ratio cross 1.0?
print("\n\n--- Interpolation: M/[(q-1)/q] = 1.0 crossing ---")
qs = sorted(M_at_16.keys())
ratios = [M_at_16[q]['ratio'] for q in qs]
# Find where ratio crosses 1.0
for i in range(len(qs) - 1):
    if (ratios[i] - 1.0) * (ratios[i+1] - 1.0) < 0:
        # Linear interpolation
        q_cross = qs[i] + (1.0 - ratios[i]) / (ratios[i+1] - ratios[i]) * (qs[i+1] - qs[i])
        print(f"  Crosses 1.0 between q={qs[i]} and q={qs[i+1]}: q_cross ≈ {q_cross:.2f}")
        results['q_cross'] = float(q_cross)

save()
print(f"\nSaved to results/sprint_090b_q4_tail_scaling.json")

from db_utils import record
record(sprint=90, model='sq_potts', q=4, n=0,
       quantity='tail_exponent_b', value=all_b[4]['b'],
       error=all_b[4]['b_err'],
       method='log_log_fit_n_ge_8',
       notes='w_tail ~ n^b power law exponent')
record(sprint=90, model='sq_potts', q=4, n=0,
       quantity='tail_exponent_b_mean', value=float(np.mean(bs)),
       error=float(np.std(bs)),
       method='cross_q_mean',
       notes='mean b across q=2,3,4,5,7')
print("Recorded to DB.")
