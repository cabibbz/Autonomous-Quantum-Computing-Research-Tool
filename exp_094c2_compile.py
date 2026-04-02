#!/usr/bin/env python3
"""Sprint 094c2: Compile all BW R² vs nA data, fit power laws, extract gamma(q).

q=5 nA=5 was infeasible (24GB memory for H build). Proceed with nA=3,4.
"""
import numpy as np
from db_utils import record
import json, time

t0 = time.time()

results = {
    'experiment': '094c2_compile',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

# Collect all data
all_points = {
    2: [],  # (nA, 1-R², alpha)
    3: [],
    4: [],
    5: [],
}

# 094a: q=2
with open("results/sprint_094a_bw_nA_q2.json") as f:
    d = json.load(f)
for k, v in d['data'].items():
    all_points[2].append((v['nA'], v['one_minus_R2'], v['alpha']))

# 094b: q=3, q=4
with open("results/sprint_094b_bw_nA_q34.json") as f:
    d = json.load(f)
for k, v in d['data'].items():
    if 'R2' not in v:
        continue
    q = v['q']
    all_points[q].append((v['nA'], v['one_minus_R2'], v['alpha']))

# 094c: q=5
with open("results/sprint_094c_bw_nA_q5_compile.json") as f:
    d = json.load(f)
for k, v in d['data'].items():
    if 'R2' not in v:
        continue
    all_points[5].append((v['nA'], v['one_minus_R2'], v['alpha']))

# Print all data
print("Complete BW R² vs nA data:")
print(f"{'q':>3} {'nA':>4} {'1-R²':>12} {'R²':>10} {'alpha':>8}")
print("-" * 42)
for q in sorted(all_points.keys()):
    pts = sorted(all_points[q])
    for nA, omr, alpha in pts:
        print(f"{q:>3} {nA:>4} {omr:>12.4e} {1-omr:>10.6f} {alpha:>8.3f}")
    print()

# Power-law fits for each q
print("\n=== Power-law fits: 1-R² = A * nA^gamma ===")
compilation = {}

for q in sorted(all_points.keys()):
    pts = sorted(all_points[q])
    nA_arr = np.array([p[0] for p in pts], dtype=float)
    omr_arr = np.array([p[1] for p in pts])

    # Log-log fit
    log_nA = np.log(nA_arr)
    log_omr = np.log(omr_arr)
    coeffs = np.polyfit(log_nA, log_omr, 1)
    gamma = coeffs[0]
    A = np.exp(coeffs[1])

    # Fit quality
    pred = A * nA_arr**gamma
    ss_res = np.sum((omr_arr - pred)**2)
    ss_tot = np.sum((omr_arr - np.mean(omr_arr))**2)
    fit_R2 = 1 - ss_res / ss_tot if ss_tot > 1e-30 else 1.0

    print(f"\nq={q}: gamma = {gamma:.2f}, A = {A:.2e}, fit R² = {fit_R2:.4f} ({len(pts)} points)")
    for nA, omr, _ in pts:
        p = A * nA**gamma
        print(f"  nA={nA}: measured {omr:.4e}, predicted {p:.4e}, ratio {omr/p:.3f}")

    compilation[q] = {
        'gamma': float(gamma), 'A': float(A), 'fit_R2': float(fit_R2),
        'n_points': len(pts),
        'nA_values': [int(x) for x in nA_arr],
        'one_minus_R2': [float(x) for x in omr_arr],
    }
    record(sprint=94, model='sq_potts', q=q, n=0,
           quantity='bw_gamma_nA', value=gamma, method='094c_powerlaw_fit')

# Summary table
print("\n\n=== GAMMA(q) SUMMARY ===")
print(f"{'q':>3} {'gamma':>8} {'A':>12} {'fit_R²':>8} {'nA range':>12}")
print("-" * 48)
for q in sorted(compilation.keys()):
    c = compilation[q]
    nA_range = f"{min(c['nA_values'])}-{max(c['nA_values'])}"
    print(f"{q:>3} {c['gamma']:>8.2f} {c['A']:>12.2e} {c['fit_R2']:>8.4f} {nA_range:>12}")

# Ratio analysis: 1-R² at fixed nA=3 and nA=4
print("\n\n=== Walking amplification: (1-R²)_q / (1-R²)_q=2 at fixed nA ===")
for fixed_nA in [3, 4]:
    print(f"\nnA={fixed_nA}:")
    vals = {}
    for q in sorted(all_points.keys()):
        for nA, omr, _ in all_points[q]:
            if nA == fixed_nA:
                vals[q] = omr
    if 2 in vals:
        base = vals[2]
        for q in sorted(vals.keys()):
            print(f"  q={q}: 1-R² = {vals[q]:.4e}, ratio to q=2: {vals[q]/base:.1f}×")

# nA=3 to nA=4 amplification ratio by q
print("\n\n=== Amplification from nA=3 to nA=4 ===")
for q in sorted(all_points.keys()):
    vals = {}
    for nA, omr, _ in all_points[q]:
        vals[nA] = omr
    if 3 in vals and 4 in vals:
        ratio = vals[4] / vals[3]
        print(f"  q={q}: (1-R²)(nA=4)/(1-R²)(nA=3) = {ratio:.1f}×")

# Alpha (BW scale) vs nA
print("\n\n=== BW scale alpha vs nA ===")
for q in sorted(all_points.keys()):
    pts = sorted(all_points[q])
    alphas = [(p[0], p[2]) for p in pts]
    print(f"  q={q}: {[(nA, f'{a:.2f}') for nA, a in alphas]}")

results['compilation'] = {f'q{q}': v for q, v in compilation.items()}
results['all_data'] = {f'q{q}': [(nA, omr, alpha) for nA, omr, alpha in pts]
                       for q, pts in all_points.items()}

with open("results/sprint_094c2_compile.json", "w") as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done!")
