#!/usr/bin/env python3
"""Sprint 084c: Quantitative analysis of entanglement spectrum across walking boundary.

Combine data from 084a and 084b. Analyze:
1. Entanglement gap Δξ scaling with q and n
2. Entropy decomposition: which spectral levels carry the walking signature
3. Normalized spectrum comparison across q (rescale by Δξ)
4. Spectral participation ratio as walking discriminator
"""
import numpy as np
import json, time

results = {
    'experiment': '084c_entspec_analysis',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'analysis': {},
}

def save():
    with open("results/sprint_084c_entspec_analysis.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# Load data from 084a and 084b
with open("results/sprint_084a_entspec_q235.json") as f:
    data_a = json.load(f)
with open("results/sprint_084b_entspec_q678.json") as f:
    data_b = json.load(f)

# Known values
Re_c = {2: 0.500, 3: 0.800, 5: 1.138, 6: 1.253, 7: 1.351, 8: 1.438}
# c_eff from prior sprints (best values, from KNOWLEDGE.md)
c_eff_best = {2: 0.500, 3: 0.893, 5: 1.152, 6: 1.115, 7: 1.059, 8: 1.062}

print("Sprint 084c: Entanglement spectrum analysis across walking boundary")
print("=" * 70)

# 1. Collect comparable data (use largest available n for each q)
print("\n1. ENTANGLEMENT GAP vs q")
print("-" * 70)
print(f"{'q':>3} {'n':>3} {'nA':>3} {'Δξ':>8} {'λ_max':>8} {'S':>8} {'n_sig':>6} {'deg_1':>6}")

all_data = {}
for q_str, qd in list(data_a['data'].items()) + list(data_b['data'].items()):
    q = int(q_str)
    # Get largest n
    max_n = max(qd['sizes'].keys(), key=lambda x: int(x))
    d = qd['sizes'][max_n]
    all_data[q] = d

    # Detect first excited degeneracy
    xi = d['top_ent_energies']
    deg1 = 1
    if len(xi) > 1:
        for j in range(1, min(10, len(xi))):
            if abs(xi[j] - xi[1]) < 0.01 * max(abs(xi[1]), 0.1):
                deg1 += 1
            else:
                break

    d['deg_1'] = deg1
    print(f"{q:3d} {d['n']:3d} {d['nA']:3d} {d['ent_gap']:8.4f} {d['lambda_max']:8.5f} "
          f"{d['S_total']:8.4f} {d['n_significant']:6d} {deg1:6d}")

# 2. Entanglement gap scaling
print("\n\n2. ENTANGLEMENT GAP SCALING")
print("-" * 70)

q_vals = sorted(all_data.keys())
gaps = [all_data[q]['ent_gap'] for q in q_vals]

# Fit Δξ vs q
from numpy.polynomial import polynomial as P
q_arr = np.array(q_vals, dtype=float)
gap_arr = np.array(gaps)

# Linear fit: Δξ = a + b*q
coeffs = np.polyfit(q_arr, gap_arr, 1)
gap_pred = np.polyval(coeffs, q_arr)
ss_res = np.sum((gap_arr - gap_pred)**2)
ss_tot = np.sum((gap_arr - np.mean(gap_arr))**2)
R2_lin = 1 - ss_res/ss_tot

# Log fit: Δξ = a + b*ln(q)
lnq = np.log(q_arr)
coeffs_ln = np.polyfit(lnq, gap_arr, 1)
gap_pred_ln = np.polyval(coeffs_ln, lnq)
ss_res_ln = np.sum((gap_arr - gap_pred_ln)**2)
R2_ln = 1 - ss_res_ln/ss_tot

print(f"Linear fit: Δξ = {coeffs[1]:.4f} + {coeffs[0]:.4f}*q, R²={R2_lin:.6f}")
print(f"Log fit:    Δξ = {coeffs_ln[1]:.4f} + {coeffs_ln[0]:.4f}*ln(q), R²={R2_ln:.6f}")

results['analysis']['gap_scaling'] = {
    'q': q_vals, 'gaps': gaps,
    'linear_fit': {'a': float(coeffs[1]), 'b': float(coeffs[0]), 'R2': float(R2_lin)},
    'log_fit': {'a': float(coeffs_ln[1]), 'b': float(coeffs_ln[0]), 'R2': float(R2_ln)},
}

# 3. Spectral participation ratio (how many eigenvalues effectively contribute)
print("\n\n3. SPECTRAL PARTICIPATION RATIO (IPR)")
print("-" * 70)
print("  PR = 1/Σ(λ_i²) — effective number of participating eigenvalues")
print("  PR/dimA — fraction of Hilbert space participating")
print(f"{'q':>3} {'n':>3} {'PR':>8} {'PR/dimA':>8} {'dimA':>6} {'Re(c)':>7} {'c_eff/Re(c)':>12}")

pr_data = {}
for q in q_vals:
    d = all_data[q]
    lambdas = np.array(d['top_eigenvalues'])
    # PR = 1/Σλ² (only using top eigenvalues, which dominate)
    pr = 1.0 / np.sum(lambdas**2)
    pr_frac = pr / d['dimA']
    ratio = c_eff_best.get(q, 0) / Re_c[q]
    print(f"{q:3d} {d['n']:3d} {pr:8.2f} {pr_frac:8.5f} {d['dimA']:6d} {Re_c[q]:7.3f} {ratio:12.4f}")
    pr_data[q] = {'PR': float(pr), 'PR_frac': float(pr_frac)}

results['analysis']['participation_ratio'] = pr_data

# 4. Entropy decomposition by spectral level
print("\n\n4. ENTROPY DECOMPOSITION BY LEVEL")
print("-" * 70)
print("  S = -Σ λ_i ln(λ_i). Group by degeneracy level.")
print("  Level 0: ground state (non-degenerate)")
print("  Level 1: first excited (q-1 fold)")
print("  Level 2: second excited, etc.")

for q in q_vals:
    d = all_data[q]
    xi = np.array(d['top_ent_energies'])
    lam = np.array(d['top_eigenvalues'])
    s_i = np.array(d['top_entropy_contributions'])
    S_total = d['S_total']

    # Group into levels by degeneracy
    levels = []
    i = 0
    while i < len(xi) and len(levels) < 6:
        group_lam = [lam[i]]
        group_s = [s_i[i]]
        group_xi = [xi[i]]
        j = i + 1
        while j < len(xi) and abs(xi[j] - xi[i]) < 0.02 * max(abs(xi[i]), 0.1):
            group_lam.append(lam[j])
            group_s.append(s_i[j])
            group_xi.append(xi[j])
            j += 1
        levels.append({
            'deg': len(group_lam),
            'xi_mean': float(np.mean(group_xi)),
            'lam_sum': float(np.sum(group_lam)),
            's_sum': float(np.sum(group_s)),
            's_frac': float(np.sum(group_s) / S_total),
        })
        i = j

    print(f"\n  q={q} (n={d['n']}, S={S_total:.4f}, c_eff/Re(c)={c_eff_best.get(q,0)/Re_c[q]:.3f}):")
    print(f"  {'Level':>6} {'deg':>4} {'ξ_mean':>8} {'Σλ':>8} {'ΣS':>8} {'%S':>6}")
    cum_s = 0
    for k, lev in enumerate(levels):
        cum_s += lev['s_frac']
        print(f"  {k:6d} {lev['deg']:4d} {lev['xi_mean']:8.3f} {lev['lam_sum']:8.5f} "
              f"{lev['s_sum']:8.5f} {100*lev['s_frac']:5.1f}%  (cum: {100*cum_s:.1f}%)")

    results['analysis'][f'levels_q{q}'] = levels

# 5. Normalized spectrum comparison
print("\n\n5. NORMALIZED ENTANGLEMENT ENERGIES (ξ/Δξ)")
print("-" * 70)
print("  Rescale all ξ by the gap Δξ = ξ₁ - ξ₀. Do the rescaled spectra match across q?")
print(f"{'q':>3} {'ξ₀/Δξ':>7} {'ξ₁/Δξ':>7} {'ξ₂/Δξ':>7} {'ξ₃/Δξ':>7} {'Δξ':>7} {'deg_1':>6}")

norm_data = {}
for q in q_vals:
    d = all_data[q]
    xi = np.array(d['top_ent_energies'])
    dxi = d['ent_gap']
    xi_norm = (xi - xi[0]) / dxi  # Normalize: ground state = 0, gap = 1

    # Find index of first non-degenerate-with-1 level
    # Level 0 is at xi_norm = 0
    # Level 1 is at xi_norm = 1 (with deg q-1)
    # What's level 2?
    deg1 = d.get('deg_1', 1)
    idx_lev2 = 1 + deg1  # first index beyond level 1
    xi_lev2 = xi_norm[idx_lev2] if idx_lev2 < len(xi_norm) else float('nan')

    # Ratio of level 2 to level 1 gap
    R21 = xi_lev2  # since level 1 is at 1.0 by construction

    print(f"{q:3d} {xi_norm[0]:7.3f} {xi_norm[1]:7.3f} "
          f"{xi_norm[min(2,len(xi_norm)-1)]:7.3f} {xi_norm[min(3,len(xi_norm)-1)]:7.3f} "
          f"{dxi:7.4f} {deg1:6d}  R₂₁={R21:.3f}")

    norm_data[q] = {
        'xi_normalized': [float(x) for x in xi_norm[:15]],
        'R21': float(R21),
    }

results['analysis']['normalized_spectrum'] = norm_data

# 6. N-scaling of entanglement gap
print("\n\n6. ENTANGLEMENT GAP SCALING WITH N (multiple sizes)")
print("-" * 70)

for q_str, qd in list(data_a['data'].items()) + list(data_b['data'].items()):
    q = int(q_str)
    sizes_sorted = sorted(qd['sizes'].keys(), key=int)
    if len(sizes_sorted) >= 2:
        gaps_n = [(int(ns), qd['sizes'][ns]['ent_gap']) for ns in sizes_sorted]
        # Fit Δξ = a + b/N
        N_arr = np.array([g[0] for g in gaps_n], dtype=float)
        gap_arr = np.array([g[1] for g in gaps_n])

        if len(N_arr) >= 2:
            # Linear fit in 1/N
            inv_N = 1.0 / N_arr
            A = np.vstack([inv_N, np.ones_like(inv_N)]).T
            (b, a), _, _, _ = np.linalg.lstsq(A, gap_arr, rcond=None)

            print(f"  q={q}: ", end='')
            for n, g in gaps_n:
                print(f"n={n}→Δξ={g:.4f}  ", end='')
            print(f"  Extrapolated Δξ(∞) = {a:.4f}, slope = {b:.4f}")

# 7. Key discriminator: entropy concentration ratio
print("\n\n7. WALKING DISCRIMINATOR: S(level≥2)/S_total")
print("-" * 70)
print("  If walking breaks in entropy but not gap/correlator, the excess entropy")
print("  must live in the TAIL of the spectrum (levels ≥ 2).")
print(f"{'q':>3} {'S_lev0':>8} {'S_lev1':>8} {'S_tail':>8} {'%tail':>7} {'c/Re(c)':>8}")

tail_data = {}
for q in q_vals:
    levels = results['analysis'].get(f'levels_q{q}', [])
    if len(levels) >= 3:
        s0 = levels[0]['s_frac']
        s1 = levels[1]['s_frac']
        s_tail = 1.0 - s0 - s1
        ratio = c_eff_best.get(q, 0) / Re_c[q]
        print(f"{q:3d} {100*s0:7.1f}% {100*s1:7.1f}% {100*s_tail:7.1f}% {100*s_tail:6.1f}% {ratio:8.3f}")
        tail_data[q] = {'s0': s0, 's1': s1, 's_tail': s_tail}

results['analysis']['tail_entropy'] = tail_data
save()

# Correlation between tail entropy and walking deviation
print("\n\n8. CORRELATION: tail entropy vs c_eff/Re(c)")
print("-" * 70)

q_list = sorted(tail_data.keys())
tails = [tail_data[q]['s_tail'] for q in q_list]
ratios = [c_eff_best[q]/Re_c[q] for q in q_list]

# Pearson correlation
from numpy import corrcoef
r = corrcoef(tails, ratios)[0, 1]
print(f"  Pearson r = {r:.4f}")
print(f"  Interpretation: {'STRONG' if abs(r) > 0.8 else 'MODERATE' if abs(r) > 0.5 else 'WEAK'} "
      f"{'negative' if r < 0 else 'positive'} correlation")
print(f"  More tail entropy ↔ {'LOWER' if r < 0 else 'HIGHER'} c_eff/Re(c)")

results['analysis']['correlation'] = {
    'pearson_r': float(r),
    'q_values': q_list,
    'tail_entropy': tails,
    'c_ratio': ratios,
}

save()
print(f"\nSaved to results/sprint_084c_entspec_analysis.json")

from db_utils import record
for q in q_vals:
    d = all_data[q]
    record(sprint=84, model='sq_potts', q=q, n=d['n'],
           quantity='ent_gap', value=d['ent_gap'],
           method='exact_diag_periodic',
           notes=f'deg1={d.get("deg_1",0)}, PR={pr_data[q]["PR"]:.2f}')
print("Recorded to DB.")
