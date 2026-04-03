#!/usr/bin/env python3
"""Sprint 098b: Pairwise convergence extrapolation and c_Casimir scope analysis.

Questions:
1. Does pairwise c_implied(N1,N2) converge to Re(c) as N grows?
2. Can we extrapolate c_implied(∞) from the pairwise trend?
3. At what N does 2-param Casimir deviate >5% from Re(c)? (finite-size scope)
4. How independent are c_Casimir and c_entropy? (use same velocity v from gap)

Uses data from Sprint 098a results.
"""
import numpy as np
import json

# Load 098a results
with open("results/sprint_098a_casimir_gpu_extended.json") as f:
    data = json.load(f)

rec_values = {2: 0.500, 3: 0.800, 5: 1.138, 7: 1.351, 8: 1.438}
c_eff_values = {2: 0.500, 3: 0.893, 5: 1.152, 7: 1.059, 8: 1.062}

print("Sprint 098b: Pairwise Convergence and Scope Analysis")
print("=" * 70)

# 1. Pairwise convergence extrapolation
print("\n1. PAIRWISE c_implied/Re(c) CONVERGENCE")
print("-" * 70)

extrapolated = {}
for q_str in ['2', '3', '5', '7', '8']:
    q = int(q_str)
    d = data['data'][q_str]
    pairs = d['pairwise']
    rec = rec_values[q]

    # Extract midpoint N and c_implied/Re(c) for each pair
    sizes = d['sizes']
    N_mid = []
    c_ratios = []
    for i, p in enumerate(pairs):
        N1, N2 = sizes[i], sizes[i+1]
        N_mid.append((N1 + N2) / 2)
        c_ratios.append(p['c_imp_over_rec'])

    N_mid = np.array(N_mid)
    c_ratios = np.array(c_ratios)

    print(f"\n  q={q} (Re(c)={rec}):")
    for i, p in enumerate(pairs):
        print(f"    N_mid={N_mid[i]:.1f}: c/Re(c)={c_ratios[i]:.4f}")

    # Try extrapolation: c/Re(c) = a + b/N² (leading FSS correction)
    if len(N_mid) >= 3:
        x = 1.0 / N_mid**2
        A = np.vstack([x, np.ones_like(x)]).T
        (slope, intercept), _, _, _ = np.linalg.lstsq(A, c_ratios, rcond=None)
        extrapolated[q] = intercept
        print(f"    Extrapolation c/Re(c)(∞) = {intercept:.4f} (fit c/Re(c) = {intercept:.4f} + {slope:.2f}/N²)")

        # Also try 1/N extrapolation
        x1 = 1.0 / N_mid
        A1 = np.vstack([x1, np.ones_like(x1)]).T
        (slope1, intercept1), _, _, _ = np.linalg.lstsq(A1, c_ratios, rcond=None)
        print(f"    Alt extrap (1/N): c/Re(c)(∞) = {intercept1:.4f}")
    elif len(N_mid) == 2:
        # Linear extrap from 2 points
        x = 1.0 / N_mid**2
        slope = (c_ratios[1] - c_ratios[0]) / (x[1] - x[0])
        intercept = c_ratios[0] - slope * x[0]
        extrapolated[q] = intercept
        print(f"    Linear extrap c/Re(c)(∞) = {intercept:.4f}")

# 2. Finite-size scope: when does 2-param c_Casimir break?
print(f"\n\n2. FINITE-SIZE SCOPE OF CASIMIR = Re(c)")
print("-" * 70)
print("   2-param fit c_implied/Re(c) vs number of data points:")
print(f"   {'q':>3} {'N_max':>6} {'c/Rec(2p)':>10} {'c/Rec(3p)':>10} {'c/Rec(pw_last)':>14} {'dev_2p%':>8}")
for q_str in ['2', '3', '5', '7', '8']:
    q = int(q_str)
    d = data['data'][q_str]
    f2 = d['fit_2param']
    pw_last = d['pairwise'][-1]['c_imp_over_rec']
    f3_str = f"{d['fit_3param']['c_imp_over_rec']:.4f}" if 'fit_3param' in d else "N/A"
    dev = abs(f2['c_imp_over_rec'] - 1) * 100
    print(f"   {q:3d} {max(d['sizes']):6d} {f2['c_imp_over_rec']:10.4f} {f3_str:>10} {pw_last:14.4f} {dev:8.1f}")

# The key: does 1/N^4 coefficient grow with q?
print(f"\n   1/N⁴ correction strength:")
print(f"   {'q':>3} {'coeff(1/N⁴)':>12} {'|coeff|/q':>10}")
for q_str in ['2', '3', '5', '7', '8']:
    q = int(q_str)
    d = data['data'][q_str]
    if 'fit_3param' in d:
        coeff = d['fit_3param']['coeff_N4']
        print(f"   {q:3d} {coeff:12.4f} {abs(coeff)/q:10.4f}")

# 3. Independence comparison: c_Casimir vs c_entropy
print(f"\n\n3. c_Casimir vs c_entropy — INDEPENDENT COMPARISON")
print("-" * 70)
print("   Both use v = gap*N/(2π*x_σ) from the SAME gap measurement.")
print("   c_Casimir = (vc from E₀ fit) / v")
print("   c_entropy = 6·ΔS / ((1+1/α)·Δln(N))  [from entropy, independent of gap]")
print(f"\n   {'q':>3} {'c_Cas':>7} {'c_eff':>7} {'Re(c)':>7} {'Cas/Rec':>8} {'eff/Rec':>8} {'Cas-eff':>8}")
for q_str in ['2', '3', '5', '7', '8']:
    q = int(q_str)
    d = data['data'][q_str]
    f2 = d['fit_2param']
    rec = rec_values[q]
    c_eff = c_eff_values[q]
    print(f"   {q:3d} {f2['c_implied']:7.3f} {c_eff:7.3f} {rec:7.3f} "
          f"{f2['c_imp_over_rec']:8.4f} {c_eff/rec:8.4f} "
          f"{f2['c_implied']-c_eff:+8.3f}")

# 4. Use 3-param fit as "corrected" estimate
print(f"\n\n4. CORRECTED c_Casimir (3-param fit, removing 1/N⁴)")
print("-" * 70)
print(f"   {'q':>3} {'c_Cas(2p)/Rec':>14} {'c_Cas(3p)/Rec':>14} {'extrap pw/Rec':>14}")
for q_str in ['2', '3', '5', '7', '8']:
    q = int(q_str)
    d = data['data'][q_str]
    f2 = d['fit_2param']
    f3 = d.get('fit_3param', {})
    ext = extrapolated.get(q, None)
    f3_str = f"{f3['c_imp_over_rec']:.4f}" if f3 else "N/A"
    ext_str = f"{ext:.4f}" if ext else "N/A"
    print(f"   {q:3d} {f2['c_imp_over_rec']:14.4f} {f3_str:>14} {ext_str:>14}")

# Summary statistics
print(f"\n\n5. SUMMARY STATISTICS")
print("-" * 70)

# 2-param fit
vals_2p = [data['data'][q]['fit_2param']['c_imp_over_rec'] for q in ['2','3','5','7','8']]
print(f"   2-param c/Re(c): mean={np.mean(vals_2p):.4f}, std={np.std(vals_2p):.4f}, "
      f"range=[{min(vals_2p):.4f}, {max(vals_2p):.4f}]")

# 3-param fit
vals_3p = [data['data'][q]['fit_3param']['c_imp_over_rec']
           for q in ['2','3','5','7','8'] if 'fit_3param' in data['data'][q]]
if vals_3p:
    print(f"   3-param c/Re(c): mean={np.mean(vals_3p):.4f}, std={np.std(vals_3p):.4f}, "
          f"range=[{min(vals_3p):.4f}, {max(vals_3p):.4f}]")

# Pairwise last (largest sizes)
vals_pw = [data['data'][q]['pairwise'][-1]['c_imp_over_rec'] for q in ['2','3','5','7','8']]
print(f"   Pairwise-last c/Re(c): mean={np.mean(vals_pw):.4f}, std={np.std(vals_pw):.4f}, "
      f"range=[{min(vals_pw):.4f}, {max(vals_pw):.4f}]")

# Extrapolated
vals_ext = [extrapolated[q] for q in [2,3,5,7,8] if q in extrapolated]
if vals_ext:
    print(f"   Extrapolated c/Re(c): mean={np.mean(vals_ext):.4f}, std={np.std(vals_ext):.4f}, "
          f"range=[{min(vals_ext):.4f}, {max(vals_ext):.4f}]")

# c_eff comparison
vals_eff = [c_eff_values[q]/rec_values[q] for q in [2,3,5,7,8]]
print(f"   c_eff/Re(c):     mean={np.mean(vals_eff):.4f}, std={np.std(vals_eff):.4f}, "
      f"range=[{min(vals_eff):.4f}, {max(vals_eff):.4f}]")

print(f"\n   ** Casimir spread ({np.std(vals_pw):.4f}) vs entropy spread ({np.std(vals_eff):.4f}) **")
print(f"   ** Casimir is {np.std(vals_eff)/np.std(vals_pw):.1f}× more consistent with Re(c) than entropy **")

# Save analysis
analysis = {
    'experiment': '098b_casimir_convergence',
    'extrapolated_c_over_rec': {str(q): float(v) for q, v in extrapolated.items()},
    'summary': {
        '2param_mean': float(np.mean(vals_2p)),
        '2param_std': float(np.std(vals_2p)),
        'pairwise_last_mean': float(np.mean(vals_pw)),
        'pairwise_last_std': float(np.std(vals_pw)),
        'c_eff_mean': float(np.mean(vals_eff)),
        'c_eff_std': float(np.std(vals_eff)),
        'consistency_ratio': float(np.std(vals_eff)/np.std(vals_pw)),
    }
}
if vals_ext:
    analysis['summary']['extrap_mean'] = float(np.mean(vals_ext))
    analysis['summary']['extrap_std'] = float(np.std(vals_ext))

with open("results/sprint_098b_casimir_convergence.json", "w") as f:
    json.dump(analysis, f, indent=2)
print(f"\nSaved to results/sprint_098b_casimir_convergence.json")

from db_utils import record
record(sprint=98, model='sq_potts', q=0, n=0,
       quantity='casimir_consistency', value=float(np.std(vals_eff)/np.std(vals_pw)),
       method='convergence_analysis',
       notes=f"Casimir {np.std(vals_pw):.4f} spread vs entropy {np.std(vals_eff):.4f} spread")
print("Recorded to DB.")
