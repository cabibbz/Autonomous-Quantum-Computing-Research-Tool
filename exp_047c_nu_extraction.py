#!/usr/bin/env python3
"""Sprint 047c: Extract ν(q=7) from MI-CV slope ratio and data collapse.

Uses n=8,12 data from experiments 047a,047b.
"""
import numpy as np
import json
from scipy.optimize import minimize_scalar, minimize

# Load data
with open('results/sprint_047a_q7_n8.json') as f:
    d8 = json.load(f)
with open('results/sprint_047b_q7_n12.json') as f:
    d12 = json.load(f)

# Extract CV vs g for both sizes
g8 = np.array([pt['g'] for pt in d8['data'] if pt['status']=='ok'])
cv8 = np.array([pt['cv'] for pt in d8['data'] if pt['status']=='ok'])
g12 = np.array([pt['g'] for pt in d12['data'] if pt['status']=='ok'])
cv12 = np.array([pt['cv'] for pt in d12['data'] if pt['status']=='ok'])

print("=== n=8 vs n=12 comparison ===")
print("  g     n=8    n=12   n12-n8  crossing?")
for g in [0.15, 0.20, 0.24, 0.26, 0.28, 0.30, 0.35]:
    i8 = np.argmin(np.abs(g8 - g))
    i12 = np.argmin(np.abs(g12 - g))
    if abs(g8[i8] - g) < 0.01 and abs(g12[i12] - g) < 0.01:
        diff = cv12[i12] - cv8[i8]
        cross = "n12<n8 (ordered)" if diff < 0 else "n12>n8 (disordered)"
        print(f"  {g:.2f}  {cv8[i8]:.4f}  {cv12[i12]:.4f}  {diff:+.4f}  {cross}")

# Find crossing point by interpolation
# Need common g values where both have data
common_g = []
diff_vals = []
for g in g12:
    i8 = np.argmin(np.abs(g8 - g))
    if abs(g8[i8] - g) < 0.011:
        common_g.append(g)
        i12 = np.where(g12 == g)[0][0]
        diff_vals.append(cv12[i12] - cv8[i8])

common_g = np.array(common_g)
diff_vals = np.array(diff_vals)

print(f"\n=== Crossing point interpolation ===")
for i in range(len(diff_vals)-1):
    if diff_vals[i] * diff_vals[i+1] < 0:  # sign change
        # Linear interpolation
        g_cross = common_g[i] - diff_vals[i] * (common_g[i+1]-common_g[i]) / (diff_vals[i+1]-diff_vals[i])
        print(f"  Crossing between g={common_g[i]:.3f} and g={common_g[i+1]:.3f}")
        print(f"  Interpolated g_c = {g_cross:.4f}")

# Slope at g_c via finite differences
# Use g=0.24 to g=0.28 (spanning g_c=0.259)
g_lo, g_hi = 0.24, 0.28
dg = g_hi - g_lo

i8_lo = np.argmin(np.abs(g8 - g_lo))
i8_hi = np.argmin(np.abs(g8 - g_hi))
i12_lo = np.argmin(np.abs(g12 - g_lo))
i12_hi = np.argmin(np.abs(g12 - g_hi))

slope8 = (cv8[i8_hi] - cv8[i8_lo]) / dg
slope12 = (cv12[i12_hi] - cv12[i12_lo]) / dg

print(f"\n=== Slope at g_c (dCV/dg over [{g_lo}, {g_hi}]) ===")
print(f"  n=8:  slope = {slope8:.4f}")
print(f"  n=12: slope = {slope12:.4f}")
print(f"  Ratio n12/n8 = {slope12/slope8:.4f}")

# ν from slope ratio: slope ~ n^(1/ν)
# slope12/slope8 = (12/8)^(1/ν)
ratio = abs(slope12 / slope8)
nu_slope = np.log(12/8) / np.log(ratio) if ratio > 1 else float('inf')
print(f"\n=== ν from slope ratio ===")
print(f"  slope_ratio = {ratio:.4f}")
print(f"  (12/8)^(1/ν) = {ratio:.4f}")
print(f"  1/ν = ln({ratio:.4f}) / ln(1.5) = {np.log(ratio)/np.log(1.5):.4f}")
print(f"  ν = {nu_slope:.2f}")

# Data collapse: CV(g, n) = F((g - g_c) * n^(1/ν))
# Quality = sum of squared differences between rescaled curves
def collapse_quality(params):
    nu, gc = params
    x8 = (g8 - gc) * 8**(1/nu)
    x12 = (g12 - gc) * 12**(1/nu)
    # Interpolate n=8 onto n=12 x-values
    total_err = 0
    count = 0
    for i, x12i in enumerate(x12):
        # Find two closest n=8 points
        in_range = (x12i >= np.min(x8)) and (x12i <= np.max(x8))
        if in_range:
            cv8_interp = np.interp(x12i, x8, cv8)
            total_err += (cv12[i] - cv8_interp)**2
            count += 1
    return total_err / max(count, 1)

print(f"\n=== Data collapse optimization ===")
# Grid search first
best_q = 1e10
best_params = (1.0, 0.26)
for nu in np.arange(0.5, 4.0, 0.1):
    for gc in np.arange(0.23, 0.30, 0.005):
        q = collapse_quality([nu, gc])
        if q < best_q:
            best_q = q
            best_params = (nu, gc)

print(f"  Grid search best: ν={best_params[0]:.1f}, g_c={best_params[1]:.3f}, quality={best_q:.6f}")

# Refine
res = minimize(collapse_quality, best_params, method='Nelder-Mead',
               options={'xatol': 0.01, 'fatol': 1e-8})
nu_opt, gc_opt = res.x
print(f"  Refined: ν={nu_opt:.2f}, g_c={gc_opt:.4f}, quality={res.fun:.6f}")

# Compare with known ν values
for nu_test, name in [(1.0, 'Ising'), (5/6, 'Potts q=3'), (2.0, 'q=5'), (0.5, 'mean-field')]:
    q = collapse_quality([nu_test, gc_opt])
    print(f"  ν={nu_test:.2f} ({name}): quality={q:.6f}, ratio to optimal={q/res.fun:.2f}")

print(f"\n=== Summary: ν(q=7) ===")
print(f"  From slope ratio (n=8,12): ν = {nu_slope:.2f}")
print(f"  From data collapse: ν = {nu_opt:.2f}, g_c = {gc_opt:.4f}")

# Save results
results = {
    'sprint': '047c',
    'crossing_analysis': {
        'common_g': common_g.tolist(),
        'diff_n12_minus_n8': diff_vals.tolist(),
    },
    'slope_analysis': {
        'g_range': [g_lo, g_hi],
        'slope_n8': float(slope8),
        'slope_n12': float(slope12),
        'ratio': float(ratio),
        'nu_slope': float(nu_slope),
    },
    'collapse_analysis': {
        'nu_optimal': float(nu_opt),
        'gc_optimal': float(gc_opt),
        'quality_optimal': float(res.fun),
    },
}
with open('results/sprint_047c_nu_extraction.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved to results/sprint_047c_nu_extraction.json")
