#!/usr/bin/env python3
"""Sprint 047d: Verify ν(q=7) with corrected slope analysis.

The crossing at g_c≈0.244 means the earlier slope window [0.24, 0.28]
spans the crossing, invalidating slope ratio. Use a narrower window
on the disordered side only.
"""
import numpy as np
import json

with open('results/sprint_047a_q7_n8.json') as f:
    d8 = json.load(f)
with open('results/sprint_047b_q7_n12.json') as f:
    d12 = json.load(f)

g8 = np.array([pt['g'] for pt in d8['data'] if pt['status']=='ok'])
cv8 = np.array([pt['cv'] for pt in d8['data'] if pt['status']=='ok'])
g12 = np.array([pt['g'] for pt in d12['data'] if pt['status']=='ok'])
cv12 = np.array([pt['cv'] for pt in d12['data'] if pt['status']=='ok'])

# Method 1: Slope on disordered side only (g=0.26 to 0.30)
print("=== Method 1: Slope on disordered side ===")
pairs = [(0.26, 0.30), (0.26, 0.28), (0.28, 0.30)]
for glo, ghi in pairs:
    i8lo = np.argmin(np.abs(g8-glo)); i8hi = np.argmin(np.abs(g8-ghi))
    i12lo = np.argmin(np.abs(g12-glo)); i12hi = np.argmin(np.abs(g12-ghi))
    s8 = (cv8[i8hi]-cv8[i8lo])/(ghi-glo)
    s12 = (cv12[i12hi]-cv12[i12lo])/(ghi-glo)
    r = s12/s8
    if r > 0:
        nu = np.log(12/8)/np.log(r)
        print(f"  [{glo},{ghi}]: s8={s8:.3f}, s12={s12:.3f}, ratio={r:.3f}, ν={nu:.2f}")
    else:
        print(f"  [{glo},{ghi}]: s8={s8:.3f}, s12={s12:.3f}, ratio={r:.3f}, ν=N/A (negative)")

# Method 2: Slope on ordered side (g=0.15 to 0.24)
print("\n=== Method 2: Slope on ordered side ===")
pairs = [(0.15, 0.24), (0.20, 0.24), (0.15, 0.20)]
for glo, ghi in pairs:
    i8lo = np.argmin(np.abs(g8-glo)); i8hi = np.argmin(np.abs(g8-ghi))
    i12lo = np.argmin(np.abs(g12-glo)); i12hi = np.argmin(np.abs(g12-ghi))
    s8 = abs((cv8[i8hi]-cv8[i8lo])/(ghi-glo))
    s12 = abs((cv12[i12hi]-cv12[i12lo])/(ghi-glo))
    r = s12/s8
    if r > 0 and r != 1:
        nu = np.log(12/8)/np.log(r)
        print(f"  [{glo},{ghi}]: |s8|={s8:.3f}, |s12|={s12:.3f}, ratio={r:.3f}, ν={nu:.2f}")
    else:
        print(f"  [{glo},{ghi}]: |s8|={s8:.3f}, |s12|={s12:.3f}, ratio={r:.3f}")

# Method 3: Maximum CV gradient near g_c
print("\n=== Method 3: Maximum CV gradient ===")
# For n=8, max gradient is between g=0.25 and g=0.26 (CV drops from 0.529 to 0.425)
max_grad_8 = abs(cv8[np.argmin(np.abs(g8-0.26))] - cv8[np.argmin(np.abs(g8-0.25))]) / 0.01
# For n=12, max gradient near transition (0.24 to 0.26)
max_grad_12 = abs(cv12[np.argmin(np.abs(g12-0.26))] - cv12[np.argmin(np.abs(g12-0.24))]) / 0.02
print(f"  n=8 max gradient (g=0.25→0.26): {max_grad_8:.2f}")
print(f"  n=12 gradient (g=0.24→0.26): {max_grad_12:.2f}")
ratio = max_grad_12/max_grad_8
print(f"  ratio: {ratio:.4f}")
if ratio > 0 and ratio != 1:
    nu = np.log(12/8)/np.log(ratio)
    print(f"  ν from max gradient ratio: {nu:.2f}")

# Method 4: Data collapse with more careful optimization
print("\n=== Method 4: Data collapse (careful) ===")
from scipy.optimize import minimize

def collapse_quality(params, sizes_g, sizes_cv):
    nu, gc = params
    if nu <= 0.1 or nu > 10: return 1e6
    all_x, all_cv = [], []
    for g_arr, cv_arr, n in sizes_g:
        x = (g_arr - gc) * n**(1/nu)
        all_x.extend(x.tolist())
        all_cv.extend(cv_arr.tolist())
    all_x = np.array(all_x); all_cv = np.array(all_cv)
    idx = np.argsort(all_x)
    all_x, all_cv = all_x[idx], all_cv[idx]

    # Quality: sum of squared CV differences for nearby x-values from different sizes
    total_err = 0; count = 0
    for i in range(len(all_x)-1):
        if abs(all_x[i+1] - all_x[i]) < 0.5:  # nearby in rescaled coord
            total_err += (all_cv[i+1] - all_cv[i])**2
            count += 1
    return total_err / max(count, 1)

sizes_g = [(g8, cv8, 8), (g12, cv12, 12)]

# Grid search
print("  Grid search...")
best_q = 1e10; best_p = (1.0, 0.26)
for nu in np.arange(0.3, 5.0, 0.05):
    for gc in np.arange(0.22, 0.30, 0.002):
        q = collapse_quality([nu, gc], sizes_g, None)
        if q < best_q:
            best_q = q; best_p = (nu, gc)
print(f"  Grid: ν={best_p[0]:.2f}, g_c={best_p[1]:.3f}, quality={best_q:.8f}")

# Refine
res = minimize(collapse_quality, best_p, args=(sizes_g, None),
               method='Nelder-Mead', options={'xatol':0.005, 'fatol':1e-10})
nu_opt, gc_opt = res.x
print(f"  Refined: ν={nu_opt:.3f}, g_c={gc_opt:.4f}, quality={res.fun:.8f}")

# Test known values
print("\n  Quality comparison:")
for nu_test, name in [(0.45, 'optimal~'), (0.5, 'mean-field'), (1.0, 'Ising'),
                       (5/6, 'q=3 Potts'), (2.0, 'q=5'), (1.5, 'test')]:
    q = collapse_quality([nu_test, gc_opt], sizes_g, None)
    print(f"    ν={nu_test:.2f} ({name}): quality={q:.8f}, ratio={q/res.fun:.2f}")

# Summary
print(f"\n=== FINAL SUMMARY: ν(q=7) ===")
print(f"  Data collapse optimal: ν = {nu_opt:.2f}, g_c = {gc_opt:.4f}")
print(f"  Crossing point: g_c ≈ 0.244 (n=8,12)")
print(f"  Qualitative: CROSSING CURVES present → standard 2nd-order (not BKT)")
print(f"  Comparison: q=4 (ν≥2.2, NO crossings) vs q=7 (ν≈{nu_opt:.1f}, crossings)")

results = {
    'sprint': '047d',
    'nu_optimal': float(nu_opt),
    'gc_optimal': float(gc_opt),
    'quality_optimal': float(res.fun),
    'crossing_gc': 0.244,
    'has_crossings': True,
}
with open('results/sprint_047d_nu_verify.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved to results/sprint_047d_nu_verify.json")
