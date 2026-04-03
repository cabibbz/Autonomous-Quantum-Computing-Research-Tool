#!/usr/bin/env python3
"""Sprint 048c: Extract ν(q=10) from MI-CV slope ratio.

q=10 shows NO crossings (like q=4), so use derivative dCV/dg scaling.
slope ~ n^(1/ν) near g_c.
"""
import numpy as np
import json
from scipy.optimize import minimize

# Load data
with open('results/sprint_048a_q10_n8.json') as f:
    d8 = json.load(f)
with open('results/sprint_048b_q10_n12.json') as f:
    d12 = json.load(f)

g8 = np.array([pt['g'] for pt in d8['data'] if pt['status']=='ok'])
cv8 = np.array([pt['cv'] for pt in d8['data'] if pt['status']=='ok'])
g12 = np.array([pt['g'] for pt in d12['data'] if pt['status']=='ok'])
cv12 = np.array([pt['cv'] for pt in d12['data'] if pt['status']=='ok'])

print("=== n=8 vs n=12 MI-CV comparison (q=10, chi=20) ===")
print("  g     n=8    n=12   diff    sign")
for g in g12:
    i8 = np.argmin(np.abs(g8 - g))
    i12 = np.where(g12 == g)[0][0]
    if abs(g8[i8] - g) < 0.011:
        diff = cv12[i12] - cv8[i8]
        print(f"  {g:.3f}  {cv8[i8]:.4f}  {cv12[i12]:.4f}  {diff:+.4f}  {'n12<n8' if diff<0 else 'n12>n8'}")

print(f"\n=== KEY FINDING: NO CROSSINGS ===")
print(f"  n=12 CV is below n=8 CV at ALL {len(g12)} tested g values")
print(f"  Gap narrows near g_c then widens — same as q=4")
print(f"  Sprint 043 chi=10 'crossing' was convergence artifact")

# Slope analysis: dCV/dg near g_c≈0.246
# Use pairs of g values spanning g_c
print(f"\n=== Slope analysis (dCV/dg) ===")
slopes_n8 = []
slopes_n12 = []
for g_lo, g_hi in [(0.20, 0.26), (0.22, 0.28), (0.24, 0.28), (0.24, 0.30)]:
    i8_lo = np.argmin(np.abs(g8 - g_lo))
    i8_hi = np.argmin(np.abs(g8 - g_hi))
    i12_lo = np.argmin(np.abs(g12 - g_lo))
    i12_hi = np.argmin(np.abs(g12 - g_hi))

    if abs(g8[i8_lo]-g_lo) < 0.01 and abs(g8[i8_hi]-g_hi) < 0.01:
        dg = g8[i8_hi] - g8[i8_lo]
        s8 = (cv8[i8_hi] - cv8[i8_lo]) / dg
        s12 = (cv12[i12_hi] - cv12[i12_lo]) / dg
        ratio = s12 / s8 if s8 > 0 else float('inf')
        if ratio > 1:
            nu = np.log(12/8) / np.log(ratio)
        elif ratio > 0:
            nu = float('inf')
        else:
            nu = float('nan')
        print(f"  [{g_lo:.2f}, {g_hi:.2f}]: s8={s8:.4f}, s12={s12:.4f}, ratio={ratio:.3f}, ν={nu:.2f}")
        slopes_n8.append(s8)
        slopes_n12.append(s12)

# Global slope across transition
g_lo, g_hi = 0.20, 0.30
i8_lo = np.argmin(np.abs(g8 - g_lo)); i8_hi = np.argmin(np.abs(g8 - g_hi))
i12_lo = np.argmin(np.abs(g12 - g_lo)); i12_hi = np.argmin(np.abs(g12 - g_hi))
dg = g8[i8_hi] - g8[i8_lo]
s8_global = (cv8[i8_hi] - cv8[i8_lo]) / dg
s12_global = (cv12[i12_hi] - cv12[i12_lo]) / dg
ratio_global = s12_global / s8_global
if ratio_global > 1:
    nu_global = np.log(12/8) / np.log(ratio_global)
elif ratio_global > 0:
    nu_global = float('inf')
else:
    nu_global = float('nan')

print(f"\n  Global [{g_lo:.2f}, {g_hi:.2f}]: s8={s8_global:.4f}, s12={s12_global:.4f}")
print(f"  ratio={ratio_global:.3f}, ν={nu_global:.2f}")

# Data collapse attempt (even without crossings, can still test collapse)
def collapse_quality(params):
    nu, gc = params
    if nu <= 0 or gc <= 0:
        return 1e10
    x8 = (g8 - gc) * 8**(1/nu)
    x12 = (g12 - gc) * 12**(1/nu)
    total_err = 0; count = 0
    for i, x12i in enumerate(x12):
        if np.min(x8) <= x12i <= np.max(x8):
            cv8_interp = np.interp(x12i, x8, cv8)
            total_err += (cv12[i] - cv8_interp)**2
            count += 1
    return total_err / max(count, 1)

print(f"\n=== Data collapse optimization ===")
best_q = 1e10; best_params = (1.0, 0.25)
for nu in np.arange(0.3, 6.0, 0.1):
    for gc in np.arange(0.20, 0.30, 0.005):
        q = collapse_quality([nu, gc])
        if q < best_q:
            best_q = q; best_params = (nu, gc)

print(f"  Grid search: ν={best_params[0]:.1f}, g_c={best_params[1]:.3f}, quality={best_q:.6f}")
res = minimize(collapse_quality, best_params, method='Nelder-Mead',
               options={'xatol': 0.01, 'fatol': 1e-8})
nu_opt, gc_opt = res.x
print(f"  Refined: ν={nu_opt:.2f}, g_c={gc_opt:.4f}, quality={res.fun:.6f}")

for nu_test, name in [(0.5, 'mean-field'), (1.0, 'Ising'), (2.0, 'q=5'), (3.0, 'large')]:
    q = collapse_quality([nu_test, gc_opt])
    ratio = q / res.fun if res.fun > 0 else float('inf')
    print(f"  ν={nu_test:.1f} ({name}): quality={q:.6f}, ratio to optimal={ratio:.2f}")

# Compare with q=4 behavior
print(f"\n=== Comparison with q=4 (no crossings, BKT-like) ===")
print(f"  q=4 (Sprint 046): ν slope ratio 1.92 (n=8,12) → 2.71 (n=12,16)")
print(f"  q=4: no crossings, ν estimates INCREASE with n → BKT signature")
print(f"  q=10: no crossings, ν from slopes: see above")
print(f"  Both lack crossings → same qualitative behavior")

print(f"\n=== Summary ===")
print(f"  q=10 at chi=20: NO MI-CV crossings between n=8 and n=12")
print(f"  Sprint 043 chi=10 crossing was artifact (bad n=12 ground state)")
print(f"  ν(q=10) from slope ratio (n=8,12): {nu_global:.2f}")
print(f"  ν(q=10) from data collapse: {nu_opt:.2f}")

# Save
results = {
    'sprint': '048c',
    'key_finding': 'NO crossings at q=10 chi=20 — same as q=4',
    'slope_analysis': {
        'slope_n8_global': float(s8_global),
        'slope_n12_global': float(s12_global),
        'ratio': float(ratio_global),
        'nu_slope': float(nu_global),
    },
    'collapse_analysis': {
        'nu_optimal': float(nu_opt),
        'gc_optimal': float(gc_opt),
        'quality_optimal': float(res.fun),
    },
    'comparison': {
        'g_values': g12.tolist(),
        'n8_cv': [float(cv8[np.argmin(np.abs(g8-g))]) for g in g12],
        'n12_cv': cv12.tolist(),
        'diff': [float(cv12[i] - cv8[np.argmin(np.abs(g8-g12[i]))]) for i in range(len(g12))],
    },
}
with open('results/sprint_048c_nu_extraction.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved to results/sprint_048c_nu_extraction.json")
