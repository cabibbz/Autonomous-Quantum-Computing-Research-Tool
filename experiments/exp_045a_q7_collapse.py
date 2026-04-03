#!/usr/bin/env python3
"""Sprint 045a: Data collapse for q=7 Potts using existing n=8,12 MI-CV data.

Test scaling ansatz: CV(g, n) = F((g - g_c) * n^{1/nu})
Jointly optimize (nu, g_c) to find the best collapse.

Existing data (chi=20):
  n=8:  g=0.20,0.25,0.30,0.35,0.40,0.50,0.60
  n=12: g=0.25,0.30,0.35

Compare nu candidates: 1 (Ising), 5/6 (Potts q=3), 2/3 (mean-field), free fit.
"""
import numpy as np
from scipy.optimize import minimize_scalar, minimize
from scipy.interpolate import interp1d
import json

# Load existing data
n8_data = {0.20: 0.5115, 0.25: 0.5292, 0.30: 0.4378, 0.35: 0.4467, 0.40: 0.3070, 0.50: 0.2365, 0.60: 0.1135}
n12_data = {0.25: 0.5079, 0.30: 0.5412, 0.35: 0.4995}

# Build arrays for each size
sizes = {8: n8_data, 12: n12_data}

# Data collapse quality: for each pair of sizes, measure how well the rescaled
# curves overlap using interpolation
def collapse_quality(nu, g_c, sizes_data):
    """Compute collapse quality by measuring overlap of rescaled curves."""
    rescaled = {}
    for n, gv_dict in sizes_data.items():
        gs = sorted(gv_dict.keys())
        cvs = [gv_dict[g] for g in gs]
        xs = [(g - g_c) * n**(1.0/nu) for g in gs]
        rescaled[n] = (np.array(xs), np.array(cvs))

    # Find overlap region in rescaled x
    all_x_ranges = [(xs.min(), xs.max()) for xs, _ in rescaled.values()]
    x_min = max(r[0] for r in all_x_ranges)
    x_max = min(r[1] for r in all_x_ranges)

    if x_min >= x_max:
        return 1e6  # no overlap

    # Interpolate each curve on common x grid
    x_grid = np.linspace(x_min, x_max, 50)
    curves = []
    for n, (xs, cvs) in rescaled.items():
        if len(xs) < 2:
            return 1e6
        f = interp1d(xs, cvs, kind='linear', fill_value='extrapolate')
        curves.append(f(x_grid))

    # Quality = variance of curves at each x point, averaged
    curves = np.array(curves)
    mean_curve = curves.mean(axis=0)
    variance = np.mean((curves - mean_curve)**2)
    # Normalize by mean CV squared to get relative quality
    mean_cv_sq = np.mean(mean_curve**2)
    if mean_cv_sq < 1e-10:
        return 1e6
    return variance / mean_cv_sq

# Test specific nu values with g_c optimization
print("=== q=7 Potts Data Collapse (n=8, n=12) ===\n")

nu_candidates = {
    'Ising (nu=1)': 1.0,
    'Potts q=3 (nu=5/6)': 5/6,
    'Mean-field (nu=2/3)': 2/3,
    'nu=1/2': 0.5,
    'nu=3/2': 1.5,
}

results = {}
for name, nu in nu_candidates.items():
    # Optimize g_c for this fixed nu
    def obj_gc(gc):
        return collapse_quality(nu, gc, sizes)

    res = minimize_scalar(obj_gc, bounds=(0.20, 0.35), method='bounded')
    q_val = res.fun
    gc_opt = res.x
    results[name] = {'nu': nu, 'g_c': gc_opt, 'quality': q_val}
    print(f"  {name:25s}: g_c={gc_opt:.4f}, quality={q_val:.6f}")

# Free joint optimization
print("\n--- Free joint optimization ---")
def obj_joint(params):
    nu, gc = params
    if nu < 0.1 or nu > 5.0 or gc < 0.15 or gc > 0.40:
        return 1e6
    return collapse_quality(nu, gc, sizes)

best = None
for nu0 in [0.5, 0.7, 0.83, 1.0, 1.2, 1.5, 2.0]:
    for gc0 in [0.24, 0.26, 0.28, 0.30]:
        res = minimize(obj_joint, [nu0, gc0], method='Nelder-Mead',
                      options={'xatol': 0.001, 'fatol': 1e-8, 'maxiter': 500})
        if best is None or res.fun < best.fun:
            best = res

nu_opt, gc_opt = best.x
q_opt = best.fun
print(f"  Optimal: nu={nu_opt:.4f}, g_c={gc_opt:.4f}, quality={q_opt:.6f}")
results['free_fit'] = {'nu': nu_opt, 'g_c': gc_opt, 'quality': q_opt}

# Rank by quality
print("\n=== RANKING ===")
ranked = sorted(results.items(), key=lambda x: x[1]['quality'])
best_q = ranked[0][1]['quality']
for i, (name, r) in enumerate(ranked):
    ratio = r['quality'] / best_q if best_q > 0 else float('inf')
    print(f"  {i+1}. {name:25s}: quality={r['quality']:.6f} (ratio={ratio:.3f}), "
          f"nu={r['nu']:.4f}, g_c={r['g_c']:.4f}")

# Slope analysis: compute MI-CV slope at crossing
print("\n=== SLOPE ANALYSIS ===")
for n, gv_dict in sorted(sizes.items()):
    gs = sorted(gv_dict.keys())
    cvs = [gv_dict[g] for g in gs]
    # Estimate slope at g nearest to g_c
    gc_est = ranked[0][1]['g_c']
    # Find two bracketing points
    for i in range(len(gs)-1):
        if gs[i] <= gc_est <= gs[i+1]:
            slope = (cvs[i+1] - cvs[i]) / (gs[i+1] - gs[i])
            print(f"  n={n:2d}: slope at g_c≈{gc_est:.3f} = {slope:.3f} (between g={gs[i]:.2f} and {gs[i+1]:.2f})")
            break

# Check: does slope ratio give 1/nu?
slopes = {}
for n, gv_dict in sorted(sizes.items()):
    gs = sorted(gv_dict.keys())
    cvs = [gv_dict[g] for g in gs]
    gc_est = ranked[0][1]['g_c']
    for i in range(len(gs)-1):
        if gs[i] <= gc_est <= gs[i+1]:
            slopes[n] = (cvs[i+1] - cvs[i]) / (gs[i+1] - gs[i])
            break

if 8 in slopes and 12 in slopes:
    ratio = slopes[12] / slopes[8]
    log_ratio = np.log(ratio) / np.log(12/8)
    print(f"\n  Slope ratio n12/n8 = {ratio:.3f}")
    print(f"  => exponent (log(ratio)/log(12/8)) = {log_ratio:.3f}")
    print(f"  Compare: 1/nu_Ising=1.0, 1/nu_Potts3=1.2, free fit 1/nu={1/nu_opt:.3f}")

# Save results
output = {
    'sprint': '045a',
    'description': 'q=7 Potts data collapse nu extraction',
    'input_data': {'n8': n8_data, 'n12': {str(k): v for k, v in n12_data.items()}},
    'collapse_results': {k: v for k, v in results.items()},
    'ranking': [(name, r['quality']) for name, r in ranked],
    'optimal_nu': nu_opt,
    'optimal_gc': gc_opt,
    'slopes': {str(k): v for k, v in slopes.items()},
}
with open('results/sprint_045a_q7_collapse.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)

print("\nResults saved to results/sprint_045a_q7_collapse.json")
