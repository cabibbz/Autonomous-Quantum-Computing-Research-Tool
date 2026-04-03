#!/usr/bin/env python3
"""Sprint 046e: q=4 Potts data collapse with 3 system sizes (n=8,12,16).

All data from direct MPS contraction, chi=20.
Joint optimization of (nu, g_c) for best scaling collapse.

Key questions:
  - Where does nu(4) fall? Between nu(3)=5/6 and nu(5)=2?
  - Does q=4 marginality (log corrections in 2D) manifest in 1D quantum?
"""
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.interpolate import interp1d
import json

# Consolidated data from exps 046a-d (all direct MPS, chi=20)
n8  = {0.20: 0.3556, 0.30: 0.3951, 0.40: 0.1739, 0.50: 0.1599, 0.60: 0.3022,
       0.70: 0.4270, 0.75: 0.4820, 0.80: 0.5329, 0.85: 0.5800, 0.88: 0.6067,
       0.90: 0.6238, 0.93: 0.6487, 0.95: 0.6646, 1.00: 0.7027, 1.05: 0.7383}
n12 = {0.20: 0.3080, 0.30: 0.3982, 0.40: 0.1711, 0.50: 0.2043, 0.60: 0.3516,
       0.70: 0.4856, 0.75: 0.5476, 0.80: 0.6065, 0.85: 0.6626, 0.88: 0.6950,
       0.90: 0.7160, 0.93: 0.7468, 0.95: 0.7669, 1.00: 0.8153, 1.05: 0.8613}
n16 = {0.40: 0.1634, 0.50: 0.2085, 0.60: 0.3515, 0.70: 0.4880, 0.80: 0.6164,
       0.85: 0.6775, 0.90: 0.7365, 0.95: 0.7935, 1.00: 0.8485, 1.05: 0.9015}

sizes_data = {8: n8, 12: n12, 16: n16}

def collapse_quality(nu, g_c, sizes_data, g_range):
    """Compute collapse quality for given (nu, g_c)."""
    rescaled = {}
    for n, gv_dict in sizes_data.items():
        gs = sorted([g for g in gv_dict.keys() if g_range[0] <= g <= g_range[1]])
        if len(gs) < 3:
            return 1e6
        cvs = [gv_dict[g] for g in gs]
        xs = [(g - g_c) * n**(1.0/nu) for g in gs]
        rescaled[n] = (np.array(xs), np.array(cvs))

    all_x_ranges = [(xs.min(), xs.max()) for xs, _ in rescaled.values()]
    x_min = max(r[0] for r in all_x_ranges)
    x_max = min(r[1] for r in all_x_ranges)
    if x_min >= x_max:
        return 1e6

    n_grid = 100
    x_grid = np.linspace(x_min, x_max, n_grid)
    curves = []
    for n, (xs, cvs) in rescaled.items():
        try:
            f = interp1d(xs, cvs, kind='linear', fill_value='extrapolate')
            curves.append(f(x_grid))
        except:
            return 1e6

    curves = np.array(curves)
    mean_curve = curves.mean(axis=0)
    variance = np.mean(np.var(curves, axis=0))
    mean_sq = np.mean(mean_curve**2)
    if mean_sq < 1e-10:
        return 1e6
    return variance / mean_sq

# ===== Analysis 1: Transition region (g=0.60 to 1.05) =====
print("=== q=4 Potts Data Collapse (n=8, 12, 16, direct MPS, chi=20) ===\n")
print("--- Analysis 1: Transition region g=0.60-1.05 ---")
g_range = (0.58, 1.07)

nu_candidates = {
    'Ising (nu=1)': 1.0,
    'Potts q=3 (nu=5/6)': 5/6,
    'nu=2/3': 2/3,
    'nu=1/2': 0.5,
    'nu=3/2': 1.5,
    'nu=2': 2.0,
    'nu=3': 3.0,
    'nu=4': 4.0,
}

results1 = {}
for name, nu in nu_candidates.items():
    def obj_gc(gc, nu=nu):
        return collapse_quality(nu, gc, sizes_data, g_range)
    res = minimize_scalar(obj_gc, bounds=(0.50, 1.10), method='bounded')
    results1[name] = {'nu': nu, 'g_c': res.x, 'quality': res.fun}
    print(f"  {name:25s}: g_c={res.x:.4f}, quality={res.fun:.6f}")

# Free joint optimization
print("\n--- Free joint optimization (transition region) ---")
best = None
for nu0 in [0.4, 0.5, 0.6, 0.7, 0.83, 1.0, 1.2, 1.5, 2.0, 3.0, 4.0, 5.0]:
    for gc0 in [0.60, 0.70, 0.80, 0.90, 1.00]:
        def obj_joint(params):
            nu, gc = params
            if nu < 0.1 or nu > 20.0 or gc < 0.30 or gc > 1.20:
                return 1e6
            return collapse_quality(nu, gc, sizes_data, g_range)
        res = minimize(obj_joint, [nu0, gc0], method='Nelder-Mead',
                      options={'xatol': 0.001, 'fatol': 1e-8, 'maxiter': 2000})
        if best is None or res.fun < best.fun:
            best = res

nu_opt, gc_opt = best.x
q_opt = best.fun
print(f"  Optimal: nu={nu_opt:.4f}, g_c={gc_opt:.4f}, quality={q_opt:.6f}")
results1['free_fit'] = {'nu': nu_opt, 'g_c': gc_opt, 'quality': q_opt}

# Rank
print("\n=== RANKING (transition region) ===")
ranked1 = sorted(results1.items(), key=lambda x: x[1]['quality'])
best_q1 = ranked1[0][1]['quality']
for i, (name, r) in enumerate(ranked1):
    ratio = r['quality'] / best_q1 if best_q1 > 0 else float('inf')
    print(f"  {i+1}. {name:25s}: quality={r['quality']:.6f} (ratio={ratio:.3f}), "
          f"nu={r['nu']:.4f}, g_c={r['g_c']:.4f}")

# ===== Analysis 2: Broader range including minimum (g=0.40 to 1.05) =====
print("\n--- Analysis 2: Broad range g=0.40-1.05 ---")
g_range2 = (0.38, 1.07)

results2 = {}
for name, nu in nu_candidates.items():
    def obj_gc(gc, nu=nu):
        return collapse_quality(nu, gc, sizes_data, g_range2)
    res = minimize_scalar(obj_gc, bounds=(0.40, 1.10), method='bounded')
    results2[name] = {'nu': nu, 'g_c': res.x, 'quality': res.fun}
    print(f"  {name:25s}: g_c={res.x:.4f}, quality={res.fun:.6f}")

best2 = None
for nu0 in [0.4, 0.5, 0.6, 0.7, 0.83, 1.0, 1.2, 1.5, 2.0, 3.0, 4.0, 5.0, 8.0, 10.0]:
    for gc0 in [0.50, 0.60, 0.70, 0.80, 0.90, 1.00]:
        def obj_joint(params):
            nu, gc = params
            if nu < 0.1 or nu > 20.0 or gc < 0.30 or gc > 1.20:
                return 1e6
            return collapse_quality(nu, gc, sizes_data, g_range2)
        res = minimize(obj_joint, [nu0, gc0], method='Nelder-Mead',
                      options={'xatol': 0.001, 'fatol': 1e-8, 'maxiter': 2000})
        if best2 is None or res.fun < best2.fun:
            best2 = res

nu_opt2, gc_opt2 = best2.x
print(f"\n  Free fit (broad): nu={nu_opt2:.4f}, g_c={gc_opt2:.4f}, quality={best2.fun:.6f}")

# ===== Crossing point analysis =====
print("\n=== CROSSING POINTS ===")
size_pairs = [(8, 12), (8, 16), (12, 16)]
crossing_gcs = {}
for n1, n2 in size_pairs:
    d1, d2 = sizes_data[n1], sizes_data[n2]
    common_gs = sorted(set(d1.keys()) & set(d2.keys()))
    deltas = [(g, d2[g] - d1[g]) for g in common_gs]
    crossings = []
    for i in range(len(deltas)-1):
        g1, delta1 = deltas[i]
        g2, delta2 = deltas[i+1]
        if delta1 * delta2 < 0:
            gc = g1 + (g2-g1) * (-delta1) / (delta2-delta1)
            crossings.append(gc)
    if crossings:
        crossing_gcs[(n1, n2)] = crossings
        for gc in crossings:
            print(f"  n={n1},{n2}: g_c={gc:.4f}")

# ===== Slope analysis =====
print("\n=== SLOPE ANALYSIS (transition region g=0.60-1.05) ===")
slopes = {}
for n, data in sorted(sizes_data.items()):
    gs = sorted([g for g in data.keys() if 0.60 <= g <= 1.05])
    cvs = [data[g] for g in gs]
    max_slope = 0
    for i in range(len(gs)-1):
        s = abs(cvs[i+1] - cvs[i]) / (gs[i+1] - gs[i])
        if s > max_slope:
            max_slope = s
            slopes[n] = {'slope': s, 'g_range': (gs[i], gs[i+1])}
    if n in slopes:
        print(f"  n={n:2d}: max |slope|={slopes[n]['slope']:.3f} "
              f"between g={slopes[n]['g_range'][0]:.2f}-{slopes[n]['g_range'][1]:.2f}")

ns = sorted(slopes.keys())
if len(ns) >= 2:
    for i in range(len(ns)-1):
        n1, n2 = ns[i], ns[i+1]
        s1, s2 = slopes[n1]['slope'], slopes[n2]['slope']
        if s1 > 0:
            ratio = s2 / s1
            exp = np.log(ratio) / np.log(n2/n1)
            print(f"  Slope ratio n{n2}/n{n1} = {ratio:.3f} => exponent = {exp:.3f} ≈ 1/nu => nu ≈ {1/exp:.2f}")
    if len(ns) >= 3:
        s_arr = np.array([slopes[n]['slope'] for n in ns])
        n_arr = np.array(ns, dtype=float)
        log_s = np.log(s_arr)
        log_n = np.log(n_arr)
        coeffs = np.polyfit(log_n, log_s, 1)
        print(f"  Power-law fit: slope ~ n^{coeffs[0]:.3f} => nu = {1/coeffs[0]:.2f}")

# ===== n=8,12 only analysis (more data points) =====
print("\n=== TWO-SIZE COLLAPSE (n=8,12 only, g=0.50-1.05) ===")
sizes_2 = {8: n8, 12: n12}
g_range3 = (0.48, 1.07)
best3 = None
for nu0 in [0.3, 0.5, 0.7, 0.83, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0]:
    for gc0 in [0.60, 0.70, 0.80, 0.90, 1.00]:
        def obj_joint(params):
            nu, gc = params
            if nu < 0.1 or nu > 20.0 or gc < 0.30 or gc > 1.20:
                return 1e6
            return collapse_quality(nu, gc, sizes_2, g_range3)
        res = minimize(obj_joint, [nu0, gc0], method='Nelder-Mead',
                      options={'xatol': 0.001, 'fatol': 1e-8, 'maxiter': 2000})
        if best3 is None or res.fun < best3.fun:
            best3 = res
print(f"  n=8,12 free fit: nu={best3.x[0]:.4f}, g_c={best3.x[1]:.4f}, quality={best3.fun:.6f}")

# Save
output = {
    'sprint': '046e',
    'description': 'q=4 Potts data collapse — nu extraction',
    'data': {
        'n8': {str(k): v for k, v in n8.items()},
        'n12': {str(k): v for k, v in n12.items()},
        'n16': {str(k): v for k, v in n16.items()},
    },
    'transition_region': {
        'g_range': [0.60, 1.05],
        'results': {k: v for k, v in results1.items()},
        'ranking': [(name, r['quality'], r['nu'], r['g_c']) for name, r in ranked1],
        'optimal_nu': float(nu_opt),
        'optimal_gc': float(gc_opt),
    },
    'broad_region': {
        'optimal_nu': float(nu_opt2),
        'optimal_gc': float(gc_opt2),
        'quality': float(best2.fun),
    },
    'two_size': {
        'optimal_nu': float(best3.x[0]),
        'optimal_gc': float(best3.x[1]),
        'quality': float(best3.fun),
    },
    'crossing_points': {f'{k[0]},{k[1]}': v for k, v in crossing_gcs.items()},
}
with open('results/sprint_046e_q4_collapse.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)

print("\nResults saved to results/sprint_046e_q4_collapse.json")
