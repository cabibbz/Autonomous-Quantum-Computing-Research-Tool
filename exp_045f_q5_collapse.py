#!/usr/bin/env python3
"""Sprint 045f: q=5 Potts data collapse with 3 system sizes (n=8,12,16).

All data from direct MPS contraction, chi=20.
Joint optimization of (nu, g_c) for best scaling collapse.

Key question: What is nu(q=5)? Compare with:
  q=2 (Ising): nu=1
  q=3 (3-state Potts): nu=5/6 ≈ 0.833
"""
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.interpolate import interp1d
import json

# Corrected consolidated data (all direct MPS, chi=20)
n8  = {0.30: 0.3732, 0.35: 0.3554, 0.38: 0.3141, 0.40: 0.2755, 0.42: 0.2325,
       0.45: 0.1700, 0.50: 0.1072, 0.55: 0.1324}
n12 = {0.30: 0.3528, 0.35: 0.3671, 0.38: 0.3399, 0.40: 0.2993, 0.42: 0.2456,
       0.45: 0.1677, 0.50: 0.1212}
n16 = {0.30: 0.4818, 0.35: 0.3559, 0.38: 0.3458, 0.40: 0.3100, 0.42: 0.2505,
       0.45: 0.1604, 0.50: 0.1227}

sizes_data = {8: n8, 12: n12, 16: n16}

# Restrict to transition region (avoid deep disordered CV minimum artifacts)
# Focus on g=0.30 to 0.50 where all three sizes have data
g_range = (0.28, 0.52)

def collapse_quality(nu, g_c, sizes_data, g_range=(0.28, 0.52)):
    """Compute collapse quality for given (nu, g_c).
    Rescale x = (g - g_c) * n^(1/nu), then measure spread of curves."""
    rescaled = {}
    for n, gv_dict in sizes_data.items():
        gs = sorted([g for g in gv_dict.keys() if g_range[0] <= g <= g_range[1]])
        if len(gs) < 3:
            return 1e6
        cvs = [gv_dict[g] for g in gs]
        xs = [(g - g_c) * n**(1.0/nu) for g in gs]
        rescaled[n] = (np.array(xs), np.array(cvs))

    # Find overlap region in rescaled x
    all_x_ranges = [(xs.min(), xs.max()) for xs, _ in rescaled.values()]
    x_min = max(r[0] for r in all_x_ranges)
    x_max = min(r[1] for r in all_x_ranges)

    if x_min >= x_max:
        return 1e6

    # Interpolate each curve on common x grid
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
    # Quality = mean variance across curves / mean^2
    mean_curve = curves.mean(axis=0)
    variance = np.mean(np.var(curves, axis=0))
    mean_sq = np.mean(mean_curve**2)
    if mean_sq < 1e-10:
        return 1e6
    return variance / mean_sq

# Test specific nu values
print("=== q=5 Potts Data Collapse (n=8, 12, 16, direct MPS, chi=20) ===\n")

nu_candidates = {
    'Ising (nu=1)': 1.0,
    'Potts q=3 (nu=5/6)': 5/6,
    'nu=2/3': 2/3,
    'nu=1/2': 0.5,
    'nu=3/2': 1.5,
    'nu=2': 2.0,
}

results = {}
for name, nu in nu_candidates.items():
    def obj_gc(gc, nu=nu):
        return collapse_quality(nu, gc, sizes_data)
    res = minimize_scalar(obj_gc, bounds=(0.25, 0.50), method='bounded')
    q_val = res.fun
    gc_opt = res.x
    results[name] = {'nu': nu, 'g_c': gc_opt, 'quality': q_val}
    print(f"  {name:25s}: g_c={gc_opt:.4f}, quality={q_val:.6f}")

# Free joint optimization with multiple starting points
print("\n--- Free joint optimization ---")
best = None
for nu0 in [0.4, 0.5, 0.6, 0.7, 0.83, 1.0, 1.2, 1.5, 2.0, 3.0]:
    for gc0 in [0.30, 0.33, 0.36, 0.40, 0.45]:
        def obj_joint(params):
            nu, gc = params
            if nu < 0.1 or nu > 10.0 or gc < 0.20 or gc > 0.55:
                return 1e6
            return collapse_quality(nu, gc, sizes_data)
        res = minimize(obj_joint, [nu0, gc0], method='Nelder-Mead',
                      options={'xatol': 0.001, 'fatol': 1e-8, 'maxiter': 1000})
        if best is None or res.fun < best.fun:
            best = res

nu_opt, gc_opt = best.x
q_opt = best.fun
print(f"  Optimal: nu={nu_opt:.4f}, g_c={gc_opt:.4f}, quality={q_opt:.6f}")
results['free_fit'] = {'nu': nu_opt, 'g_c': gc_opt, 'quality': q_opt}

# Rank
print("\n=== RANKING ===")
ranked = sorted(results.items(), key=lambda x: x[1]['quality'])
best_q = ranked[0][1]['quality']
for i, (name, r) in enumerate(ranked):
    ratio = r['quality'] / best_q if best_q > 0 else float('inf')
    print(f"  {i+1}. {name:25s}: quality={r['quality']:.6f} (ratio={ratio:.3f}), "
          f"nu={r['nu']:.4f}, g_c={r['g_c']:.4f}")

# Crossing point analysis
print("\n=== CROSSING POINTS ===")
size_pairs = [(8, 12), (8, 16), (12, 16)]
crossing_gcs = {}
for n1, n2 in size_pairs:
    d1, d2 = sizes_data[n1], sizes_data[n2]
    common_gs = sorted(set(d1.keys()) & set(d2.keys()))
    deltas = [(g, d2[g] - d1[g]) for g in common_gs]
    for i in range(len(deltas)-1):
        g1, delta1 = deltas[i]
        g2, delta2 = deltas[i+1]
        if delta1 * delta2 < 0:
            gc = g1 + (g2-g1) * (-delta1) / (delta2-delta1)
            if gc > 0.28 and gc < 0.50:  # only in transition region
                crossing_gcs[(n1, n2)] = gc
                print(f"  n={n1},{n2}: g_c={gc:.4f} (between g={g1:.2f} and g={g2:.2f})")
                break

# Extract nu from crossing point pairs
if len(crossing_gcs) >= 2:
    print("\n=== NU FROM CROSSING SHIFT ===")
    pairs = list(crossing_gcs.items())
    for i in range(len(pairs)):
        for j in range(i+1, len(pairs)):
            (n1a, n2a), gc_a = pairs[i]
            (n1b, n2b), gc_b = pairs[j]
            if gc_a != gc_b:
                # The crossing point shift relates to 1/nu
                print(f"  ({n1a},{n2a}) gc={gc_a:.4f} vs ({n1b},{n2b}) gc={gc_b:.4f}")

# Slope analysis
print("\n=== SLOPE AT TRANSITION ===")
slopes = {}
gc_est = ranked[0][1]['g_c']
for n, data in sorted(sizes_data.items()):
    gs = sorted(data.keys())
    cvs = [data[g] for g in gs]
    # Find steepest descent region
    max_slope = 0
    for i in range(len(gs)-1):
        if gs[i] >= 0.30 and gs[i+1] <= 0.50:
            s = abs(cvs[i+1] - cvs[i]) / (gs[i+1] - gs[i])
            if s > max_slope:
                max_slope = s
                slopes[n] = {'slope': s, 'g_range': (gs[i], gs[i+1])}
    if n in slopes:
        print(f"  n={n:2d}: max |slope|={slopes[n]['slope']:.3f} between g={slopes[n]['g_range'][0]:.2f}-{slopes[n]['g_range'][1]:.2f}")

# Slope ratio → 1/nu
ns = sorted(slopes.keys())
if len(ns) >= 2:
    for i in range(len(ns)-1):
        n1, n2 = ns[i], ns[i+1]
        s1, s2 = slopes[n1]['slope'], slopes[n2]['slope']
        if s1 > 0:
            ratio = s2 / s1
            exp = np.log(ratio) / np.log(n2/n1)
            print(f"  Slope ratio n{n2}/n{n1} = {ratio:.3f} => exponent = {exp:.3f} ≈ 1/nu")
    # All three
    if len(ns) >= 3:
        s_arr = np.array([slopes[n]['slope'] for n in ns])
        n_arr = np.array(ns, dtype=float)
        # Power law fit: slope ~ n^(1/nu)
        from numpy.polynomial import polynomial as P
        log_s = np.log(s_arr)
        log_n = np.log(n_arr)
        coeffs = np.polyfit(log_n, log_s, 1)
        print(f"  Power-law fit: slope ~ n^{coeffs[0]:.3f} => nu = {1/coeffs[0]:.3f}")

# Save
output = {
    'sprint': '045f',
    'description': 'q=5 Potts data collapse — nu extraction',
    'data': {
        'n8': {str(k): v for k, v in n8.items()},
        'n12': {str(k): v for k, v in n12.items()},
        'n16': {str(k): v for k, v in n16.items()},
    },
    'collapse_results': {k: v for k, v in results.items()},
    'ranking': [(name, r['quality'], r['nu'], r['g_c']) for name, r in ranked],
    'optimal_nu': float(nu_opt),
    'optimal_gc': float(gc_opt),
    'crossing_points': {f'{k[0]},{k[1]}': v for k, v in crossing_gcs.items()},
}
with open('results/sprint_045f_q5_collapse.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)

print("\nResults saved to results/sprint_045f_q5_collapse.json")
