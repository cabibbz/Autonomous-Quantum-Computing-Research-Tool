#!/usr/bin/env python3
"""Sprint 046f: q=4 Potts constrained collapse — fix g_c near expected value.

The free collapse finds g_c≈0.53, far from expected g_c≈0.89. This may be
fitting to the CV minimum, not the phase transition. Test with g_c constrained.

Also compute slopes in the transition region (g=0.80-1.00) near true g_c.
"""
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.interpolate import interp1d
import json

# Consolidated data (all direct MPS, chi=20)
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

    x_grid = np.linspace(x_min, x_max, 100)
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

# ===== 1. Constrained collapse: fix g_c, optimize nu =====
print("=== CONSTRAINED COLLAPSE: Fixed g_c, optimize nu ===")
print("(Using g=0.70-1.05 transition region)\n")
g_range = (0.68, 1.07)

gc_tests = [0.70, 0.75, 0.80, 0.85, 0.89, 0.90, 0.95, 1.00]
constrained_results = {}
for gc in gc_tests:
    best_q = 1e6; best_nu = None
    for nu0 in np.arange(0.3, 10.0, 0.1):
        q = collapse_quality(nu0, gc, sizes_data, g_range)
        if q < best_q:
            best_q = q
            best_nu = nu0
    # Refine
    res = minimize_scalar(lambda nu: collapse_quality(nu, gc, sizes_data, g_range),
                         bounds=(max(0.1, best_nu-0.5), best_nu+0.5), method='bounded')
    constrained_results[gc] = {'nu': res.x, 'quality': res.fun}
    print(f"  g_c={gc:.2f}: best nu={res.x:.3f}, quality={res.fun:.6f}")

# ===== 2. Slope analysis near expected g_c ≈ 0.89 =====
print("\n=== SLOPE ANALYSIS near g_c ≈ 0.89 ===")
# Compute local slopes in g=0.80-1.00
for n, data in sorted(sizes_data.items()):
    gs = sorted([g for g in data.keys() if 0.78 <= g <= 1.02])
    if len(gs) >= 2:
        cvs = [data[g] for g in gs]
        slopes = []
        for i in range(len(gs)-1):
            s = (cvs[i+1] - cvs[i]) / (gs[i+1] - gs[i])
            slopes.append((gs[i], gs[i+1], s))
        # Average slope near g_c=0.89
        near_gc = [(g1, g2, s) for g1, g2, s in slopes if g1 >= 0.83 and g2 <= 0.97]
        if near_gc:
            avg_s = np.mean([s for _, _, s in near_gc])
            print(f"  n={n:2d}: avg slope near g_c = {avg_s:.4f}")
            for g1, g2, s in near_gc:
                print(f"    g={g1:.2f}-{g2:.2f}: slope={s:.4f}")

# Slope ratios near g_c
print("\n--- Slope ratios near g_c ---")
# For each size pair, compare slopes in similar g ranges
for n_pair in [(8, 12), (12, 16), (8, 16)]:
    n1, n2 = n_pair
    d1, d2 = sizes_data[n1], sizes_data[n2]
    common_gs = sorted(set(d1.keys()) & set(d2.keys()))
    gc_region = [g for g in common_gs if 0.83 <= g <= 0.97]
    if len(gc_region) >= 2:
        slopes1 = []; slopes2 = []
        for i in range(len(gc_region)-1):
            g1, g2 = gc_region[i], gc_region[i+1]
            s1 = (d1[g2] - d1[g1]) / (g2 - g1)
            s2 = (d2[g2] - d2[g1]) / (g2 - g1)
            slopes1.append(s1); slopes2.append(s2)
        avg_s1, avg_s2 = np.mean(slopes1), np.mean(slopes2)
        ratio = avg_s2 / avg_s1
        if ratio > 0:
            exp = np.log(ratio) / np.log(n2/n1)
            print(f"  n={n1},{n2}: slope ratio = {ratio:.4f}, exponent = {exp:.3f} => nu ≈ {1/exp:.2f}")

# ===== 3. CV gap analysis =====
print("\n=== CV GAP: n=12-n=8 and n=16-n=12 ===")
print(f"{'g':>6s} {'CV8':>7s} {'CV12':>7s} {'CV16':>7s} {'d12-8':>7s} {'d16-12':>7s} {'d16-8':>7s}")
common_g = sorted(set(n8.keys()) & set(n12.keys()) & set(n16.keys()))
for g in common_g:
    cv8, cv12, cv16 = n8[g], n12[g], n16[g]
    print(f"  {g:.2f} {cv8:.4f} {cv12:.4f} {cv16:.4f} {cv12-cv8:+.4f} {cv16-cv12:+.4f} {cv16-cv8:+.4f}")

# ===== 4. Derivative scaling (numerical dCV/dg at g_c) =====
print("\n=== DERIVATIVE at g_c ≈ 0.89 ===")
gc_est = 0.89
for n, data in sorted(sizes_data.items()):
    gs = sorted(data.keys())
    cvs = [data[g] for g in gs]
    # Find bracketing points
    for i in range(len(gs)-1):
        if gs[i] <= gc_est <= gs[i+1]:
            deriv = (cvs[i+1] - cvs[i]) / (gs[i+1] - gs[i])
            print(f"  n={n:2d}: dCV/dg|_gc ≈ {deriv:.4f} (between g={gs[i]:.2f} and g={gs[i+1]:.2f})")
            break

# Save
output = {
    'sprint': '046f',
    'description': 'q=4 Potts constrained collapse and slope analysis',
    'constrained_collapse': {str(gc): v for gc, v in constrained_results.items()},
}
with open('results/sprint_046f_q4_constrained.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)

print("\nResults saved to results/sprint_046f_q4_constrained.json")
