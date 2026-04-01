"""
Sprint 038a: Data collapse of TFIM MI-CV curves.

Tests whether MI-CV obeys finite-size scaling: CV(h, n) = F((h - h_c) * n^{1/nu}).
If curves for different n collapse onto one universal function, MI-CV is in Ising
universality class.

Consolidates all prior TFIM data (sprints 036-037) and:
1. Tests collapse with nu=1.0 (exact Ising) vs nu=0.755 (finite-size fit)
2. Quantifies collapse quality via residual spread after interpolation
3. Optimizes nu to find best-collapse value
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
import json, time

t0 = time.time()

# ========== Consolidate all TFIM MI-CV data ==========
# Sources: 036a (n=8,16), 036a2 (n=32), 037a (n=8,12,16), 037a2 (n=24,32), 037a3 (n=16,24,32)

all_data = {}  # {n: {h: cv}}

# Load 036a
with open('results/sprint_036a_tfim_micv.json') as f:
    d036a = json.load(f)
for n_str, points in d036a['data'].items():
    n = int(n_str)
    if n not in all_data:
        all_data[n] = {}
    for h_str, vals in points.items():
        h = float(h_str)
        all_data[n][h] = vals['mi_cv']

# Load 036a2
with open('results/sprint_036a2_tfim_n32_50.json') as f:
    d036a2 = json.load(f)
for n_str, points in d036a2['data'].items():
    n = int(n_str)
    if n not in all_data:
        all_data[n] = {}
    for h_str, vals in points.items():
        h = float(h_str)
        all_data[n][h] = vals['mi_cv']

# Load 037a
with open('results/sprint_037a_tfim_crossing.json') as f:
    d037a = json.load(f)
for n_str, points in d037a['sizes'].items():
    n = int(n_str)
    if n not in all_data:
        all_data[n] = {}
    for h_str, vals in points.items():
        h = float(h_str)
        all_data[n][h] = vals['cv']

# Load 037a2
with open('results/sprint_037a2_tfim_large.json') as f:
    d037a2 = json.load(f)
for n_str, points in d037a2['sizes'].items():
    n = int(n_str)
    if n not in all_data:
        all_data[n] = {}
    for h_str, vals in points.items():
        h = float(h_str)
        all_data[n][h] = vals['cv']

# Load 037a3
with open('results/sprint_037a3_tfim_fill.json') as f:
    d037a3 = json.load(f)
for key, vals in d037a3['new_points'].items():
    n = vals['n']
    h = vals['h_J']
    if n not in all_data:
        all_data[n] = {}
    all_data[n][h] = vals['cv']

# Print consolidated data summary
print("=== Sprint 038a: TFIM MI-CV Data Collapse ===\n")
print("Consolidated data:")
for n in sorted(all_data.keys()):
    h_vals = sorted(all_data[n].keys())
    print(f"  n={n:3d}: {len(h_vals)} points, h=[{min(h_vals):.2f}, {max(h_vals):.2f}]")

# ========== Prepare scaled data ==========
h_c = 1.0  # exact TFIM critical point

def scale_data(nu, h_c=1.0):
    """Scale h -> x = (h - h_c) * n^{1/nu} for all sizes."""
    scaled = {}
    for n in sorted(all_data.keys()):
        h_arr = np.array(sorted(all_data[n].keys()))
        cv_arr = np.array([all_data[n][h] for h in sorted(all_data[n].keys())])
        x_arr = (h_arr - h_c) * n**(1.0/nu)
        scaled[n] = (x_arr, cv_arr)
    return scaled

def collapse_quality(nu, h_c=1.0, sizes=None):
    """Measure collapse quality: lower = better.

    For each pair of sizes, interpolate both onto common x range,
    compute mean squared difference. Return average across all pairs.
    """
    scaled = scale_data(nu, h_c)
    if sizes is None:
        sizes = sorted(scaled.keys())

    total_msd = 0.0
    n_pairs = 0

    for i in range(len(sizes)):
        for j in range(i+1, len(sizes)):
            n1, n2 = sizes[i], sizes[j]
            x1, cv1 = scaled[n1]
            x2, cv2 = scaled[n2]

            # Find overlapping x range
            x_min = max(x1.min(), x2.min())
            x_max = min(x1.max(), x2.max())

            if x_max <= x_min:
                continue

            # Interpolate both onto common grid
            f1 = interp1d(x1, cv1, kind='linear', fill_value='extrapolate')
            f2 = interp1d(x2, cv2, kind='linear', fill_value='extrapolate')

            x_common = np.linspace(x_min, x_max, 50)
            cv1_interp = f1(x_common)
            cv2_interp = f2(x_common)

            msd = np.mean((cv1_interp - cv2_interp)**2)
            total_msd += msd
            n_pairs += 1

    if n_pairs == 0:
        return 1e10
    return total_msd / n_pairs

# ========== Test specific nu values ==========
print("\n=== Collapse Quality (lower = better) ===")
test_nus = [0.5, 0.6, 0.7, 0.755, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0]
sizes_to_use = [n for n in sorted(all_data.keys()) if n >= 8]

results = {
    'experiment': '038a',
    'h_c': h_c,
    'sizes': {str(n): {str(h): cv for h, cv in sorted(all_data[n].items())}
              for n in sorted(all_data.keys())},
    'collapse_quality': {},
}

for nu in test_nus:
    q = collapse_quality(nu, h_c, sizes_to_use)
    print(f"  nu={nu:.3f}: quality = {q:.6f}")
    results['collapse_quality'][str(nu)] = float(q)

# ========== Optimize nu ==========
print("\n=== Optimizing nu ===")
opt = minimize_scalar(lambda nu: collapse_quality(nu, h_c, sizes_to_use),
                      bounds=(0.3, 3.0), method='bounded')
nu_opt = opt.x
q_opt = opt.fun
print(f"  Best nu = {nu_opt:.4f}, quality = {q_opt:.6f}")
results['nu_optimal'] = float(nu_opt)
results['quality_optimal'] = float(q_opt)

# Compare to key values
q_ising = collapse_quality(1.0, h_c, sizes_to_use)
q_fit = collapse_quality(0.755, h_c, sizes_to_use)
print(f"  Ising nu=1.0: quality = {q_ising:.6f}")
print(f"  Fit nu=0.755: quality = {q_fit:.6f}")
print(f"  Ratio Ising/optimal = {q_ising/q_opt:.2f}")
print(f"  Ratio fit/optimal = {q_fit/q_opt:.2f}")
results['quality_ising'] = float(q_ising)
results['quality_fit'] = float(q_fit)

# ========== Print collapsed curves for best nu ==========
print(f"\n=== Collapsed Curves (nu={nu_opt:.3f}) ===")
scaled_opt = scale_data(nu_opt, h_c)
for n in sorted(scaled_opt.keys()):
    x, cv = scaled_opt[n]
    print(f"\nn={n}:")
    for xi, cvi in zip(x, cv):
        print(f"  x={xi:8.3f}  CV={cvi:.4f}")

# ========== Also print for nu=1.0 ==========
print(f"\n=== Collapsed Curves (nu=1.0, Ising) ===")
scaled_ising = scale_data(1.0, h_c)
for n in sorted(scaled_ising.keys()):
    x, cv = scaled_ising[n]
    print(f"\nn={n}:")
    for xi, cvi in zip(x, cv):
        print(f"  x={xi:8.3f}  CV={cvi:.4f}")

# ========== Size-pair collapse matrix ==========
print("\n=== Pairwise Collapse Quality (nu_opt) ===")
sizes = sorted(all_data.keys())
pair_quality = {}
for i in range(len(sizes)):
    for j in range(i+1, len(sizes)):
        q = collapse_quality(nu_opt, h_c, [sizes[i], sizes[j]])
        pair_quality[f"{sizes[i]}_{sizes[j]}"] = float(q)
        print(f"  ({sizes[i]:2d},{sizes[j]:2d}): {q:.6f}")
results['pairwise_quality'] = pair_quality

# ========== Also optimize h_c jointly ==========
print("\n=== Joint optimization of h_c and nu ===")
from scipy.optimize import minimize

def joint_quality(params):
    nu, hc = params
    if nu < 0.2 or nu > 5.0 or hc < 0.5 or hc > 1.5:
        return 1e10
    return collapse_quality(nu, hc, sizes_to_use)

res = minimize(joint_quality, [nu_opt, 1.0], method='Nelder-Mead',
               options={'xatol': 0.001, 'fatol': 1e-8})
nu_joint, hc_joint = res.x
q_joint = res.fun
print(f"  Best: nu={nu_joint:.4f}, h_c={hc_joint:.4f}, quality={q_joint:.6f}")
results['joint_optimal'] = {'nu': float(nu_joint), 'h_c': float(hc_joint),
                            'quality': float(q_joint)}

# Compare: is Ising (nu=1, hc=1) close to joint optimum?
q_exact_ising = collapse_quality(1.0, 1.0, sizes_to_use)
print(f"  Exact Ising (nu=1, hc=1): quality = {q_exact_ising:.6f}")
print(f"  Ratio exact_Ising/joint_optimal = {q_exact_ising/q_joint:.2f}")

# Save results
results['total_runtime'] = time.time() - t0
with open('results/sprint_038a_data_collapse.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nTotal runtime: {time.time()-t0:.1f}s")
print("Results saved to results/sprint_038a_data_collapse.json")
