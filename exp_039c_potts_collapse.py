"""
Sprint 039c: Potts MI-CV data collapse — ν=5/6 vs ν=1.

Test whether Potts MI-CV collapse quality favors ν=5/6 (3-state Potts exact)
over ν=1 (Ising exact). This is THE key test: can MI-CV distinguish
universality classes?

Uses all available Potts MI-CV data: n=8 (13 pts), n=12 (7 pts), n=16 (6 pts),
plus n=24 if available.
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar, minimize
import json, time, os

t0 = time.time()

# ====== Consolidate all Potts MI-CV data ======
all_data = {}

# n=8 from 038b2
with open('results/sprint_038b2_potts_fast.json') as f:
    d8 = json.load(f)
all_data[8] = {}
for g_str, info in d8['sizes']['8'].items():
    all_data[8][float(g_str)] = info['cv']

# n=12 from 039b (includes 038b2 originals)
with open('results/sprint_039b_potts_n12_fill.json') as f:
    d12 = json.load(f)
all_data[12] = {}
for g_str, info in d12['n12_complete'].items():
    all_data[12][float(g_str)] = info['cv']

# n=16 from 039a
with open('results/sprint_039a_potts_n16.json') as f:
    d16 = json.load(f)
all_data[16] = {}
for g_str, info in d16['data'].items():
    all_data[16][float(g_str)] = info['cv']

# n=24 if available
if os.path.exists('results/sprint_039b2_potts_n24.json'):
    with open('results/sprint_039b2_potts_n24.json') as f:
        d24 = json.load(f)
    if len(d24.get('data', {})) >= 2:
        all_data[24] = {}
        for g_str, info in d24['data'].items():
            all_data[24][float(g_str)] = info['cv']
        print(f"Loaded n=24: {len(all_data[24])} points")

print("=== Sprint 039c: Potts Data Collapse ===\n")
print("Data summary:")
for n in sorted(all_data.keys()):
    g_vals = sorted(all_data[n].keys())
    print(f"  n={n:3d}: {len(g_vals)} points, g=[{g_vals[0]:.2f}, {g_vals[-1]:.2f}]")

# ====== Collapse quality function ======
def collapse_quality(nu, g_c, sizes):
    """Mean squared difference between interpolated curves after rescaling."""
    total_msd = 0.0
    n_pairs = 0
    for i in range(len(sizes)):
        for j in range(i+1, len(sizes)):
            n1, n2 = sizes[i], sizes[j]
            if n1 not in all_data or n2 not in all_data:
                continue

            g_arr1 = np.array(sorted(all_data[n1].keys()))
            cv_arr1 = np.array([all_data[n1][g] for g in sorted(all_data[n1].keys())])
            x_arr1 = (g_arr1 - g_c) * n1**(1.0/nu)

            g_arr2 = np.array(sorted(all_data[n2].keys()))
            cv_arr2 = np.array([all_data[n2][g] for g in sorted(all_data[n2].keys())])
            x_arr2 = (g_arr2 - g_c) * n2**(1.0/nu)

            x_min = max(x_arr1.min(), x_arr2.min())
            x_max = min(x_arr1.max(), x_arr2.max())
            if x_max <= x_min:
                continue

            try:
                f1 = interp1d(x_arr1, cv_arr1, kind='linear')
                f2 = interp1d(x_arr2, cv_arr2, kind='linear')
                x_common = np.linspace(x_min, x_max, 50)
                msd = np.mean((f1(x_common) - f2(x_common))**2)
                total_msd += msd
                n_pairs += 1
            except:
                continue
    return total_msd / n_pairs if n_pairs > 0 else 1e10

# ====== Test key ν values ======
print("\n=== Collapse quality at fixed g_c=1.0 ===")
sizes_all = sorted(all_data.keys())

nu_tests = {
    '5/6 (Potts exact)': 5/6,
    '1.0 (Ising exact)': 1.0,
    '0.63 (3D Ising)': 0.63,
    '0.67 (XY)': 0.6717,
}

results = {
    'experiment': '039c',
    'model': 'q=3 Potts data collapse',
    'sizes_available': sizes_all,
    'points_per_size': {str(n): len(all_data[n]) for n in sizes_all},
}

# Fixed g_c=1.0
print(f"\nSizes: {sizes_all}")
for label, nu in nu_tests.items():
    q = collapse_quality(nu, 1.0, sizes_all)
    print(f"  ν={label}: quality = {q:.6f}")

# Optimize ν with g_c=1.0 fixed
opt = minimize_scalar(lambda nu: collapse_quality(nu, 1.0, sizes_all),
                      bounds=(0.3, 3.0), method='bounded')
nu_opt_gc1 = opt.x
q_opt_gc1 = opt.fun
print(f"\n  Optimal ν (g_c=1 fixed): {nu_opt_gc1:.4f}, quality = {q_opt_gc1:.6f}")

# Quality ratio
q_potts = collapse_quality(5/6, 1.0, sizes_all)
q_ising = collapse_quality(1.0, 1.0, sizes_all)
print(f"\n  Potts/Ising quality ratio: {q_potts/q_ising:.3f} (< 1 means Potts wins)")

results['fixed_gc1'] = {
    'nu_opt': float(nu_opt_gc1),
    'quality_opt': float(q_opt_gc1),
    'quality_potts_5_6': float(q_potts),
    'quality_ising_1': float(q_ising),
    'potts_ising_ratio': float(q_potts / q_ising) if q_ising > 0 else float('inf'),
}

# ====== Joint optimization (ν and g_c free) ======
print("\n=== Joint optimization (ν and g_c both free) ===")

size_sets = {
    'all': sizes_all,
}
if 16 in all_data:
    size_sets['large (12+)'] = [n for n in sizes_all if n >= 12]
if 24 in all_data:
    size_sets['largest (16+)'] = [n for n in sizes_all if n >= 16]

results['joint_opt'] = {}

for label, sizes in size_sets.items():
    if len(sizes) < 2:
        continue

    # Start from ν=5/6, g_c=1.0
    res = minimize(lambda p: collapse_quality(p[0], p[1], sizes),
                   [5/6, 1.0], method='Nelder-Mead',
                   options={'xatol': 0.001, 'fatol': 1e-8})
    nu_j, gc_j = res.x
    q_j = res.fun

    # Also try starting from ν=1
    res2 = minimize(lambda p: collapse_quality(p[0], p[1], sizes),
                    [1.0, 1.0], method='Nelder-Mead',
                    options={'xatol': 0.001, 'fatol': 1e-8})
    if res2.fun < q_j:
        nu_j, gc_j = res2.x
        q_j = res2.fun

    q_potts_j = collapse_quality(5/6, gc_j, sizes)
    q_ising_j = collapse_quality(1.0, gc_j, sizes)

    print(f"\n{label} ({sizes}):")
    print(f"  Joint opt: ν={nu_j:.4f}, g_c={gc_j:.4f}, quality={q_j:.6f}")
    print(f"  At g_c={gc_j:.4f}: Potts(5/6)={q_potts_j:.6f}, Ising(1)={q_ising_j:.6f}")
    print(f"  Potts/Ising ratio: {q_potts_j/q_ising_j:.3f}")

    results['joint_opt'][label] = {
        'sizes': sizes,
        'nu_opt': float(nu_j),
        'gc_opt': float(gc_j),
        'quality_opt': float(q_j),
        'quality_potts': float(q_potts_j),
        'quality_ising': float(q_ising_j),
        'ratio': float(q_potts_j / q_ising_j) if q_ising_j > 0 else float('inf'),
    }

# ====== Pairwise crossing points ======
print("\n=== Crossing point estimates (g where CV(n1) = CV(n2)) ===")
results['crossings'] = {}

for i in range(len(sizes_all)):
    for j in range(i+1, len(sizes_all)):
        n1, n2 = sizes_all[i], sizes_all[j]
        # Find common g values where curves cross
        g_common = sorted(set(all_data[n1].keys()) & set(all_data[n2].keys()))
        if len(g_common) < 2:
            continue

        # Interpolate and find crossing
        g_arr = np.array(g_common)
        cv1 = np.array([all_data[n1][g] for g in g_common])
        cv2 = np.array([all_data[n2][g] for g in g_common])
        diff = cv1 - cv2  # positive when n1 > n2 (ordered side)

        for k in range(len(diff) - 1):
            if diff[k] > 0 and diff[k+1] < 0:
                # Linear interpolation for crossing
                g_cross = g_arr[k] + diff[k] * (g_arr[k+1] - g_arr[k]) / (diff[k] - diff[k+1])
                label = f"({n1},{n2})"
                print(f"  {label}: g_c ≈ {g_cross:.4f}")
                results['crossings'][label] = float(g_cross)

# ====== ν from crossing points (compare with collapse) ======
if len(results['crossings']) >= 2:
    print("\n=== ν from crossing-point pairs ===")
    cross_items = sorted(results['crossings'].items())
    results['nu_from_crossings'] = {}
    for i in range(len(cross_items)):
        for j in range(i+1, len(cross_items)):
            label1, gc1 = cross_items[i]
            label2, gc2 = cross_items[j]
            # Extract sizes from labels
            n1a = int(label1.split(',')[0].strip('('))
            n1b = int(label1.split(',')[1].strip(')'))
            n2a = int(label2.split(',')[0].strip('('))
            n2b = int(label2.split(',')[1].strip(')'))
            # Use geometric mean of pair sizes
            L1 = np.sqrt(n1a * n1b)
            L2 = np.sqrt(n2a * n2b)
            if abs(gc1 - gc2) > 1e-6 and L1 != L2:
                # g_c(L) = g_c_inf + a * L^{-1/nu}
                # Two-point: log(g_c1 - g_c2) / log(L2/L1) ~ -1/nu...
                # This is approximate. Just report crossing points.
                pass

# ====== Slope at transition ======
print("\n=== Transition slope at g=1.0 ===")
results['slopes'] = {}
for n in sorted(all_data.keys()):
    g_arr = np.array(sorted(all_data[n].keys()))
    cv_arr = np.array([all_data[n][g] for g in sorted(all_data[n].keys())])
    below = g_arr[g_arr < 1.0]
    above = g_arr[g_arr > 1.0]
    if len(below) > 0 and len(above) > 0:
        g_lo = below[-1]; g_hi = above[0]
        cv_lo = all_data[n][g_lo]; cv_hi = all_data[n][g_hi]
        slope = (cv_hi - cv_lo) / (g_hi - g_lo)
        print(f"  n={n:3d}: slope = {slope:.2f}")
        results['slopes'][str(n)] = float(slope)

# Slope scaling exponent
if len(results['slopes']) >= 2:
    ns = np.array([int(k) for k in sorted(results['slopes'].keys())])
    slopes = np.array([results['slopes'][str(n)] for n in ns])
    if len(ns) >= 2:
        log_n = np.log(ns)
        log_s = np.log(slopes)
        alpha = np.polyfit(log_n, log_s, 1)[0]
        print(f"  Slope ~ n^{alpha:.2f}")
        results['slope_exponent'] = float(alpha)

results['total_runtime'] = time.time() - t0
with open('results/sprint_039c_potts_collapse.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nTotal runtime: {time.time()-t0:.1f}s")
