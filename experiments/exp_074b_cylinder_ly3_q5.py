#!/usr/bin/env python3
"""Sprint 074b: q=5 Ly=3 cylinder gap crossings.

Ly=3 cylinder for q=5: dim = 5^{3*Lx} = 125^Lx.
Lx=3: dim=1953125 (~2M, feasible with GPU).
Lx=4: dim=244140625 (~244M, likely infeasible).
Try Lx=3 and Lx=4 timing. If only Lx=3 available, use gap minimum
and compare with Ly=2 gap shape to estimate g_c.
"""
import numpy as np
import json
import time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from gpu_utils import eigsh


def build_H_cylinder(Lx, Ly, q, g):
    """Build hybrid Hamiltonian on Lx x Ly cylinder (open x, periodic y)."""
    n = Lx * Ly
    dim = q ** n
    potts_2site = np.zeros(q ** 2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q ** 2, q ** 2), format='csr')
    X = np.zeros((q, q))
    for s in range(q):
        X[(s + 1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)
    H = csr_matrix((dim, dim))
    bonds = set()
    for y in range(Ly):
        for x in range(Lx):
            site = y * Lx + x
            if x + 1 < Lx:
                nbr = y * Lx + (x + 1)
                bonds.add((min(site, nbr), max(site, nbr)))
            nbr_v = ((y + 1) % Ly) * Lx + x
            if nbr_v != site:
                bonds.add((min(site, nbr_v), max(site, nbr_v)))
    for (i, j) in bonds:
        if j == i + 1:
            left = q ** i
            right = q ** (n - j - 1)
            op = potts_op
            if left > 1:
                op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1:
                op = sp_kron(op, sp_eye(right), format='csr')
        else:
            indices = np.arange(dim)
            si = (indices // (q ** i)) % q
            sj = (indices // (q ** j)) % q
            diag_vals = (si == sj).astype(float)
            op = diags(diag_vals, 0, shape=(dim, dim), format='csr')
        H = H - op
    for i in range(n):
        left = q ** i
        right = q ** (n - i - 1)
        op = XpXd.copy()
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H


def get_gap(Lx, Ly, q, g):
    """Get energy gap."""
    H = build_H_cylinder(Lx, Ly, q, g)
    dim = H.shape[0]
    if dim < 500:
        evals = np.linalg.eigvalsh(H.toarray())
        return evals[1] - evals[0], evals[0]
    k = min(6, dim - 1)
    evals, _ = eigsh(H, k=k, which='SA')
    evals = np.sort(evals)
    return evals[1] - evals[0], evals[0]


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 074b: q=5 Ly=3 cylinder gap×Lx crossings", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '074b_cylinder_ly3_q5',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': 5, 'Ly': 3,
}

# Timing tests
print("\nTiming test: Lx=2 (dim=15625)...", flush=True)
t0 = time.time()
gap2, E02 = get_gap(2, 3, 5, 0.9)
dt2 = time.time() - t0
print(f"  Lx=2: {dt2:.1f}s per point", flush=True)
results['timing_Lx2'] = float(dt2)

print("Timing test: Lx=3 (dim=1953125)...", flush=True)
t0 = time.time()
gap3, E03 = get_gap(3, 3, 5, 0.9)
dt3 = time.time() - t0
print(f"  Lx=3: {dt3:.1f}s per point", flush=True)
results['timing_Lx3'] = float(dt3)

# Lx=4 at 244M is certainly infeasible (q=3 Ly=4 Lx=4 at 43M took 838s)
print("Lx=4 (dim=244M) — skipped, infeasible based on 074a timing", flush=True)
results['timing_Lx4'] = 'skipped_infeasible'

# Scan range: between Ly=2 (0.714) and 2D (1.588)
# Estimate: q=5 Ly=2 is 0.714; estimate Ly=3 ~ 0.714 * 1.35 ≈ 0.96
g_range = np.linspace(0.4, 1.5, 45)

# Include Lx=2 for crossing with Lx=3
Lx_list = [2, 3]

gap_data = {}
for Lx in Lx_list:
    n = Lx * 3
    dim = 5 ** n
    g_scan = g_range
    t0 = time.time()
    data = []
    for gi, g in enumerate(g_scan):
        gap, E0 = get_gap(Lx, 3, 5, g)
        data.append({'g': round(float(g), 4), 'gap': float(gap),
                     'gap_Lx': float(gap * Lx), 'E0': float(E0)})
        if dim > 100000 and (gi + 1) % 5 == 0:
            elapsed = time.time() - t0
            est_total = elapsed * len(g_scan) / (gi + 1)
            print(f"  Lx={Lx}: {gi+1}/{len(g_scan)} ({elapsed:.0f}s, est {est_total:.0f}s total)", flush=True)
    dt = time.time() - t0
    gapLx = [d['gap_Lx'] for d in data]
    print(f"Lx={Lx} (n={n}, dim={dim}): {dt:.1f}s, gap×Lx=[{min(gapLx):.4f}, {max(gapLx):.4f}]", flush=True)
    gap_data[Lx] = data

results['gap_data'] = {str(k): v for k, v in gap_data.items()}

# Find crossings if 2+ sizes
print("\n--- Crossings ---", flush=True)
crossings = []
Lx_sorted = sorted(gap_data.keys())
if len(Lx_sorted) >= 2:
    for i in range(len(Lx_sorted) - 1):
        Lx1, Lx2 = Lx_sorted[i], Lx_sorted[i + 1]
        g1_g = np.array([d['g'] for d in gap_data[Lx1]])
        g1_v = np.array([d['gap_Lx'] for d in gap_data[Lx1]])
        g2_g = np.array([d['g'] for d in gap_data[Lx2]])
        g2_v = np.array([d['gap_Lx'] for d in gap_data[Lx2]])
        g_common = np.linspace(max(g1_g[0], g2_g[0]), min(g1_g[-1], g2_g[-1]), 200)
        v1 = np.interp(g_common, g1_g, g1_v)
        v2 = np.interp(g_common, g2_g, g2_v)
        diff = v1 - v2
        found = False
        for j in range(len(diff) - 1):
            if diff[j] * diff[j + 1] < 0:
                g_cross = g_common[j] - diff[j] * (g_common[j + 1] - g_common[j]) / (diff[j + 1] - diff[j])
                gLx_at = np.interp(g_cross, g1_g, g1_v)
                print(f"  ({Lx1},{Lx2}): g_c = {g_cross:.4f}, gap×Lx = {gLx_at:.4f}", flush=True)
                crossings.append({
                    'Lx_pair': [Lx1, Lx2], 'g_c': float(g_cross),
                    'gap_Lx_at_crossing': float(gLx_at)
                })
                found = True
                break
        if not found:
            print(f"  ({Lx1},{Lx2}): NO crossing", flush=True)
else:
    print("  Only one Lx available — using gap minimum as proxy", flush=True)

results['crossings'] = crossings

if crossings:
    gc_vals = [c['g_c'] for c in crossings]
    gc_mean = np.mean(gc_vals)
    results['gc_mean'] = float(gc_mean)

# Gap minimum analysis (works even with single Lx)
gaps = [d['gap'] for d in gap_data[Lx_list[0]]]
g_vals = [d['g'] for d in gap_data[Lx_list[0]]]
min_idx = np.argmin(gaps)
g_min = g_vals[min_idx]
gap_min = gaps[min_idx]
print(f"\nGap minimum at g = {g_min:.4f}, gap = {gap_min:.6f} (Lx={Lx_list[0]})", flush=True)
results['gap_min_proxy'] = {'g': float(g_min), 'gap': float(gap_min), 'Lx': Lx_list[0]}

# Also compute d²E/dg² to find peak (another g_c proxy)
E0_vals = np.array([d['E0'] for d in gap_data[Lx_list[0]]])
g_arr = np.array([d['g'] for d in gap_data[Lx_list[0]]])
dg = g_arr[1] - g_arr[0]
d2E = np.gradient(np.gradient(E0_vals, dg), dg)
peak_idx = np.argmin(d2E)  # most negative d²E/dg²
g_peak = g_arr[peak_idx]
print(f"d²E/dg² peak at g = {g_peak:.4f} (Lx={Lx_list[0]})", flush=True)
results['d2E_peak'] = {'g': float(g_peak), 'Lx': Lx_list[0]}

# Ly convergence for q=5
print(f"\n{'='*60}", flush=True)
print("Ly convergence for q=5", flush=True)
print(f"{'='*60}", flush=True)
gc_1d = 0.441
gc_ly2 = 0.714
gc_2d = 1.588

gc_ly3 = results.get('gc_mean', g_min)  # Use crossing if available, else gap min
method = 'crossing' if crossings else 'gap_minimum'

print(f"  1D:    g_c = {gc_1d:.3f}", flush=True)
print(f"  Ly=2:  g_c = {gc_ly2:.3f}  ({(gc_ly2-gc_1d)/(gc_2d-gc_1d)*100:.1f}% to 2D)", flush=True)
print(f"  Ly=3:  g_c ~ {gc_ly3:.3f}  ({(gc_ly3-gc_1d)/(gc_2d-gc_1d)*100:.1f}% to 2D, method={method})", flush=True)
print(f"  2D:    g_c = {gc_2d:.3f}", flush=True)

results['ly_convergence'] = {
    'gc_1d': gc_1d, 'gc_ly2': gc_ly2, 'gc_ly3_est': gc_ly3,
    'gc_2d': gc_2d, 'method': method,
    'progress_frac': float((gc_ly3 - gc_1d) / (gc_2d - gc_1d)),
}

# Save
with open("results/sprint_074b_cylinder_ly3_q5.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_074b_cylinder_ly3_q5.json", flush=True)

from db_utils import record
record(sprint=74, model='hybrid_cyl', q=5, n=3,
       quantity='gc_gap', value=gc_ly3,
       method=f'gapLx_{method}_Ly3',
       notes=f'Ly=3 cylinder, {len(crossings)} crossings, method={method}')

print(f"\n{'='*60}", flush=True)
print("SUMMARY", flush=True)
print(f"{'='*60}", flush=True)
print(f"g_c(q=5, Ly=3 cyl) ~ {gc_ly3:.4f} (method={method})", flush=True)
print(f"Gap minimum: g={g_min:.4f}", flush=True)
print(f"d²E peak: g={g_peak:.4f}", flush=True)
