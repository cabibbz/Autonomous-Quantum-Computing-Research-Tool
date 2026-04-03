#!/usr/bin/env python3
"""Sprint 073a: q=3 Ly=3 cylinder gap crossings.

Ly=3 cylinder for q=3: dim = 3^{3*Lx} = 27^Lx.
Lx=3: dim=19683, Lx=4: 531441, Lx=5: 14348907.
Lx=3,4 easy. Lx=5 feasible with GPU (~14M, similar to q=5 n=10).
Compare g_c(Ly=3) with Ly=2 (0.565) and 2D torus (1.267).
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
print("Sprint 073a: q=3 Ly=3 cylinder gap×Lx crossings", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '073a_cylinder_ly3_q3',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': 3, 'Ly': 3,
}

# Timing test
print("\nTiming test: Lx=4 (dim=531441)...", flush=True)
t0 = time.time()
gap4, E04 = get_gap(4, 3, 3, 0.8)
dt4 = time.time() - t0
print(f"  Lx=4: {dt4:.1f}s per point", flush=True)

print("Timing test: Lx=5 (dim=14348907)...", flush=True)
t0 = time.time()
try:
    gap5, E05 = get_gap(5, 3, 3, 0.8)
    dt5 = time.time() - t0
    print(f"  Lx=5: {dt5:.1f}s per point", flush=True)
    feasible_5 = dt5 < 120  # need ~30 g-points
except Exception as e:
    print(f"  Lx=5 failed: {e}", flush=True)
    dt5 = None
    feasible_5 = False

# g_c expected between Ly=2 (0.565) and 2D (1.267)
# Ly=2→3 jump for q=2 was 0.451→0.655 (45% increase)
# For q=3: 0.565 * 1.45 ≈ 0.82
# Scan range centered on that estimate
g_range = np.linspace(0.4, 1.2, 41)

Lx_list = [3, 4]
if feasible_5:
    Lx_list.append(5)
    # Use coarser grid for Lx=5
    g_range_5 = np.linspace(0.5, 1.1, 25)

gap_data = {}
for Lx in Lx_list:
    n = Lx * 3
    dim = 3 ** n
    g_scan = g_range_5 if (Lx == 5 and feasible_5) else g_range
    t0 = time.time()
    data = []
    for gi, g in enumerate(g_scan):
        gap, E0 = get_gap(Lx, 3, 3, g)
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

# Find crossings
print("\n--- Crossings ---", flush=True)
crossings = []
Lx_sorted = sorted(gap_data.keys())
for i in range(len(Lx_sorted) - 1):
    Lx1, Lx2 = Lx_sorted[i], Lx_sorted[i + 1]
    g1_g = np.array([d['g'] for d in gap_data[Lx1]])
    g1_v = np.array([d['gap_Lx'] for d in gap_data[Lx1]])
    g2_g = np.array([d['g'] for d in gap_data[Lx2]])
    g2_v = np.array([d['gap_Lx'] for d in gap_data[Lx2]])
    # Interpolate both to common grid
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

results['crossings'] = crossings

if crossings:
    gc_vals = [c['g_c'] for c in crossings]
    gc_mean = np.mean(gc_vals)
    results['gc_mean'] = float(gc_mean)
    results['gc_spread'] = float(max(gc_vals) - min(gc_vals)) if len(gc_vals) > 1 else 0.0

# Ly convergence for q=3
print(f"\n{'='*60}", flush=True)
print("Ly convergence for q=3", flush=True)
print(f"{'='*60}", flush=True)
gc_1d = 0.333
gc_ly2 = 0.565
gc_2d = 1.267
gc_ly3 = results.get('gc_mean', None)

print(f"  1D:    g_c = {gc_1d:.3f}", flush=True)
print(f"  Ly=2:  g_c = {gc_ly2:.3f}  ({(gc_ly2-gc_1d)/(gc_2d-gc_1d)*100:.1f}% to 2D)", flush=True)
if gc_ly3:
    print(f"  Ly=3:  g_c = {gc_ly3:.3f}  ({(gc_ly3-gc_1d)/(gc_2d-gc_1d)*100:.1f}% to 2D)", flush=True)
    print(f"  2D:    g_c = {gc_2d:.3f}", flush=True)
    progress = (gc_ly3 - gc_1d) / (gc_2d - gc_1d)
    results['ly_convergence'] = {
        'gc_1d': gc_1d, 'gc_ly2': gc_ly2, 'gc_ly3': gc_ly3, 'gc_2d': gc_2d,
        'progress_frac': float(progress),
    }

# Save
with open("results/sprint_073a_cylinder_ly3_q3.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_073a_cylinder_ly3_q3.json", flush=True)

from db_utils import record
if gc_ly3:
    record(sprint=73, model='hybrid_cyl', q=3, n=3,
           quantity='gc_gap', value=gc_ly3,
           method='gapLx_crossing_Ly3',
           notes=f'Ly=3 cylinder, {len(crossings)} crossing pairs')

print(f"\n{'='*60}", flush=True)
print("SUMMARY", flush=True)
print(f"{'='*60}", flush=True)
print(f"g_c(q=3, Ly=3 cyl) = {gc_ly3}", flush=True)
print(f"Crossings: {len(crossings)} pairs", flush=True)
for c in crossings:
    print(f"  ({c['Lx_pair'][0]},{c['Lx_pair'][1]}): {c['g_c']:.4f}", flush=True)
