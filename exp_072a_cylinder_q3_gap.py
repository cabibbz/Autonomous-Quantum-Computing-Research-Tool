#!/usr/bin/env python3
"""Sprint 072a: q=3 Ly=2 cylinder gap×Lx crossings.

q=3 is exactly solvable in 1D (g_c=1/3). Fill the gap in cylinder data
between q=2 (g_c_cyl=0.451) and q=5 (g_c_cyl=0.714).

Sizes: Lx=3,4,5,6 (dim = 3^6=729, 3^8=6561, 3^10=59049, 3^12=531441).
All feasible on CPU/GPU within time limits.
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

    # Potts coupling: delta(s_i, s_j) diagonal for adjacent sites
    potts_2site = np.zeros(q ** 2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q ** 2, q ** 2), format='csr')

    # Clock field: X + X†
    X = np.zeros((q, q))
    for s in range(q):
        X[(s + 1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)

    H = csr_matrix((dim, dim))

    # Build bond list: x-direction open, y-direction periodic
    bonds = set()
    for y in range(Ly):
        for x in range(Lx):
            site = y * Lx + x
            # x-direction (open)
            if x + 1 < Lx:
                nbr = y * Lx + (x + 1)
                bonds.add((min(site, nbr), max(site, nbr)))
            # y-direction (periodic)
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
print("Sprint 072a: q=3 Ly=2 cylinder gap×Lx crossings", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '072a_cylinder_q3_gap',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': 3, 'Ly': 2,
}

# Timing test first
print("\nTiming test: Lx=6 (dim=531441) single point...", flush=True)
t0 = time.time()
gap_test, E0_test = get_gap(6, 2, 3, 0.6)
dt_test = time.time() - t0
print(f"  {dt_test:.1f}s per point. gap={gap_test:.6f}", flush=True)

# Determine if Lx=7 (dim=3^14=4.8M) is feasible
feasible_7 = False
if dt_test < 5:
    print("\nTiming test: Lx=7 (dim=4782969) single point...", flush=True)
    t0 = time.time()
    try:
        gap7, E07 = get_gap(7, 2, 3, 0.6)
        dt7 = time.time() - t0
        print(f"  {dt7:.1f}s per point. gap={gap7:.6f}", flush=True)
        if dt7 < 60:
            feasible_7 = True
    except Exception as e:
        print(f"  Failed: {e}", flush=True)

# Scan g range — centered around expected g_c
# 1D g_c = 0.333, cyl/1D for q=2 is 1.80 → predict g_c(cyl,q=3) ≈ 0.60
# cyl/1D for q=5 is 1.62 → q=3 might be ~1.7 → predict 0.57
g_range = np.linspace(0.3, 0.9, 61)

Lx_list = [3, 4, 5, 6]
if feasible_7:
    Lx_list.append(7)

gap_data = {}
for Lx in Lx_list:
    n = Lx * 2
    dim = 3 ** n
    t0 = time.time()
    data = []
    for gi, g in enumerate(g_range):
        gap, E0 = get_gap(Lx, 2, 3, g)
        data.append({'g': round(float(g), 4), 'gap': float(gap),
                     'gap_Lx': float(gap * Lx), 'E0': float(E0)})
        if dim > 100000 and (gi + 1) % 15 == 0:
            print(f"  Lx={Lx}: {gi+1}/{len(g_range)} done", flush=True)
    dt = time.time() - t0
    gapLx = [d['gap_Lx'] for d in data]
    print(f"Lx={Lx} (dim={dim}): {dt:.1f}s, gap×Lx=[{min(gapLx):.4f}, {max(gapLx):.4f}]", flush=True)
    gap_data[Lx] = data

results['gap_data'] = {str(k): v for k, v in gap_data.items()}

# Find crossings
print("\n--- Crossings ---", flush=True)
crossings = []
Lx_sorted = sorted(gap_data.keys())
for i in range(len(Lx_sorted) - 1):
    Lx1, Lx2 = Lx_sorted[i], Lx_sorted[i + 1]
    g1 = np.array([d['gap_Lx'] for d in gap_data[Lx1]])
    g2_raw = gap_data[Lx2]
    g2_g = np.array([d['g'] for d in g2_raw])
    g2_v = np.array([d['gap_Lx'] for d in g2_raw])
    # Interpolate to common grid
    g_common = np.array([d['g'] for d in gap_data[Lx1]])
    g2_interp = np.interp(g_common, g2_g, g2_v)
    diff = g1 - g2_interp
    found = False
    for j in range(len(diff) - 1):
        if diff[j] * diff[j + 1] < 0:
            g_cross = g_common[j] - diff[j] * (g_common[j + 1] - g_common[j]) / (diff[j + 1] - diff[j])
            gLx_at_cross = np.interp(g_cross, g_common, g1)
            print(f"  ({Lx1},{Lx2}): g_c = {g_cross:.4f}, gap×Lx = {gLx_at_cross:.4f}", flush=True)
            crossings.append({
                'Lx_pair': [Lx1, Lx2], 'g_c': float(g_cross),
                'gap_Lx_at_crossing': float(gLx_at_cross)
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
    cyl_1d_ratio = gc_mean / (1.0 / 3.0)
    results['cyl_1d_ratio'] = float(cyl_1d_ratio)
    print(f"\n  Mean g_c = {gc_mean:.4f}", flush=True)
    print(f"  Cyl/1D ratio = {cyl_1d_ratio:.3f}  [q=2: 1.80, q=5: 1.62]", flush=True)

# Save
with open("results/sprint_072a_cylinder_q3_gap.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_072a_cylinder_q3_gap.json", flush=True)

from db_utils import record
if 'gc_mean' in results:
    record(sprint=72, model='hybrid_cyl', q=3, n=2,
           quantity='gc_gap', value=results['gc_mean'],
           method='gapLx_crossing_Ly2',
           notes=f"Ly=2 cylinder, {len(crossings)} crossing pairs")

print(f"\n{'='*60}", flush=True)
print("SUMMARY", flush=True)
print(f"{'='*60}", flush=True)
print(f"g_c(q=3, Ly=2 cyl) = {results.get('gc_mean', 'N/A')}", flush=True)
print(f"1D g_c = 0.333 (exact)", flush=True)
print(f"Cyl/1D = {results.get('cyl_1d_ratio', 'N/A')}", flush=True)
print(f"Crossings: {len(crossings)} pairs", flush=True)
for c in crossings:
    print(f"  ({c['Lx_pair'][0]},{c['Lx_pair'][1]}): {c['g_c']:.4f}", flush=True)
