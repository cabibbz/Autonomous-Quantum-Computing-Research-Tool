#!/usr/bin/env python3
"""Sprint 072c: q=2 Ly=3 cylinder gap crossings.

Ly=3 cylinder for q=2: dim = 2^{3*Lx} = 8^Lx.
Lx=3: dim=512, Lx=4: 4096, Lx=5: 32768, Lx=6: 262144, Lx=7: 2097152.
All feasible. Compare g_c(Ly=3) with Ly=2 (0.451) and 2D torus (0.771).
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
print("Sprint 072c: q=2 Ly=3 cylinder gap×Lx crossings", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '072c_cylinder_ly3_q2',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': 2, 'Ly': 3,
}

# Timing test for Lx=7 (dim=2M)
print("\nTiming test: Lx=7 (dim=2097152)...", flush=True)
t0 = time.time()
gap7, E07 = get_gap(7, 3, 2, 0.6)
dt7 = time.time() - t0
print(f"  {dt7:.1f}s per point", flush=True)

# Try Lx=8 (dim=16M) if Lx=7 is fast
feasible_8 = False
if dt7 < 30:
    print("Timing test: Lx=8 (dim=16777216)...", flush=True)
    t0 = time.time()
    try:
        gap8, E08 = get_gap(8, 3, 2, 0.6)
        dt8 = time.time() - t0
        print(f"  {dt8:.1f}s per point", flush=True)
        if dt8 < 120:
            feasible_8 = True
    except Exception as e:
        print(f"  Failed: {e}", flush=True)

# Expected g_c: between Ly=2 (0.451) and 2D (0.771)
# Ly=2 z=3, Ly=3 z≈3.33 (most sites have 4 nbrs, edge sites 3), 2D z=4
# Linear interpolation in z: g_c ~ 0.45 + (0.77-0.45)*(3.33-3)/(4-3) ≈ 0.56
g_range = np.linspace(0.35, 0.85, 51)

Lx_list = [3, 4, 5, 6, 7]
if feasible_8:
    Lx_list.append(8)

gap_data = {}
for Lx in Lx_list:
    n = Lx * 3
    dim = 2 ** n
    t0 = time.time()
    data = []
    for gi, g in enumerate(g_range):
        gap, E0 = get_gap(Lx, 3, 2, g)
        data.append({'g': round(float(g), 4), 'gap': float(gap),
                     'gap_Lx': float(gap * Lx), 'E0': float(E0)})
        if dim > 100000 and (gi + 1) % 10 == 0:
            elapsed = time.time() - t0
            est_total = elapsed * len(g_range) / (gi + 1)
            print(f"  Lx={Lx}: {gi+1}/{len(g_range)} ({elapsed:.0f}s, est {est_total:.0f}s total)", flush=True)
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
    g_common = np.array([d['g'] for d in gap_data[Lx1]])
    g1 = np.array([d['gap_Lx'] for d in gap_data[Lx1]])
    g2_g = np.array([d['g'] for d in gap_data[Lx2]])
    g2_v = np.array([d['gap_Lx'] for d in gap_data[Lx2]])
    g2_interp = np.interp(g_common, g2_g, g2_v)
    diff = g1 - g2_interp
    found = False
    for j in range(len(diff) - 1):
        if diff[j] * diff[j + 1] < 0:
            g_cross = g_common[j] - diff[j] * (g_common[j + 1] - g_common[j]) / (diff[j + 1] - diff[j])
            gLx_at = np.interp(g_cross, g_common, g1)
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

# Compare Ly convergence
print(f"\n{'='*60}", flush=True)
print("Ly convergence for q=2", flush=True)
print(f"{'='*60}", flush=True)
gc_ly2 = 0.451
gc_2d = 0.771
gc_ly3 = results.get('gc_mean', None)

print(f"  Ly=2 (ladder): g_c = {gc_ly2:.3f}", flush=True)
if gc_ly3:
    print(f"  Ly=3:          g_c = {gc_ly3:.3f}", flush=True)
    print(f"  2D (L→∞):     g_c = {gc_2d:.3f}", flush=True)
    progress = (gc_ly3 - gc_ly2) / (gc_2d - gc_ly2)
    print(f"  Ly=3 is {progress*100:.1f}% of the way from Ly=2 to 2D", flush=True)
    results['ly_convergence'] = {
        'gc_ly2': gc_ly2, 'gc_ly3': gc_ly3, 'gc_2d': gc_2d,
        'progress_frac': float(progress),
    }

# Save
with open("results/sprint_072c_cylinder_ly3_q2.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_072c_cylinder_ly3_q2.json", flush=True)

from db_utils import record
if gc_ly3:
    record(sprint=72, model='hybrid_cyl', q=2, n=3,
           quantity='gc_gap', value=gc_ly3,
           method='gapLx_crossing_Ly3',
           notes=f'Ly=3 cylinder, {len(crossings)} crossing pairs')

print(f"\n{'='*60}", flush=True)
print("SUMMARY", flush=True)
print(f"{'='*60}", flush=True)
print(f"g_c(q=2, Ly=3 cyl) = {gc_ly3}", flush=True)
print(f"Crossings: {len(crossings)} pairs", flush=True)
for c in crossings:
    print(f"  ({c['Lx_pair'][0]},{c['Lx_pair'][1]}): {c['g_c']:.4f}", flush=True)
