#!/usr/bin/env python3
"""Sprint 075a: q=2 Ly=5 cylinder gap crossings.

Ly=5 cylinder for q=2: dim = 2^{5*Lx} = 32^Lx.
Lx=3: dim=32768 (trivial).
Lx=4: dim=1,048,576 (~1M, feasible with GPU).
Lx=5: dim=33,554,432 (~33M, feasible with GPU if <300s/pt).
Lx=6: dim=1,073,741,824 (~1B, infeasible).

Prediction from exponential fit: g_c(Ly=5) ≈ 0.745 (95% to 2D).
Known: g_c(Ly=4) = 0.688, g_c(2D) = 0.771.
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
    # Potts coupling: delta(s_i, s_j)
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

    # Build bond list: open x, periodic y
    bonds = set()
    for y in range(Ly):
        for x in range(Lx):
            site = y * Lx + x
            # Horizontal bond (open BC in x)
            if x + 1 < Lx:
                nbr = y * Lx + (x + 1)
                bonds.add((min(site, nbr), max(site, nbr)))
            # Vertical bond (periodic BC in y)
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
            # Non-adjacent bond (vertical wrap-around)
            indices = np.arange(dim)
            si = (indices // (q ** i)) % q
            sj = (indices // (q ** j)) % q
            diag_vals = (si == sj).astype(float)
            op = diags(diag_vals, 0, shape=(dim, dim), format='csr')
        H = H - op

    # Transverse field
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
print("Sprint 075a: q=2 Ly=5 cylinder gap×Lx crossings", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '075a_cylinder_ly5_q2',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': 2, 'Ly': 5,
}

# Timing tests
for Lx_test in [3, 4, 5]:
    n = Lx_test * 5
    dim = 2 ** n
    print(f"\nTiming test: Lx={Lx_test} (n={n}, dim={dim:,})...", flush=True)
    t0 = time.time()
    gap, E0 = get_gap(Lx_test, 5, 2, 0.7)
    dt = time.time() - t0
    print(f"  Lx={Lx_test}: {dt:.1f}s per point", flush=True)
    results[f'timing_Lx{Lx_test}'] = float(dt)
    if dt > 250:
        print(f"  TOO SLOW — skipping Lx={Lx_test} in full scan", flush=True)
        break

# Determine feasible Lx values
Lx_list = []
for Lx in [3, 4, 5]:
    key = f'timing_Lx{Lx}'
    if key in results and results[key] < 250:
        Lx_list.append(Lx)
    elif key not in results:
        break

print(f"\nFeasible Lx values: {Lx_list}", flush=True)

# Scan range: between Ly=4 (0.688) and 2D (0.771)
# Wider range to be safe
g_range = np.linspace(0.5, 0.9, 40)

gap_data = {}
for Lx in Lx_list:
    n = Lx * 5
    dim = 2 ** n
    t0 = time.time()
    data = []
    for gi, g in enumerate(g_range):
        gap, E0 = get_gap(Lx, 5, 2, g)
        data.append({'g': round(float(g), 4), 'gap': float(gap),
                     'gap_Lx': float(gap * Lx), 'E0': float(E0)})
        if dim > 100000 and (gi + 1) % 5 == 0:
            elapsed = time.time() - t0
            est_total = elapsed * len(g_range) / (gi + 1)
            print(f"  Lx={Lx}: {gi+1}/{len(g_range)} ({elapsed:.0f}s, est {est_total:.0f}s total)", flush=True)
    dt = time.time() - t0
    gapLx = [d['gap_Lx'] for d in data]
    print(f"Lx={Lx} (n={n}, dim={dim:,}): {dt:.1f}s total, gap×Lx=[{min(gapLx):.4f}, {max(gapLx):.4f}]", flush=True)
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
    g_common = np.linspace(max(g1_g[0], g2_g[0]), min(g1_g[-1], g2_g[-1]), 500)
    v1 = np.interp(g_common, g1_g, g1_v)
    v2 = np.interp(g_common, g2_g, g2_v)
    diff = v1 - v2
    found = False
    for j in range(len(diff) - 1):
        if diff[j] * diff[j + 1] < 0:
            g_cross = g_common[j] - diff[j] * (g_common[j + 1] - g_common[j]) / (diff[j + 1] - diff[j])
            gLx_at = float(np.interp(g_cross, g1_g, g1_v))
            print(f"  ({Lx1},{Lx2}): g_c = {g_cross:.4f}, gap×Lx = {gLx_at:.4f}", flush=True)
            crossings.append({
                'Lx_pair': [Lx1, Lx2], 'g_c': float(g_cross),
                'gap_Lx_at_crossing': gLx_at
            })
            found = True
            break
    if not found:
        print(f"  ({Lx1},{Lx2}): NO crossing found", flush=True)

results['crossings'] = crossings

if crossings:
    gc_vals = [c['g_c'] for c in crossings]
    gc_mean = np.mean(gc_vals)
    results['gc_mean'] = float(gc_mean)
    print(f"\nMean g_c = {gc_mean:.4f}", flush=True)

# Ly convergence for q=2
print(f"\n{'='*60}", flush=True)
print("Ly convergence for q=2 (complete series)", flush=True)
print(f"{'='*60}", flush=True)
gc_series = {
    1: 0.250,
    2: 0.451,
    3: 0.655,
    4: 0.688,
}
gc_2d = 0.771
gc_ly5 = results.get('gc_mean', None)
if gc_ly5:
    gc_series[5] = gc_ly5

for Ly, gc in sorted(gc_series.items()):
    pct = (gc - 0.250) / (gc_2d - 0.250) * 100
    print(f"  Ly={Ly}: g_c = {gc:.3f}  ({pct:.1f}% to 2D)", flush=True)
print(f"  2D:   g_c = {gc_2d:.3f}  (100%)", flush=True)

if gc_ly5:
    pct5 = (gc_ly5 - 0.250) / (gc_2d - 0.250) * 100
    results['progress_to_2D'] = float(pct5)
    # Check prediction
    predicted = 0.745
    error = abs(gc_ly5 - predicted) / predicted * 100
    print(f"\n  Prediction: {predicted:.3f}, Measured: {gc_ly5:.3f}, Error: {error:.1f}%", flush=True)
    results['prediction_error_pct'] = float(error)

# Save
with open("results/sprint_075a_cylinder_ly5_q2.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_075a_cylinder_ly5_q2.json", flush=True)

from db_utils import record
if gc_ly5:
    record(sprint=75, model='hybrid_cyl', q=2, n=5,
           quantity='gc_gap', value=gc_ly5,
           method=f'gapLx_crossing_Ly5',
           notes=f'Ly=5 cylinder, {len(crossings)} crossings')

print(f"\n{'='*60}", flush=True)
print("SUMMARY", flush=True)
print(f"{'='*60}", flush=True)
if gc_ly5:
    print(f"g_c(q=2, Ly=5 cyl) = {gc_ly5:.4f}", flush=True)
    print(f"Crossings: {len(crossings)}", flush=True)
    for c in crossings:
        print(f"  ({c['Lx_pair'][0]},{c['Lx_pair'][1]}): g_c={c['g_c']:.4f}", flush=True)
else:
    print("No crossings found — check data", flush=True)
