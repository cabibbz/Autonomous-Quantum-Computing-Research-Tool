#!/usr/bin/env python3
"""Sprint 073b: q=2 Ly=4 cylinder gap crossings.

Ly=4 cylinder for q=2: dim = 2^{4*Lx} = 16^Lx.
Lx=3: dim=4096, Lx=4: 65536, Lx=5: 1048576, Lx=6: 16777216.
All feasible (Lx=6 similar to q=2 Ly=3 Lx=7 which ran at 2M dim).
Test if Ly convergence to 2D accelerates (77.7% at Ly=3, predict ~85-90%).
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
print("Sprint 073b: q=2 Ly=4 cylinder gap×Lx crossings", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '073b_cylinder_ly4_q2',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': 2, 'Ly': 4,
}

# Timing tests
for test_Lx, test_dim in [(5, 1048576), (6, 16777216)]:
    print(f"\nTiming test: Lx={test_Lx} (dim={test_dim})...", flush=True)
    t0 = time.time()
    try:
        gap_t, E0_t = get_gap(test_Lx, 4, 2, 0.7)
        dt = time.time() - t0
        print(f"  {dt:.1f}s per point", flush=True)
        results[f'timing_Lx{test_Lx}'] = float(dt)
    except Exception as e:
        print(f"  Failed: {e}", flush=True)
        results[f'timing_Lx{test_Lx}'] = 'failed'

# Expected g_c: Ly=3 was 0.655, 2D is 0.771
# Ly=2→3 jump: 0.451→0.655 = +0.204
# Ly=3→4 jump should be smaller. Estimate ~0.72
g_range = np.linspace(0.45, 0.85, 41)

# Determine max Lx based on timing
max_Lx = 5
if isinstance(results.get('timing_Lx6'), float) and results['timing_Lx6'] < 100:
    max_Lx = 6

Lx_list = [3, 4, 5]
if max_Lx >= 6:
    Lx_list.append(6)

gap_data = {}
for Lx in Lx_list:
    n = Lx * 4
    dim = 2 ** n
    g_scan = g_range
    # Coarser grid for large dims
    if dim > 5000000:
        g_scan = np.linspace(0.50, 0.80, 25)
    t0 = time.time()
    data = []
    for gi, g in enumerate(g_scan):
        gap, E0 = get_gap(Lx, 4, 2, g)
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

# Full Ly convergence table for q=2
print(f"\n{'='*60}", flush=True)
print("Ly convergence for q=2", flush=True)
print(f"{'='*60}", flush=True)
gc_1d = 0.250
gc_ly2 = 0.451
gc_ly3 = 0.655
gc_2d = 0.771
gc_ly4 = results.get('gc_mean', None)

ly_data = [
    (1, gc_1d), (2, gc_ly2), (3, gc_ly3),
]
if gc_ly4:
    ly_data.append((4, gc_ly4))

for Ly, gc in ly_data:
    progress = (gc - gc_1d) / (gc_2d - gc_1d) * 100
    print(f"  Ly={Ly}: g_c = {gc:.4f}  ({progress:.1f}% to 2D)", flush=True)
print(f"  2D:  g_c = {gc_2d:.4f}", flush=True)

if gc_ly4:
    results['ly_convergence'] = {
        'gc_1d': gc_1d, 'gc_ly2': gc_ly2, 'gc_ly3': gc_ly3,
        'gc_ly4': gc_ly4, 'gc_2d': gc_2d,
        'progress_ly4': float((gc_ly4 - gc_1d) / (gc_2d - gc_1d)),
    }
    # Check if convergence is exponential: gc(Ly) = gc_2d - A*exp(-Ly/B)
    # From 3 points (Ly=2,3,4) fit A, B
    from scipy.optimize import curve_fit
    def exp_conv(Ly, A, B):
        return gc_2d - A * np.exp(-np.array(Ly) / B)
    try:
        Lys = [2, 3, 4]
        gcs = [gc_ly2, gc_ly3, gc_ly4]
        popt, pcov = curve_fit(exp_conv, Lys, gcs, p0=[0.5, 1.5])
        A_fit, B_fit = popt
        gc_ly5_pred = exp_conv(5, A_fit, B_fit)
        print(f"\nExponential fit: g_c(Ly) = {gc_2d:.3f} - {A_fit:.3f}*exp(-Ly/{B_fit:.2f})", flush=True)
        print(f"  Prediction: g_c(Ly=5) = {gc_ly5_pred:.4f}", flush=True)
        results['exp_fit'] = {'A': float(A_fit), 'B': float(B_fit),
                              'gc_ly5_pred': float(gc_ly5_pred)}
    except Exception as e:
        print(f"Exponential fit failed: {e}", flush=True)

# Save
with open("results/sprint_073b_cylinder_ly4_q2.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_073b_cylinder_ly4_q2.json", flush=True)

from db_utils import record
if gc_ly4:
    record(sprint=73, model='hybrid_cyl', q=2, n=4,
           quantity='gc_gap', value=gc_ly4,
           method='gapLx_crossing_Ly4',
           notes=f'Ly=4 cylinder, {len(crossings)} crossing pairs')

print(f"\n{'='*60}", flush=True)
print("SUMMARY", flush=True)
print(f"{'='*60}", flush=True)
print(f"g_c(q=2, Ly=4 cyl) = {gc_ly4}", flush=True)
print(f"Crossings: {len(crossings)} pairs", flush=True)
for c in crossings:
    print(f"  ({c['Lx_pair'][0]},{c['Lx_pair'][1]}): {c['g_c']:.4f}", flush=True)
