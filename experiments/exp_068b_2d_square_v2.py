#!/usr/bin/env python3
"""Sprint 068b (v2): 2D phase transition on square lattices — q=2,3,5

Use gap*L scaling on square L×L periodic lattices.
GPU-accelerated for large Hilbert spaces. Fewer g points for big systems.

Lattices:
  q=2: L=2 (16), L=3 (512), L=4 (65536) — all CPU
  q=3: L=2 (81), L=3 (19683) — CPU
  q=5: L=2 (625), L=3 (1953125) — L=3 needs GPU, fewer points
"""
import numpy as np
import json
import time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh

# Check GPU availability
try:
    import cupy as cp
    from cupyx.scipy.sparse import csr_matrix as cp_csr
    from cupyx.scipy.sparse.linalg import eigsh as cp_eigsh
    HAS_GPU = True
    print("GPU available: NVIDIA TITAN RTX")
except ImportError:
    HAS_GPU = False
    print("GPU not available, using CPU only")

GPU_THRESHOLD = 100000  # Use GPU above this dim


def build_H_2d(Lx, Ly, q, g):
    """Build 2D hybrid Hamiltonian on Lx x Ly periodic square lattice."""
    n = Lx * Ly
    dim = q**n

    potts_2site = np.zeros(q**2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q**2, q**2), format='csr')

    X = np.zeros((q, q))
    for s in range(q):
        X[(s + 1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)

    H = csr_matrix((dim, dim))

    # Enumerate bonds (periodic square lattice, no duplicates)
    bonds = []
    for y in range(Ly):
        for x in range(Lx):
            site = y * Lx + x
            nbr_h = y * Lx + (x + 1) % Lx
            if nbr_h != site:
                bonds.append((min(site, nbr_h), max(site, nbr_h)))
            nbr_v = ((y + 1) % Ly) * Lx + x
            if nbr_v != site:
                bonds.append((min(site, nbr_v), max(site, nbr_v)))
    bonds = list(set(bonds))

    for (i, j) in bonds:
        if j == i + 1:
            left = q**i
            right = q**(n - i - 2)
            op = potts_op
            if left > 1:
                op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1:
                op = sp_kron(op, sp_eye(right), format='csr')
        else:
            diag_vals = np.zeros(dim)
            for idx in range(dim):
                si = (idx // (q**i)) % q
                sj = (idx // (q**j)) % q
                diag_vals[idx] = 1.0 if si == sj else 0.0
            op = diags(diag_vals, 0, shape=(dim, dim), format='csr')
        H = H - op

    for i in range(n):
        left = q**i
        right = q**(n - i - 1)
        op = XpXd.copy()
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op

    return H


def eigsolve(H, k=4):
    """Solve for lowest k eigenvalues, using GPU if beneficial."""
    dim = H.shape[0]
    if dim < 500:
        evals = np.linalg.eigvalsh(H.toarray())[:k]
        return evals

    if HAS_GPU and dim > GPU_THRESHOLD:
        try:
            H_gpu = cp_csr(H)
            evals_gpu, _ = cp_eigsh(H_gpu, k=k, which='SA')
            return np.sort(cp.asnumpy(evals_gpu))
        except Exception as e:
            print(f"    GPU failed ({e}), falling back to CPU")

    evals, _ = eigsh(H, k=k, which='SA')
    return np.sort(evals)


def scan_gap_square(L, q, g_values):
    """Scan energy gap for L×L lattice."""
    n = L * L
    data = []
    for gi, g in enumerate(g_values):
        t0 = time.time()
        H = build_H_2d(L, L, q, g)
        evals = eigsolve(H, k=4)
        gap = evals[1] - evals[0]
        gap2 = evals[2] - evals[0] if len(evals) > 2 else None
        dt = time.time() - t0
        data.append({
            'g': round(g, 5),
            'gap': round(float(gap), 10),
            'gap_x_L': round(float(gap * L), 8),
            'gap2': round(float(gap2), 10) if gap2 is not None else None,
            'E0_per_site': round(float(evals[0] / n), 8),
            'time_s': round(dt, 2)
        })
        if (gi + 1) % 10 == 0:
            print(f"      {gi+1}/{len(g_values)} done ({dt:.1f}s/pt)")
    return data


def find_crossing(data1, data2):
    """Find g where gap*L crosses."""
    from scipy.interpolate import interp1d
    g1 = [d['g'] for d in data1]
    y1 = [d['gap_x_L'] for d in data1]
    g2 = [d['g'] for d in data2]
    y2 = [d['gap_x_L'] for d in data2]

    g_min = max(min(g1), min(g2))
    g_max = min(max(g1), max(g2))
    g_grid = np.linspace(g_min, g_max, 500)

    f1 = interp1d(g1, y1, kind='cubic', fill_value='extrapolate')
    f2 = interp1d(g2, y2, kind='cubic', fill_value='extrapolate')

    diff = f1(g_grid) - f2(g_grid)
    crossings = []
    for i in range(len(diff) - 1):
        if diff[i] * diff[i + 1] < 0:
            gc = g_grid[i] - diff[i] * (g_grid[i + 1] - g_grid[i]) / (diff[i + 1] - diff[i])
            crossings.append(round(float(gc), 4))
    return crossings


# ---- MAIN ----
results = {
    'experiment': '068b_2d_square_v2',
    'description': '2D gap*L on square lattices for q=2,3,5 (GPU accelerated)',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
}

# === q=2: L=2,3,4 ===
print("=" * 60)
print("q=2: Square lattices L=2,3,4")
print("=" * 60)

q = 2
g_scan = np.linspace(0.2, 1.2, 40)
q2 = {}
for L in [2, 3, 4]:
    dim = q**(L*L)
    print(f"\n  L={L} (dim={dim})...")
    t0 = time.time()
    data = scan_gap_square(L, q, g_scan)
    dt = time.time() - t0
    q2[L] = data
    print(f"    Total: {dt:.1f}s")

print("\n  gap*L crossings (q=2):")
for L1, L2 in [(2, 3), (3, 4), (2, 4)]:
    cr = find_crossing(q2[L1], q2[L2])
    print(f"    L={L1} vs L={L2}: {cr}")
    results[f'q2_crossing_{L1}_{L2}'] = cr

results['q2_data'] = {str(L): v for L, v in q2.items()}

# === q=3: L=2,3 ===
print("\n" + "=" * 60)
print("q=3: Square lattices L=2,3")
print("=" * 60)

q = 3
g_scan_q3 = np.linspace(0.3, 1.8, 40)
q3 = {}
for L in [2, 3]:
    dim = q**(L*L)
    print(f"\n  L={L} (dim={dim})...")
    t0 = time.time()
    data = scan_gap_square(L, q, g_scan_q3)
    dt = time.time() - t0
    q3[L] = data
    print(f"    Total: {dt:.1f}s")

cr = find_crossing(q3[2], q3[3])
print(f"\n  gap*L crossing (q=3): L=2 vs L=3: {cr}")
results['q3_crossing_2_3'] = cr
results['q3_data'] = {str(L): v for L, v in q3.items()}

# === q=5: L=2,3 (KEY TEST) ===
print("\n" + "=" * 60)
print("q=5: Square lattices L=2,3  *** KEY TEST ***")
print("=" * 60)

q = 5
# L=2: dim=625 (trivial), L=3: dim=1953125 (~2M, GPU)
# Use FEWER points for L=3
g_scan_q5_small = np.linspace(0.3, 3.0, 40)
g_scan_q5_large = np.linspace(0.5, 2.5, 20)  # Fewer points for L=3

q5 = {}
print(f"\n  L=2 (dim={5**4}=625)...")
t0 = time.time()
q5[2] = scan_gap_square(2, q, g_scan_q5_small)
print(f"    Total: {time.time()-t0:.1f}s")

print(f"\n  L=3 (dim={5**9}=1953125, {'GPU' if HAS_GPU else 'CPU'})...")
print(f"    This will take a while ({len(g_scan_q5_large)} g points)...")
t0 = time.time()
q5[3] = scan_gap_square(3, q, g_scan_q5_large)
dt = time.time() - t0
print(f"    Total: {dt:.1f}s ({dt/len(g_scan_q5_large):.1f}s per g point)")

cr = find_crossing(q5[2], q5[3])
print(f"\n  gap*L crossing (q=5): L=2 vs L=3: {cr}")
results['q5_crossing_2_3'] = cr
results['q5_data'] = {str(L): v for L, v in q5.items()}

# === Analysis ===
print("\n" + "=" * 60)
print("ANALYSIS: First-order vs continuous")
print("=" * 60)

print("\nGap minimum comparison:")
for q_val, qdata in [(2, q2), (3, q3), (5, q5)]:
    print(f"\n  q={q_val}:")
    for L in sorted(qdata.keys()):
        data = qdata[L]
        gaps = [d['gap'] for d in data]
        min_gap = min(gaps)
        min_g = data[np.argmin(gaps)]['g']
        print(f"    L={L}: gap_min = {min_gap:.6e} at g={min_g:.3f}")

print("\n2D vs 1D critical points:")
gc_1d = {2: 0.250, 3: 0.333, 5: 0.441}
summary = {}
for q_val in [2, 3, 5]:
    cr = results.get(f'q{q_val}_crossing_2_3', [])
    gc_2d = cr[0] if cr else None
    ratio = gc_2d / gc_1d[q_val] if gc_2d else None
    print(f"  q={q_val}: 1D g_c={gc_1d[q_val]:.3f}, 2D g_c≈{gc_2d}, ratio={ratio:.2f}" if ratio else f"  q={q_val}: no crossing found")
    summary[f'q{q_val}'] = {'gc_1d': gc_1d[q_val], 'gc_2d': gc_2d, 'ratio': ratio}

results['summary'] = summary

# Additional: check gap*L at crossing (should be ~O(1))
print("\nGap*L at crossing point:")
for q_val, qdata in [(2, q2), (3, q3), (5, q5)]:
    cr = results.get(f'q{q_val}_crossing_2_3', [])
    if cr:
        gc = cr[0]
        for L, data in sorted(qdata.items()):
            g_arr = np.array([d['g'] for d in data])
            gL_arr = np.array([d['gap_x_L'] for d in data])
            from scipy.interpolate import interp1d
            f = interp1d(g_arr, gL_arr, kind='cubic', fill_value='extrapolate')
            val = float(f(gc))
            print(f"  q={q_val}, L={L}: gap*L(g_c={gc:.3f}) = {val:.4f}")

# Save
with open("results/sprint_068b_2d_square.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_068b_2d_square.json")
