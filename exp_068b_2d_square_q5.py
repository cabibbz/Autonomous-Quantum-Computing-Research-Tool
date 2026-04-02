#!/usr/bin/env python3
"""Sprint 068b: 2D phase transition on square lattices — q=2,3,5

Key insight from 068a: strip geometries (2×N) are quasi-1D.
Need SQUARE lattices (L×L) with gap*L scaling for 2D physics.

In 2D, gap ~ L^{-z} with z=1 at Lorentz-invariant critical point.
So gap*L should be scale-invariant at g_c.

Lattices:
  q=2: 2×2 (16), 3×3 (512), 4×4 (65536) — all trivial
  q=3: 2×2 (81), 3×3 (19683) — CPU feasible
  q=5: 2×2 (625), 3×3 (~2M) — CPU/GPU feasible

KEY QUESTION: Does q=5 show continuous transition (gap*L crossing)
or first-order (gap minimum shrinking exponentially)?
"""
import numpy as np
import json
import time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh


def build_H_2d(Lx, Ly, q, g):
    """Build 2D hybrid Hamiltonian on Lx x Ly periodic square lattice."""
    n = Lx * Ly
    dim = q**n

    # Potts coupling: delta(s_i, s_j)
    potts_2site = np.zeros(q**2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q**2, q**2), format='csr')

    # Clock field: X + X†
    X = np.zeros((q, q))
    for s in range(q):
        X[(s + 1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)

    H = csr_matrix((dim, dim))

    # Enumerate bonds on periodic square lattice
    bonds = []
    for y in range(Ly):
        for x in range(Lx):
            site = y * Lx + x
            # Horizontal
            nbr = y * Lx + (x + 1) % Lx
            if nbr != site:
                bonds.append((min(site, nbr), max(site, nbr)))
            # Vertical
            nbr = ((y + 1) % Ly) * Lx + x
            if nbr != site:
                bonds.append((min(site, nbr), max(site, nbr)))
    bonds = list(set(bonds))  # Remove duplicates

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


def scan_gap_square(L, q, g_values):
    """Scan energy gap for L×L lattice. Return gap*L."""
    n = L * L
    dim = q**n
    data = []
    for g in g_values:
        t0 = time.time()
        H = build_H_2d(L, L, q, g)
        k_eig = min(6, dim - 2)
        if dim < 500:
            evals = np.linalg.eigvalsh(H.toarray())[:6]
        else:
            evals, _ = eigsh(H, k=k_eig, which='SA')
            evals = np.sort(evals)
        gap = evals[1] - evals[0]
        gap2 = evals[2] - evals[0] if len(evals) > 2 else None
        dt = time.time() - t0
        data.append({
            'g': round(g, 5),
            'gap': round(float(gap), 10),
            'gap_x_L': round(float(gap * L), 8),
            'gap2': round(float(gap2), 10) if gap2 is not None else None,
            'gap2_x_L': round(float(gap2 * L), 8) if gap2 is not None else None,
            'E0_per_site': round(float(evals[0] / n), 8),
            'time_s': round(dt, 2)
        })
    return data


def find_crossing_gapL(data_small, data_large, L_small, L_large):
    """Find g where gap*L crosses between two square lattice sizes."""
    from scipy.interpolate import interp1d
    g1 = [d['g'] for d in data_small]
    y1 = [d['gap_x_L'] for d in data_small]
    g2 = [d['g'] for d in data_large]
    y2 = [d['gap_x_L'] for d in data_large]

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
    'experiment': '068b_2d_square_q5',
    'description': '2D gap*L on square lattices for q=2,3,5',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
}

# === q=2: 2×2, 3×3, 4×4 ===
print("=" * 60)
print("q=2: Square lattices 2×2, 3×3, 4×4")
print("=" * 60)

q = 2
g_scan = np.linspace(0.2, 1.2, 50)
q2 = {}
for L in [2, 3, 4]:
    dim = q**(L*L)
    print(f"\n  L={L} ({L}×{L}, dim={dim})...")
    t0 = time.time()
    data = scan_gap_square(L, q, g_scan)
    dt = time.time() - t0
    q2[L] = data
    # Report gap*L at a few g values
    mid = len(data) // 2
    print(f"    gap*L at g={data[mid]['g']:.2f}: {data[mid]['gap_x_L']:.4f}")
    print(f"    Total: {dt:.1f}s")

print("\n  gap*L crossings (q=2 square):")
for L1, L2 in [(2, 3), (3, 4), (2, 4)]:
    cr = find_crossing_gapL(q2[L1], q2[L2], L1, L2)
    print(f"    L={L1} vs L={L2}: {cr}")
    results[f'q2_crossing_{L1}_{L2}'] = cr

results['q2_data'] = {str(L): v for L, v in q2.items()}

# === q=3: 2×2, 3×3 ===
print("\n" + "=" * 60)
print("q=3: Square lattices 2×2, 3×3")
print("=" * 60)

q = 3
g_scan = np.linspace(0.2, 1.8, 50)
q3 = {}
for L in [2, 3]:
    dim = q**(L*L)
    print(f"\n  L={L} ({L}×{L}, dim={dim})...")
    t0 = time.time()
    data = scan_gap_square(L, q, g_scan)
    dt = time.time() - t0
    q3[L] = data
    mid = len(data) // 2
    print(f"    gap*L at g={data[mid]['g']:.2f}: {data[mid]['gap_x_L']:.4f}")
    print(f"    Total: {dt:.1f}s")

cr = find_crossing_gapL(q3[2], q3[3], 2, 3)
print(f"\n  gap*L crossing (q=3, L=2 vs L=3): {cr}")
results['q3_crossing_2_3'] = cr
results['q3_data'] = {str(L): v for L, v in q3.items()}

# === q=5: 2×2, 3×3 (KEY TEST) ===
print("\n" + "=" * 60)
print("q=5: Square lattices 2×2, 3×3  ** KEY TEST **")
print("=" * 60)

q = 5
g_scan_q5 = np.linspace(0.3, 3.0, 50)
q5 = {}
for L in [2, 3]:
    dim = q**(L*L)
    print(f"\n  L={L} ({L}×{L}, dim={dim})...")
    t0 = time.time()
    if dim > 500000:
        # 3×3 at q=5: dim = 5^9 = 1953125, use eigsh
        print(f"    Large system (dim={dim}), using sparse eigsh...")
    data = scan_gap_square(L, q, g_scan_q5)
    dt = time.time() - t0
    mid = len(data) // 2
    print(f"    gap*L at g={data[mid]['g']:.2f}: {data[mid]['gap_x_L']:.4f}")
    print(f"    Total: {dt:.1f}s")
    q5[L] = data

cr = find_crossing_gapL(q5[2], q5[3], 2, 3)
print(f"\n  gap*L crossing (q=5, L=2 vs L=3): {cr}")
results['q5_crossing_2_3'] = cr
results['q5_data'] = {str(L): v for L, v in q5.items()}

# === Analysis: continuous vs first-order ===
print("\n" + "=" * 60)
print("FIRST-ORDER TEST: Gap minimum scaling")
print("=" * 60)

for q_val, qdata in [(2, q2), (3, q3), (5, q5)]:
    print(f"\n  q={q_val}:")
    for L, data in sorted(qdata.items()):
        gaps = [d['gap'] for d in data]
        min_gap = min(gaps)
        min_g = data[np.argmin(gaps)]['g']
        print(f"    L={L}: gap_min = {min_gap:.6e} at g={min_g:.3f}")

print("\n  First-order signature: gap_min ~ exp(-α·L²)")
print("  Continuous signature: gap_min ~ L^{-z}")

# === Comparison: 2D vs 1D g_c ===
print("\n" + "=" * 60)
print("2D vs 1D critical points")
print("=" * 60)
gc_1d = {2: 0.250, 3: 0.333, 5: 0.441}
print(f"  q=2: 1D g_c=0.250, 2D crossings: L=2vs3={results.get('q2_crossing_2_3')}, L=3vs4={results.get('q2_crossing_3_4')}")
print(f"  q=3: 1D g_c=0.333, 2D crossing: L=2vs3={results.get('q3_crossing_2_3')}")
print(f"  q=5: 1D g_c=0.441, 2D crossing: L=2vs3={results.get('q5_crossing_2_3')}")
print(f"\n  Ratio 2D/1D should be ~2 (doubling coordination number z=2→4)")

for q_val in [2, 3, 5]:
    cr_key = f'q{q_val}_crossing_2_3'
    cr_val = results.get(cr_key, [])
    if cr_val:
        ratio = cr_val[0] / gc_1d[q_val]
        print(f"    q={q_val}: ratio = {ratio:.2f}")

# Save
with open("results/sprint_068b_2d_square.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_068b_2d_square.json")
