#!/usr/bin/env python3
"""Sprint 068a: 2D hybrid Hamiltonian — validation at q=2,3

Build Potts-clock hybrid on 2D square lattice (periodic BC):
  H = -sum_<ij> delta(s_i, s_j) - g * sum_i (X_i + X_i†)

Validation:
- q=2: 2D quantum Ising. g_c should be higher than 1D (0.25)
  because z=4 neighbors doubles coupling strength vs z=2.
  Expected: g_c(2D) ≈ 0.5 (double 1D by mean-field argument).
- q=3: Should have continuous transition in 2D (maps to 3D classical).

Method: Energy gap Δ·N on periodic 2D lattices. Gap crossing = g_c.
Lattices: 2×2, 2×3, 3×3 (q=2,3). Also 3×4 for q=3.
"""
import numpy as np
import json
import time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh


def build_H_2d(Lx, Ly, q, g, periodic=True):
    """Build 2D hybrid Hamiltonian on Lx x Ly square lattice.

    Site index: flat = y * Lx + x
    Bonds: horizontal (x,y)-(x+1,y) and vertical (x,y)-(x,y+1)
    """
    n = Lx * Ly
    dim = q**n

    # Single-site operators
    # Potts 2-site diagonal: delta(s_i, s_j)
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

    # Enumerate bonds
    bonds = []
    for y in range(Ly):
        for x in range(Lx):
            site = y * Lx + x
            # Horizontal bond
            if periodic or x + 1 < Lx:
                nbr_x = (x + 1) % Lx
                nbr = y * Lx + nbr_x
                if nbr != site:
                    bonds.append((site, nbr))
            # Vertical bond
            if periodic or y + 1 < Ly:
                nbr_y = (y + 1) % Ly
                nbr = nbr_y * Lx + x
                if nbr != site:
                    bonds.append((site, nbr))

    # Coupling terms: -delta(s_i, s_j) for each bond
    for (i, j) in bonds:
        if j == i + 1:
            # Adjacent sites in flat index — use kron
            left = q**i
            right = q**(n - i - 2)
            op = potts_op
            if left > 1:
                op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1:
                op = sp_kron(op, sp_eye(right), format='csr')
        else:
            # Non-adjacent: build diagonal manually
            diag_vals = np.zeros(dim)
            for idx in range(dim):
                si = (idx // (q**i)) % q
                sj = (idx // (q**j)) % q
                diag_vals[idx] = 1.0 if si == sj else 0.0
            op = diags(diag_vals, 0, shape=(dim, dim), format='csr')
        H = H - op

    # Field terms: -g * (X_i + X_i†) for each site
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


def scan_gap(Lx, Ly, q, g_values):
    """Scan energy gap across g values."""
    n = Lx * Ly
    dim = q**n
    data = []
    for g in g_values:
        t0 = time.time()
        H = build_H_2d(Lx, Ly, q, g)
        k_eig = min(4, dim - 2)
        if dim < 500:
            # Full diag for small systems
            evals = np.linalg.eigvalsh(H.toarray())
        else:
            evals, _ = eigsh(H, k=k_eig, which='SA')
            evals = np.sort(evals)
        gap = evals[1] - evals[0]
        dt = time.time() - t0
        data.append({
            'g': round(g, 5),
            'E0': round(float(evals[0]), 8),
            'E1': round(float(evals[1]), 8),
            'gap': round(float(gap), 8),
            'gap_x_N': round(float(gap * n), 6),
            'dim': dim,
            'time_s': round(dt, 2)
        })
    return data


def find_crossing(data1, data2, Lx1, Ly1, Lx2, Ly2):
    """Find g where gap*N crosses between two lattice sizes."""
    n1, n2 = Lx1 * Ly1, Lx2 * Ly2
    g1 = [d['g'] for d in data1]
    gN1 = [d['gap_x_N'] for d in data1]
    g2 = [d['g'] for d in data2]
    gN2 = [d['gap_x_N'] for d in data2]

    # Interpolate to common g grid
    from scipy.interpolate import interp1d
    g_common = sorted(set(g1) & set(g2))
    if len(g_common) < 3:
        g_min = max(min(g1), min(g2))
        g_max = min(max(g1), max(g2))
        g_common = np.linspace(g_min, g_max, 200)

    f1 = interp1d(g1, gN1, kind='cubic', fill_value='extrapolate')
    f2 = interp1d(g2, gN2, kind='cubic', fill_value='extrapolate')

    diff = f1(g_common) - f2(g_common)
    # Find sign changes
    crossings = []
    for i in range(len(diff) - 1):
        if diff[i] * diff[i + 1] < 0:
            # Linear interpolation for crossing
            g_cross = g_common[i] - diff[i] * (g_common[i + 1] - g_common[i]) / (diff[i + 1] - diff[i])
            crossings.append(round(float(g_cross), 4))
    return crossings


# ---- Main ----
results = {
    'experiment': '068a_2d_validation',
    'description': '2D hybrid Hamiltonian validation at q=2,3',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
}

# === q=2 (2D Ising) ===
print("=" * 60)
print("q=2: 2D Quantum Ising on square lattice")
print("=" * 60)

q = 2
# Lattices: 2×2 (4 sites, dim=16), 2×3 (6, dim=64), 3×3 (9, dim=512), 2×4 (8, dim=256)
lattices_q2 = [(2, 2), (2, 3), (2, 4), (3, 3)]
g_scan = np.linspace(0.1, 1.5, 60)

q2_data = {}
for Lx, Ly in lattices_q2:
    n = Lx * Ly
    dim = q**n
    label = f"{Lx}x{Ly}"
    print(f"\n  {label} (n={n}, dim={dim})...")
    t0 = time.time()
    data = scan_gap(Lx, Ly, q, g_scan)
    dt = time.time() - t0
    q2_data[label] = data

    # Find gap minimum
    gaps = [d['gap'] for d in data]
    min_idx = np.argmin(gaps)
    print(f"    Gap minimum: {gaps[min_idx]:.6f} at g={data[min_idx]['g']:.3f}")
    print(f"    Time: {dt:.1f}s")

# Find crossings
print("\n  Gap*N crossings (q=2):")
crossing_pairs = [
    ('2x2', '2x3', 2, 2, 2, 3),
    ('2x3', '2x4', 2, 3, 2, 4),
    ('2x3', '3x3', 2, 3, 3, 3),
    ('2x4', '3x3', 2, 4, 3, 3),
]
q2_crossings = {}
for label1, label2, Lx1, Ly1, Lx2, Ly2 in crossing_pairs:
    crossings = find_crossing(q2_data[label1], q2_data[label2], Lx1, Ly1, Lx2, Ly2)
    q2_crossings[f"{label1}_{label2}"] = crossings
    print(f"    {label1} vs {label2}: {crossings}")

results['q2'] = {
    'lattices': {k: v for k, v in q2_data.items()},
    'crossings': q2_crossings,
    'note': '1D g_c(q=2) = 0.25. 2D should be higher due to z=4 neighbors.'
}

# === q=3 ===
print("\n" + "=" * 60)
print("q=3: 2D Potts on square lattice")
print("=" * 60)

q = 3
# 2×2 (81), 2×3 (729), 2×4 (6561), 3×3 (19683)
lattices_q3 = [(2, 2), (2, 3), (2, 4), (3, 3)]
g_scan_q3 = np.linspace(0.1, 2.0, 60)

q3_data = {}
for Lx, Ly in lattices_q3:
    n = Lx * Ly
    dim = q**n
    label = f"{Lx}x{Ly}"
    print(f"\n  {label} (n={n}, dim={dim})...")
    t0 = time.time()
    data = scan_gap(Lx, Ly, q, g_scan_q3)
    dt = time.time() - t0
    q3_data[label] = data

    gaps = [d['gap'] for d in data]
    min_idx = np.argmin(gaps)
    print(f"    Gap minimum: {gaps[min_idx]:.6f} at g={data[min_idx]['g']:.3f}")
    print(f"    Time: {dt:.1f}s")

# Find crossings
print("\n  Gap*N crossings (q=3):")
q3_crossings = {}
for label1, label2, Lx1, Ly1, Lx2, Ly2 in crossing_pairs:
    crossings = find_crossing(q3_data[label1], q3_data[label2], Lx1, Ly1, Lx2, Ly2)
    q3_crossings[f"{label1}_{label2}"] = crossings
    print(f"    {label1} vs {label2}: {crossings}")

results['q3'] = {
    'lattices': {k: v for k, v in q3_data.items()},
    'crossings': q3_crossings,
    'note': '1D g_c(q=3) = 0.333. 2D should be higher.'
}

# === Summary ===
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)

# q=2 expected: g_c(1D) = 0.25, g_c(2D) should be ~0.5 (z doubles)
# q=3 expected: g_c(1D) = 0.333, g_c(2D) should be ~0.66
print("\nq=2 1D g_c=0.250. 2D crossings suggest g_c ≈", q2_crossings)
print("q=3 1D g_c=0.333. 2D crossings suggest g_c ≈", q3_crossings)

# Save results
with open("results/sprint_068a_2d_validation.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_068a_2d_validation.json")
