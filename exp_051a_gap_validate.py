#!/usr/bin/env python3
"""Sprint 051a: Validate energy gap method on q=2,3 Potts (known g_c).

Energy gap Δ = E₁ - E₀ closes as 1/N at criticality.
So Δ·N is scale-invariant at g_c → curves for different N cross there.

q=2: g_c = 0.25 (our Potts convention, = J/4)
q=3: g_c = 1/3  (self-duality, = J/3)

Use exact diag for small n (faster, no χ issues).
"""
import numpy as np, json, time, warnings
warnings.filterwarnings('ignore')
from scipy.sparse.linalg import eigsh
from scipy.sparse import kron, eye, csr_matrix

def potts_hamiltonian(q, n, g, J=1.0):
    """Build H = -J Σ δ(s_i,s_j) - g Σ (X_i + X_i†) as sparse matrix."""
    dim = q**n
    # Single-site operators
    X = np.zeros((q, q), dtype=complex)
    for a in range(q):
        X[(a+1) % q, a] = 1.0
    Xphc = X + X.conj().T  # X + X†

    # Projectors for delta coupling
    projectors = []
    for a in range(q):
        P = np.zeros((q, q)); P[a, a] = 1.0
        projectors.append(csr_matrix(P))

    H = csr_matrix((dim, dim), dtype=complex)
    I = csr_matrix(np.eye(q))

    # Nearest-neighbor coupling: -J Σ_a P_a ⊗ P_a
    for i in range(n - 1):
        for a in range(q):
            op = csr_matrix(np.eye(1))
            for j in range(n):
                if j == i:
                    op = kron(op, projectors[a])
                elif j == i + 1:
                    op = kron(op, projectors[a])
                else:
                    op = kron(op, I)
            H += -J * op

    # Transverse field: -g Σ (X + X†)
    Xphc_sp = csr_matrix(Xphc)
    for i in range(n):
        op = csr_matrix(np.eye(1))
        for j in range(n):
            if j == i:
                op = kron(op, Xphc_sp)
            else:
                op = kron(op, I)
        H += -g * op

    return H

def energy_gap(q, n, g, J=1.0):
    """Compute E₁ - E₀ using sparse diag."""
    H = potts_hamiltonian(q, n, g, J)
    # Need 2 lowest eigenvalues
    vals = eigsh(H.real if q == 2 else H, k=min(4, H.shape[0]-1), which='SA', return_eigenvectors=False)
    vals = np.sort(vals.real)
    return vals[1] - vals[0]

# Timing test
print("=== Timing test ===", flush=True)
t0 = time.time()
gap = energy_gap(2, 8, 0.25)
dt = time.time() - t0
print(f"q=2 n=8: gap={gap:.6f}, t={dt:.1f}s", flush=True)

t0 = time.time()
gap = energy_gap(3, 6, 0.33)
dt = time.time() - t0
print(f"q=3 n=6: gap={gap:.6f}, t={dt:.1f}s", flush=True)

# Check if q=3 n=8 (3^8=6561) is feasible
t0 = time.time()
gap = energy_gap(3, 8, 0.33)
dt = time.time() - t0
print(f"q=3 n=8: gap={gap:.6f}, t={dt:.1f}s (dim=6561)", flush=True)

results = {}

# === q=2 Potts: g_c = 0.25 ===
print("\n=== q=2 Potts (g_c = 0.25) ===", flush=True)
g_vals_q2 = np.linspace(0.10, 0.45, 36)
for n in [4, 6, 8, 10]:
    print(f"  n={n}...", end="", flush=True)
    t0 = time.time()
    gaps = []
    for g in g_vals_q2:
        gap = energy_gap(2, n, g)
        gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * n)})
    dt = time.time() - t0
    print(f" done ({dt:.1f}s)", flush=True)
    results[f'q2_n{n}'] = gaps

# Find crossings for q=2
print("\n  q=2 crossings (Δ·N):", flush=True)
for i, n1 in enumerate([4, 6, 8]):
    for n2 in [6, 8, 10][i:]:
        if n1 == n2: continue
        g1_data = results[f'q2_n{n1}']
        g2_data = results[f'q2_n{n2}']
        for j in range(len(g_vals_q2) - 1):
            diff_j = g1_data[j]['gap_x_n'] - g2_data[j]['gap_x_n']
            diff_j1 = g1_data[j+1]['gap_x_n'] - g2_data[j+1]['gap_x_n']
            if diff_j * diff_j1 < 0:
                # Linear interpolation
                g_cross = g_vals_q2[j] - diff_j * (g_vals_q2[j+1] - g_vals_q2[j]) / (diff_j1 - diff_j)
                print(f"    n={n1},{n2}: g_cross = {g_cross:.4f}", flush=True)

# Save incrementally
with open('results/sprint_051a_gap_validate.json', 'w') as f:
    json.dump({'sprint': '051a', 'method': 'exact_diag_gap', 'data': results}, f, indent=2)

# === q=3 Potts: g_c = 1/3 ===
print("\n=== q=3 Potts (g_c = 1/3 ≈ 0.333) ===", flush=True)
g_vals_q3 = np.linspace(0.15, 0.55, 41)
for n in [4, 6, 8]:
    dim = 3**n
    print(f"  n={n} (dim={dim})...", end="", flush=True)
    t0 = time.time()
    gaps = []
    for g in g_vals_q3:
        gap = energy_gap(3, n, g)
        gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * n)})
    dt = time.time() - t0
    print(f" done ({dt:.1f}s)", flush=True)
    results[f'q3_n{n}'] = gaps

# Find crossings for q=3
print("\n  q=3 crossings (Δ·N):", flush=True)
for n1, n2 in [(4,6), (4,8), (6,8)]:
    g1_data = results[f'q3_n{n1}']
    g2_data = results[f'q3_n{n2}']
    for j in range(len(g_vals_q3) - 1):
        diff_j = g1_data[j]['gap_x_n'] - g2_data[j]['gap_x_n']
        diff_j1 = g1_data[j+1]['gap_x_n'] - g2_data[j+1]['gap_x_n']
        if diff_j * diff_j1 < 0:
            g_cross = g_vals_q3[j] - diff_j * (g_vals_q3[j+1] - g_vals_q3[j]) / (diff_j1 - diff_j)
            print(f"    n={n1},{n2}: g_cross = {g_cross:.4f}", flush=True)

# Save final
with open('results/sprint_051a_gap_validate.json', 'w') as f:
    json.dump({'sprint': '051a', 'method': 'exact_diag_gap',
               'q2_gc_exact': 0.25, 'q3_gc_exact': 1/3,
               'data': results}, f, indent=2)

print("\nDone! Results saved.", flush=True)
