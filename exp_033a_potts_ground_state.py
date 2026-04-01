"""
Sprint 033a: 3-State Potts Model Ground State Entanglement Sweep

H = -J Σ δ(σ_i,σ_{i+1}) - h Σ (P_i + P_i†)
S₃ symmetry, critical point at self-dual h/J, second-order (c=4/5 CFT for q=3).

System: n=6 (dim=729), n=8 (dim=6561) for scaling check.
"""

import numpy as np
from scipy import linalg as la
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import eigsh
import json, time

t0 = time.time()

# ---- Sparse Potts Hamiltonian ----
def sparse_potts_hamiltonian(n, h_over_J, boundary='open'):
    """Build Potts Hamiltonian using sparse matrices for efficiency."""
    d = 3
    dim = d**n

    # Single-site operators (sparse)
    P = csr_matrix(np.array([[0,0,1],[1,0,0],[0,1,0]], dtype=float))  # |0>->|1>->|2>->|0>
    Pd = P.T.tocsr()
    I_d = sp_eye(d, format='csr')

    # Potts delta: Σ_a |aa><aa|
    delta = csr_matrix((d*d, d*d))
    for a in range(d):
        idx = a * d + a
        row = csr_matrix(([1.0], ([idx], [idx])), shape=(d*d, d*d))
        delta = delta + row

    H = csr_matrix((dim, dim))

    # Interaction: -J δ(σ_i, σ_{i+1})
    n_bonds = n - 1 if boundary == 'open' else n
    for bond in range(n_bonds):
        j = (bond + 1) % n
        if j == bond + 1:
            prefix = sp_eye(d**bond, format='csr') if bond > 0 else sp_eye(1, format='csr')
            suffix = sp_eye(d**(n - bond - 2), format='csr') if bond + 2 < n else sp_eye(1, format='csr')
            H = H - sp_kron(sp_kron(prefix, delta), suffix, format='csr')
        else:
            # PBC wrap: sites 0 and n-1
            for a in range(d):
                proj = csr_matrix(([1.0], ([a], [a])), shape=(d, d))
                middle = sp_eye(d**(n-2), format='csr')
                H = H - sp_kron(sp_kron(proj, middle), proj, format='csr')

    # Transverse field: -h (P + P†)
    PPd = P + Pd
    for site in range(n):
        prefix = sp_eye(d**site, format='csr') if site > 0 else sp_eye(1, format='csr')
        suffix = sp_eye(d**(n - site - 1), format='csr') if site < n - 1 else sp_eye(1, format='csr')
        H = H - h_over_J * sp_kron(sp_kron(prefix, PPd), suffix, format='csr')

    return H

def partial_trace_general(rho, dims, keep):
    n = len(dims)
    trace_out = [i for i in range(n) if i not in keep]
    shape = list(dims) + list(dims)
    rho_tensor = rho.reshape(shape)
    for q in sorted(trace_out, reverse=True):
        n_remaining = rho_tensor.ndim // 2
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q + n_remaining)
    dim_keep = int(np.prod([dims[k] for k in keep]))
    return rho_tensor.reshape(dim_keep, dim_keep)

def von_neumann_entropy(rho):
    eigvals = la.eigvalsh(rho)
    eigvals = eigvals[eigvals > 1e-15]
    return -np.sum(eigvals * np.log2(eigvals))

def mutual_information(rho_full, dims, i, j):
    rho_i = partial_trace_general(rho_full, dims, [i])
    rho_j = partial_trace_general(rho_full, dims, [j])
    rho_ij = partial_trace_general(rho_full, dims, [i, j])
    return von_neumann_entropy(rho_i) + von_neumann_entropy(rho_j) - von_neumann_entropy(rho_ij)

def tripartite_info(rho_full, dims, i, j, k):
    rho_i = partial_trace_general(rho_full, dims, [i])
    rho_j = partial_trace_general(rho_full, dims, [j])
    rho_k = partial_trace_general(rho_full, dims, [k])
    rho_ij = partial_trace_general(rho_full, dims, [i, j])
    rho_ik = partial_trace_general(rho_full, dims, [i, k])
    rho_jk = partial_trace_general(rho_full, dims, [j, k])
    rho_ijk = partial_trace_general(rho_full, dims, [i, j, k])
    S = lambda r: von_neumann_entropy(r)
    return S(rho_i) + S(rho_j) + S(rho_k) - S(rho_ij) - S(rho_ik) - S(rho_jk) + S(rho_ijk)

def analyze_point(n, d, h_J, boundary='open'):
    """Full entanglement analysis at one coupling value."""
    dims = [d] * n
    n_A = n // 2
    H = sparse_potts_hamiltonian(n, h_J, boundary)

    # Ground state via sparse diag (get 5 lowest for gap analysis)
    n_eig = min(6, d**n - 2)
    eigvals_H, eigvecs_H = eigsh(H, k=n_eig, which='SA')
    sort_idx = np.argsort(eigvals_H)
    eigvals_H = eigvals_H[sort_idx]
    eigvecs_H = eigvecs_H[:, sort_idx]
    psi_gs = eigvecs_H[:, 0]
    e0 = eigvals_H[0]

    # Gap to first state above ground manifold
    tol = 1e-6
    n_degen = 1
    for k in range(1, len(eigvals_H)):
        if eigvals_H[k] - e0 > tol:
            gap = eigvals_H[k] - e0
            n_degen = k
            break
    else:
        gap = 0.0

    rho = np.outer(psi_gs, psi_gs.conj())
    rho_A = partial_trace_general(rho, dims, list(range(n_A)))
    S_half = von_neumann_entropy(rho_A)

    # Schmidt spectrum
    rho_eigvals = la.eigvalsh(rho_A)
    rho_eigvals = rho_eigvals[rho_eigvals > 1e-15]
    schmidt_rank = len(rho_eigvals)

    # Pairwise MI
    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            mi_vals.append(mutual_information(rho, dims, i, j))
    mi_mean = np.mean(mi_vals)
    mi_cv = np.std(mi_vals) / mi_mean if mi_mean > 1e-10 else 0

    # I3 for consecutive triple
    i3_012 = tripartite_info(rho, dims, 0, 1, 2)

    return {
        'h_J': float(h_J), 'e0': float(e0), 'gap': float(gap),
        'n_degen': int(n_degen), 'S': float(S_half),
        'schmidt_rank': int(schmidt_rank),
        'MI_mean': float(mi_mean), 'MI_cv': float(mi_cv),
        'I3_012': float(i3_012),
        'spectrum': sorted([float(x) for x in rho_eigvals], reverse=True),
        'energies': [float(x) for x in eigvals_H[:6]]
    }

# ---- n=6 sweep ----
print("=== Sprint 033a: 3-State Potts Ground State (n=6) ===")
n6 = 6; d = 3

h_values_6 = np.concatenate([
    np.linspace(0.01, 0.10, 5),
    np.linspace(0.12, 0.40, 15),
    np.linspace(0.45, 1.0, 8),
    [1.5, 2.0, 3.0]
])

results_6 = []
print("\n  h/J     S       gap     deg  rank  MI_cv   I3(012)  time")
for h_J in h_values_6:
    t1 = time.time()
    r = analyze_point(n6, d, h_J)
    dt = time.time() - t1
    results_6.append(r)
    print(f"  {r['h_J']:5.3f}  {r['S']:6.3f}  {r['gap']:7.4f}  {r['n_degen']:2d}   "
          f"{r['schmidt_rank']:2d}   {r['MI_cv']:6.3f}  {r['I3_012']:+7.4f}  {dt:.1f}s")

# ---- Find transition features ----
print("\n=== Transition Analysis (n=6) ===")
S_vals = [r['S'] for r in results_6]
h_vals = [r['h_J'] for r in results_6]
gaps = [r['gap'] for r in results_6]

# I3 sign change
for i in range(len(results_6)-1):
    if results_6[i]['I3_012'] > 0 and results_6[i+1]['I3_012'] <= 0:
        print(f"I3 sign change: h/J ∈ [{results_6[i]['h_J']:.3f}, {results_6[i+1]['h_J']:.3f}]")
    if results_6[i]['I3_012'] >= 0 and results_6[i+1]['I3_012'] < 0:
        print(f"I3 goes negative: h/J ∈ [{results_6[i]['h_J']:.3f}, {results_6[i+1]['h_J']:.3f}]")

# Gap minimum (for non-degenerate region)
nondegen = [(r['h_J'], r['gap']) for r in results_6 if r['n_degen'] == 1]
if nondegen:
    min_gap = min(nondegen, key=lambda x: x[1])
    print(f"Min gap (non-degenerate): Δ={min_gap[1]:.4f} at h/J={min_gap[0]:.3f}")

# Steepest entropy drop
dS = [(S_vals[i+1]-S_vals[i])/(h_vals[i+1]-h_vals[i]) for i in range(len(S_vals)-1)]
idx_steep = np.argmin(dS)
print(f"Steepest dS/d(h/J) = {dS[idx_steep]:.2f} at h/J ≈ {(h_vals[idx_steep]+h_vals[idx_steep+1])/2:.3f}")

# ---- n=8 check at key points ----
print("\n=== n=8 Finite-Size Check ===")
n8 = 8
key_h8 = [0.01, 0.15, 0.25, 0.50, 1.0]

results_8 = []
for h_J in key_h8:
    t1 = time.time()
    r = analyze_point(n8, d, h_J)
    dt = time.time() - t1
    results_8.append(r)
    print(f"  h/J={r['h_J']:.2f}: S={r['S']:.4f}, gap={r['gap']:.4f}, "
          f"degen={r['n_degen']}, MI_cv={r['MI_cv']:.3f} [{dt:.1f}s]")

# Compare n=6 vs n=8
print("\n=== Size Comparison ===")
for h_J in key_h8:
    r6 = next((r for r in results_6 if abs(r['h_J'] - h_J) < 0.02), None)
    r8 = next((r for r in results_8 if abs(r['h_J'] - h_J) < 0.02), None)
    if r6 and r8:
        print(f"  h/J={h_J:.2f}: S(n=6)={r6['S']:.4f}, S(n=8)={r8['S']:.4f}, "
              f"ratio={r8['S']/r6['S']:.3f}" if r6['S'] > 0.01 else
              f"  h/J={h_J:.2f}: S(n=6)={r6['S']:.4f}, S(n=8)={r8['S']:.4f}")

# ---- Save ----
all_results = {
    'experiment': '033a',
    'description': '3-state Potts ground state entanglement sweep',
    'parameters': {'n_values': [6, 8], 'local_dim': 3, 'boundary': 'open'},
    'n6_sweep': results_6,
    'n8_key_points': results_8,
    'runtime_seconds': time.time() - t0
}

with open('results/sprint_033a_potts_ground_state.json', 'w') as f:
    json.dump(all_results, f, indent=2)

print(f"\nSaved. Total runtime: {time.time()-t0:.1f}s")
