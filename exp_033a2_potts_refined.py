"""
Sprint 033a (refined): Find actual Potts critical point and entanglement signatures.

The first run showed entropy monotonically decreasing from log2(3) in ordered phase.
Need: (1) denser scan in transition region, (2) proper gap (to excited states above
the 3-fold degenerate ground manifold), (3) n=8 comparison for finite-size effects.
"""

import numpy as np
from scipy import linalg as la
import json, time

t0 = time.time()

def potts_operators():
    d = 3
    I3 = np.eye(d)
    P = np.zeros((d, d))
    P[1, 0] = 1; P[2, 1] = 1; P[0, 2] = 1
    Pd = P.T
    return I3, P, Pd

def potts_hamiltonian(n, h_over_J, boundary='open'):
    d = 3
    I3, P, Pd = potts_operators()
    dim = d**n
    H = np.zeros((dim, dim))

    # Potts interaction: -J δ(σ_i, σ_{i+1})
    n_bonds = n - 1 if boundary == 'open' else n
    for bond in range(n_bonds):
        i = bond
        j = (bond + 1) % n
        if j == i + 1:
            # Adjacent sites
            prefix = np.eye(d**i) if i > 0 else np.eye(1)
            suffix = np.eye(d**(n - i - 2)) if i + 2 < n else np.eye(1)
            delta_ij = np.zeros((d*d, d*d))
            for a in range(d):
                idx = a * d + a
                delta_ij[idx, idx] = 1.0
            H -= np.kron(np.kron(prefix, delta_ij), suffix)
        else:
            # PBC wrap
            for a in range(d):
                proj_0 = np.zeros((d, d)); proj_0[a, a] = 1
                proj_n = np.zeros((d, d)); proj_n[a, a] = 1
                middle = np.eye(d**(n-2))
                H -= np.kron(np.kron(proj_0, middle), proj_n)

    # Transverse field: -h (P_i + P_i†)
    for site in range(n):
        prefix = np.eye(d**site) if site > 0 else np.eye(1)
        suffix = np.eye(d**(n - site - 1)) if site < n - 1 else np.eye(1)
        H -= h_over_J * np.kron(np.kron(prefix, P + Pd), suffix)

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

# ---- Sweep for n=6 with denser grid in transition region ----
n = 6
d = 3
dims = [d] * n
n_A = n // 2

# Dense sweep focused on transition region (0.01 to 0.6) + coarse above
h_values = np.concatenate([
    np.linspace(0.01, 0.08, 8),    # deep ordered
    np.linspace(0.10, 0.40, 16),   # transition region
    np.linspace(0.45, 0.80, 8),    # above transition
    np.linspace(1.0, 3.0, 5)       # deep disordered
])

print(f"=== Sprint 033a (refined): 3-State Potts Transition (n={n}) ===")
print(f"Sweeping {len(h_values)} points from {h_values[0]:.3f} to {h_values[-1]:.1f}")

results = []
for h_J in h_values:
    H = potts_hamiltonian(n, h_J, boundary='open')
    eigvals_H, eigvecs_H = la.eigh(H)
    psi_gs = eigvecs_H[:, 0]
    e0 = eigvals_H[0]

    # Proper gap: to first state ABOVE the ground state manifold
    # (skip near-degenerate states within tolerance)
    tol = 1e-8
    for k in range(1, len(eigvals_H)):
        if eigvals_H[k] - e0 > tol:
            gap = eigvals_H[k] - e0
            n_degen = k  # number of quasi-degenerate ground states
            break

    rho = np.outer(psi_gs, psi_gs.conj())
    rho_A = partial_trace_general(rho, dims, list(range(n_A)))
    S_half = von_neumann_entropy(rho_A)

    # Schmidt rank (number of significant eigenvalues)
    rho_eigvals = la.eigvalsh(rho_A)
    schmidt_rank = np.sum(rho_eigvals > 1e-10)

    # MI for all pairs
    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            mi_vals.append(mutual_information(rho, dims, i, j))
    mi_mean = np.mean(mi_vals)
    mi_cv = np.std(mi_vals) / mi_mean if mi_mean > 1e-10 else 0

    # I3 for consecutive triple
    i3_012 = tripartite_info(rho, dims, 0, 1, 2)

    results.append({
        'h_J': float(h_J), 'S': float(S_half), 'gap': float(gap),
        'n_degen': int(n_degen), 'schmidt_rank': int(schmidt_rank),
        'MI_mean': float(mi_mean), 'MI_cv': float(mi_cv), 'I3_012': float(i3_012)
    })

# Print key results
print("\n  h/J     S       gap     degen  rank  MI_cv   I3(012)")
print("  " + "-"*65)
for r in results:
    print(f"  {r['h_J']:5.3f}  {r['S']:6.3f}  {r['gap']:7.4f}  {r['n_degen']:3d}    "
          f"{r['schmidt_rank']:2d}   {r['MI_cv']:6.3f}  {r['I3_012']:+6.3f}")

# ---- Find transition features ----
print("\n=== Transition Analysis ===")

# Where does I3 change sign?
for i in range(len(results)-1):
    if results[i]['I3_012'] > 0 and results[i+1]['I3_012'] <= 0:
        print(f"I3 sign change: h/J ∈ [{results[i]['h_J']:.3f}, {results[i+1]['h_J']:.3f}]")

# Where is the gap minimum (excluding degeneracy region)?
gaps = [r['gap'] for r in results if r['n_degen'] == 1]
h_nondegen = [r['h_J'] for r in results if r['n_degen'] == 1]
if gaps:
    idx_min = np.argmin(gaps)
    print(f"Gap minimum (non-degenerate): Δ={gaps[idx_min]:.4f} at h/J={h_nondegen[idx_min]:.3f}")

# Where does entropy have steepest drop?
S_vals = [r['S'] for r in results]
h_vals = [r['h_J'] for r in results]
dS_dh = [(S_vals[i+1]-S_vals[i])/(h_vals[i+1]-h_vals[i]) for i in range(len(results)-1)]
idx_steep = np.argmin(dS_dh)
print(f"Steepest entropy drop: dS/d(h/J)={dS_dh[idx_steep]:.2f} at h/J≈{(h_vals[idx_steep]+h_vals[idx_steep+1])/2:.3f}")

# Where does Schmidt rank change?
for i in range(len(results)-1):
    if results[i]['schmidt_rank'] != results[i+1]['schmidt_rank']:
        print(f"Schmidt rank change: {results[i]['schmidt_rank']}→{results[i+1]['schmidt_rank']} "
              f"at h/J∈[{results[i]['h_J']:.3f},{results[i+1]['h_J']:.3f}]")

# ---- Now do n=8 for key points to check finite-size scaling ----
print("\n=== n=8 comparison (key points only) ===")
n8 = 8
dims8 = [d] * n8
n_A8 = n8 // 2
key_h = [0.01, 0.10, 0.20, 0.30, 0.50, 1.0]

t_n8 = time.time()
# Check timing for n=8 (3^8 = 6561)
H8_test = potts_hamiltonian(n8, 0.5, 'open')
print(f"n=8 Hamiltonian shape: {H8_test.shape}, build time: {time.time()-t_n8:.1f}s")

n8_results = []
for h_J in key_h:
    t1 = time.time()
    H8 = potts_hamiltonian(n8, h_J, 'open')
    eigvals8, eigvecs8 = la.eigh(H8)
    psi8 = eigvecs8[:, 0]
    rho8 = np.outer(psi8, psi8.conj())
    rho_A8 = partial_trace_general(rho8, dims8, list(range(n_A8)))
    S8 = von_neumann_entropy(rho_A8)

    for k in range(1, len(eigvals8)):
        if eigvals8[k] - eigvals8[0] > 1e-8:
            gap8 = eigvals8[k] - eigvals8[0]
            n_degen8 = k
            break

    dt = time.time() - t1
    n8_results.append({'h_J': h_J, 'S': S8, 'gap': gap8, 'n_degen': n_degen8})
    print(f"  h/J={h_J:.2f}: S={S8:.4f}, gap={gap8:.4f}, degen={n_degen8} [{dt:.1f}s]")

# ---- Save all results ----
all_results = {
    'experiment': '033a_refined',
    'description': '3-state Potts ground state entanglement sweep (refined)',
    'n6_sweep': results,
    'n8_key_points': n8_results,
    'runtime_seconds': time.time() - t0
}

with open('results/sprint_033a_potts_ground_state.json', 'w') as f:
    json.dump(all_results, f, indent=2)

print(f"\nResults saved. Total runtime: {time.time()-t0:.1f}s")
