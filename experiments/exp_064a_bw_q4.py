"""
Sprint 064a: BW Locality for q=4 Potts at g_c=0.392

Extends Sprint 032-034 BW analysis to q=4 Potts (novel CFT regime boundary).
System: n=8, subsystem A = left 4 sites, dim_A = 4^4 = 256.
Full dim = 4^8 = 65536. Sparse eigsh for ground state.

Key questions:
- Does BW locality decrease from q=3's 76.5%?
- What is the coupling profile (Unruh-like gradient)?
- Which envelope (sin_inv, linear, uniform) fits best?
"""

import numpy as np
from scipy import linalg as la
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import eigsh
import json, time

t0 = time.time()

def sparse_potts_hamiltonian(n, q, g, boundary='open'):
    """H = -J sum delta(s_i,s_j) - g sum (X_i + X_i^dag)"""
    dim = q**n
    # Cyclic permutation X
    X = np.zeros((q, q))
    for a in range(q):
        X[(a+1) % q, a] = 1.0
    X = csr_matrix(X)
    Xd = X.T.tocsr()
    XXd = csr_matrix(X + Xd)

    # Delta interaction: sum_a |a><a| x |a><a|
    delta = csr_matrix((q*q, q*q))
    for a in range(q):
        idx = a * q + a
        delta = delta + csr_matrix(([1.0], ([idx], [idx])), shape=(q*q, q*q))

    H = csr_matrix((dim, dim))

    # Bond terms: -J delta(s_i, s_{i+1})
    n_bonds = n - 1 if boundary == 'open' else n
    for bond in range(n_bonds):
        j = (bond + 1) % n
        if j == bond + 1:
            prefix = sp_eye(q**bond) if bond > 0 else sp_eye(1)
            suffix = sp_eye(q**(n - bond - 2)) if bond + 2 < n else sp_eye(1)
            H = H - sp_kron(sp_kron(prefix, delta, format='csr'), suffix, format='csr')
        else:
            for a in range(q):
                proj = csr_matrix(([1.0], ([a], [a])), shape=(q, q))
                middle = sp_eye(q**(n-2))
                H = H - sp_kron(sp_kron(proj, middle, format='csr'), proj, format='csr')

    # Site terms: -g (X + X^dag)
    for site in range(n):
        prefix = sp_eye(q**site) if site > 0 else sp_eye(1)
        suffix = sp_eye(q**(n - site - 1)) if site < n - 1 else sp_eye(1)
        H = H - g * sp_kron(sp_kron(prefix, XXd, format='csr'), suffix, format='csr')

    return H

def partial_trace_via_reshape(psi, n, q, keep_sites):
    """Efficient partial trace: reshape state vector, no full rho needed."""
    dims = [q] * n
    n_A = len(keep_sites)
    trace_sites = [i for i in range(n) if i not in keep_sites]

    # Reshape psi into tensor
    psi_tensor = psi.reshape(dims)

    # Move keep sites to front, trace sites to back
    all_sites = list(keep_sites) + list(trace_sites)
    psi_tensor = np.transpose(psi_tensor, all_sites)

    dim_A = q**n_A
    dim_B = q**(n - n_A)
    psi_matrix = psi_tensor.reshape(dim_A, dim_B)

    rho_A = psi_matrix @ psi_matrix.conj().T
    return rho_A

def entanglement_hamiltonian(rho_A):
    eigvals, eigvecs = la.eigh(rho_A)
    eigvals = np.maximum(eigvals, 1e-15)
    H_E = -eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.T
    return np.real(H_E)

def potts_subsystem_operators(n_A, q):
    """Build Potts-type operators on subsystem A: delta bonds + (X+X^dag) sites."""
    dim = q**n_A

    X_single = np.zeros((q, q))
    for a in range(q):
        X_single[(a+1) % q, a] = 1.0
    XXd_single = X_single + X_single.T

    delta_pair = np.zeros((q*q, q*q))
    for a in range(q):
        idx = a * q + a
        delta_pair[idx, idx] = 1.0

    ops = {}

    # Bond operators: delta(s_i, s_{i+1})
    for i in range(n_A - 1):
        prefix = np.eye(q**i) if i > 0 else np.eye(1)
        suffix = np.eye(q**(n_A - i - 2)) if i + 2 < n_A else np.eye(1)
        op = np.kron(np.kron(prefix, delta_pair), suffix)
        ops[f'delta_{i}_{i+1}'] = op

    # Site operators: X_i + X_i^dag
    for i in range(n_A):
        prefix = np.eye(q**i) if i > 0 else np.eye(1)
        suffix = np.eye(q**(n_A - i - 1)) if i < n_A - 1 else np.eye(1)
        op = np.kron(np.kron(prefix, XXd_single), suffix)
        ops[f'XXd_{i}'] = op

    return ops

def compute_potts_locality(H_E, n_A, q):
    """Decompose H_E into Potts-type terms via least squares. Return locality fraction."""
    dim = q**n_A
    H_E_tl = H_E - np.trace(H_E) / dim * np.eye(dim)
    total_norm_sq = la.norm(H_E_tl, 'fro')**2
    if total_norm_sq < 1e-15:
        return 0.0, {}, 0.0

    ops = potts_subsystem_operators(n_A, q)
    labels = list(ops.keys())

    op_vecs = []
    for label in labels:
        O = ops[label]
        O_tl = O - np.trace(O) / dim * np.eye(dim)
        op_vecs.append(O_tl.ravel())

    A_mat = np.array(op_vecs).T
    b_vec = H_E_tl.ravel()

    x_opt, _, _, _ = la.lstsq(A_mat, b_vec)

    potts_fit = A_mat @ x_opt
    fit_norm_sq = np.dot(potts_fit, potts_fit)
    locality = float(np.real(fit_norm_sq / np.dot(b_vec, b_vec)))

    coeffs = {labels[i]: float(np.real(x_opt[i])) for i in range(len(labels))}
    return locality, coeffs, total_norm_sq

def bw_hamiltonian(n_A, q, g, envelope_type='sin_inv'):
    """BW prediction: H_BW = sum beta(i) h_i with position-dependent temperature."""
    dim = q**n_A
    L_total = 2 * n_A

    def get_beta(pos, is_bond=True):
        if is_bond:
            dist = n_A - 1 - pos
        else:
            dist = n_A - 0.5 - pos
        if envelope_type == 'sin_inv':
            return np.sin(np.pi * (dist + 0.5) / L_total)
        elif envelope_type == 'linear':
            return dist + 0.5
        elif envelope_type == 'uniform':
            return 1.0

    ops = potts_subsystem_operators(n_A, q)
    H_BW = np.zeros((dim, dim))

    for i in range(n_A - 1):
        beta = get_beta(i, is_bond=True)
        O = ops[f'delta_{i}_{i+1}']
        O_tl = O - np.trace(O) / dim * np.eye(dim)
        H_BW -= beta * O_tl

    for i in range(n_A):
        beta = get_beta(i, is_bond=False)
        O = ops[f'XXd_{i}']
        O_tl = O - np.trace(O) / dim * np.eye(dim)
        H_BW -= g * beta * O_tl

    return H_BW

def compute_bw_metrics(H_E, H_BW):
    dim = H_E.shape[0]
    H_E_tl = H_E - np.trace(H_E) / dim * np.eye(dim)
    H_BW_tl = H_BW - np.trace(H_BW) / dim * np.eye(dim)

    denom = np.trace(H_BW_tl.T @ H_BW_tl)
    if denom < 1e-15:
        return 0.0, 0.0, 0.0
    alpha = np.trace(H_E_tl.T @ H_BW_tl) / denom
    residual = H_E_tl - float(alpha) * H_BW_tl
    resid_frac = la.norm(residual, 'fro') / la.norm(H_E_tl, 'fro')
    var_captured = 1 - resid_frac**2

    fidelity = np.real(np.trace(H_E_tl.T @ H_BW_tl)) / (la.norm(H_E_tl, 'fro') * la.norm(H_BW_tl, 'fro'))

    return float(fidelity), float(var_captured), float(np.real(alpha))

# ---- Main ----
q = 4
n = 8
n_A = n // 2
g_c = 0.392  # from energy gap method (Sprint 051)

print(f"=== Sprint 064a: BW Locality for q={q} Potts at g_c={g_c} ===")
print(f"n={n}, n_A={n_A}, dim_full={q**n}, dim_A={q**n_A}")

# Sweep g values through the transition
g_values = np.array([0.05, 0.10, 0.20, 0.30, 0.35, 0.38, 0.39, 0.392, 0.40, 0.42, 0.45, 0.50, 0.60, 0.80, 1.0])

results_list = []
envelope_types = ['sin_inv', 'linear', 'uniform']

for g in g_values:
    t1 = time.time()

    H = sparse_potts_hamiltonian(n, q, g, 'open')
    eigvals_H, eigvecs_H = eigsh(H, k=2, which='SA')
    sort_idx = np.argsort(eigvals_H)
    psi_gs = eigvecs_H[:, sort_idx[0]]
    gap = eigvals_H[sort_idx[1]] - eigvals_H[sort_idx[0]]

    # Partial trace via reshape (efficient)
    rho_A = partial_trace_via_reshape(psi_gs, n, q, list(range(n_A)))
    rho_A = (rho_A + rho_A.T) / 2

    # Entropy
    ev = la.eigvalsh(rho_A)
    ev_sig = ev[ev > 1e-15]
    S = -np.sum(ev_sig * np.log2(ev_sig))

    # Entanglement Hamiltonian
    H_E = entanglement_hamiltonian(rho_A)

    # Potts locality
    locality, coeffs, total_norm_sq = compute_potts_locality(H_E, n_A, q)

    # BW envelope comparison
    env_results = {}
    for etype in envelope_types:
        H_BW = bw_hamiltonian(n_A, q, g, etype)
        fid, var_cap, alpha = compute_bw_metrics(H_E, H_BW)
        env_results[etype] = {'fidelity': float(fid), 'variance_captured': float(var_cap), 'alpha': float(alpha)}

    best_env = max(envelope_types, key=lambda e: env_results[e]['variance_captured'])
    dt = time.time() - t1

    result = {
        'g': float(g), 'entropy_bits': float(S), 'gap': float(gap),
        'potts_locality': float(locality), 'coefficients': coeffs,
        'envelope_results': env_results, 'best_envelope': best_env,
        'time_s': float(dt)
    }
    results_list.append(result)

    loc_pct = locality * 100
    best_var = env_results[best_env]['variance_captured'] * 100
    print(f"  g={g:.3f}: S={S:.3f} bits, gap={gap:.4f}, locality={loc_pct:.1f}%, "
          f"BW[{best_env}]={best_var:.1f}% [{dt:.1f}s]")

# ---- Summary ----
print("\n=== Summary ===")
print("  g      Locality%  BW_sin%  BW_lin%  BW_uni%  best")
for r in results_list:
    loc = r['potts_locality'] * 100
    si = r['envelope_results']['sin_inv']['variance_captured'] * 100
    li = r['envelope_results']['linear']['variance_captured'] * 100
    un = r['envelope_results']['uniform']['variance_captured'] * 100
    print(f"  {r['g']:5.3f}  {loc:7.1f}     {si:5.1f}    {li:5.1f}    {un:5.1f}    {r['best_envelope']}")

# Peak locality
peak = max(results_list, key=lambda r: r['potts_locality'])
print(f"\nPeak locality: {peak['potts_locality']*100:.1f}% at g={peak['g']:.3f}")

# Coupling profile at g_c
gc_result = [r for r in results_list if abs(r['g'] - g_c) < 0.001][0]
print(f"\nCoupling profile at g_c={g_c}:")
for label, c in sorted(gc_result['coefficients'].items()):
    if abs(c) > 1e-6:
        print(f"  {label}: {c:+.6f}")

# Compare to q=3
print(f"\n=== Comparison ===")
print(f"q=3 Potts (S_3): peak BW locality = 76.5% (Sprint 033)")
print(f"q=4 Potts (S_4?): peak BW locality = {peak['potts_locality']*100:.1f}%")

# ---- Save ----
results = {
    'experiment': '064a',
    'description': f'BW locality for q={q} Potts at g_c={g_c}',
    'parameters': {'n': n, 'q': q, 'n_A': n_A, 'g_c': g_c},
    'sweep': results_list,
    'peak_locality': {'g': peak['g'], 'locality': peak['potts_locality']},
    'runtime_seconds': time.time() - t0
}

with open('results/sprint_064a_bw_q4.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nSaved. Runtime: {time.time()-t0:.1f}s")
