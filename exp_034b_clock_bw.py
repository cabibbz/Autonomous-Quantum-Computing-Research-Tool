"""
Sprint 034b: BW Locality for Z₃ Clock Model — Isolating Group Size from Local Dimension

The Z₃ clock model on qutrits (d=3):
  H = -J Σ (Z_i Z_{i+1}† + h.c.) - h Σ (X_i + X_i†)
where Z = diag(1, ω, ω²), X = cyclic permutation |k⟩ → |k+1 mod 3⟩.

This has Z₃ symmetry (order 3) at d=3, vs Potts which has S₃ symmetry (order 6) at d=3.
From 034a: Z₃ has MORE invariant operators (125 vs 69 at n_A=3), so:
  PREDICTION: Z₃ BW locality < S₃ BW locality (76.5%)

If confirmed: group size matters (within fixed d), but d dominates across models.
"""

import numpy as np
from scipy import linalg as la
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import eigsh
import json, time

t0 = time.time()

# ---- Z₃ Clock Hamiltonian ----
omega = np.exp(2j * np.pi / 3)

def clock_Z():
    """Z operator: Z|k⟩ = ω^k |k⟩"""
    return np.diag([1.0, omega, omega**2])

def clock_X():
    """X operator: X|k⟩ = |k+1 mod 3⟩"""
    return np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=complex)

def sparse_clock_hamiltonian(n, h_over_J, boundary='open'):
    """
    H = -J Σ (Z_i Z_{i+1}† + Z_i† Z_{i+1}) - h Σ (X_i + X_i†)
    """
    d = 3
    dim = d**n

    Z = clock_Z()
    Zd = Z.conj().T  # Z†
    X = clock_X()
    Xd = X.conj().T

    # ZZ† + Z†Z bond operator (Hermitian)
    ZZd = np.kron(Z, Zd) + np.kron(Zd, Z)
    ZZd_sparse = csr_matrix(ZZd)

    # X + X† site operator (Hermitian)
    XXd = X + Xd
    XXd_sparse = csr_matrix(XXd)

    H = csr_matrix((dim, dim), dtype=complex)

    # Bond terms
    n_bonds = n - 1 if boundary == 'open' else n
    for bond in range(n_bonds):
        j = (bond + 1) % n
        if j == bond + 1:
            prefix = sp_eye(d**bond, format='csr') if bond > 0 else sp_eye(1, format='csr')
            suffix = sp_eye(d**(n - bond - 2), format='csr') if bond + 2 < n else sp_eye(1, format='csr')
            H = H - sp_kron(sp_kron(prefix, ZZd_sparse), suffix, format='csr')

    # Site terms
    for site in range(n):
        prefix = sp_eye(d**site, format='csr') if site > 0 else sp_eye(1, format='csr')
        suffix = sp_eye(d**(n - site - 1), format='csr') if site < n - 1 else sp_eye(1, format='csr')
        H = H - h_over_J * sp_kron(sp_kron(prefix, XXd_sparse), suffix, format='csr')

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

def entanglement_hamiltonian(rho_A):
    eigvals, eigvecs = la.eigh(rho_A)
    eigvals = np.maximum(eigvals, 1e-15)
    H_E = -eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.T
    return np.real(H_E)

# ---- Clock operator basis for subsystem A ----
def clock_subsystem_operators(n_A):
    """
    Build all clock-type operators on subsystem A (n_A sites, local dim 3):
    - (Z_i Z_{i+1}† + h.c.) for adjacent bonds
    - (X_i + X_i†) for each site
    """
    d = 3
    dim = d**n_A

    Z = clock_Z()
    Zd = Z.conj().T
    X = clock_X()
    Xd = X.conj().T
    I_d = np.eye(d, dtype=complex)

    ZZd_pair = np.kron(Z, Zd) + np.kron(Zd, Z)  # Hermitian
    XXd_single = X + Xd  # Hermitian

    ops = {}

    # Bond operators: Z_i Z_{i+1}† + h.c.
    for i in range(n_A - 1):
        prefix = np.eye(d**i, dtype=complex) if i > 0 else np.eye(1, dtype=complex)
        suffix = np.eye(d**(n_A - i - 2), dtype=complex) if i + 2 < n_A else np.eye(1, dtype=complex)
        op = np.kron(np.kron(prefix, ZZd_pair), suffix)
        ops[f'ZZd_{i}_{i+1}'] = np.real(op)  # Should be real for Hermitian

    # Site operators: X + X†
    for i in range(n_A):
        prefix = np.eye(d**i, dtype=complex) if i > 0 else np.eye(1, dtype=complex)
        suffix = np.eye(d**(n_A - i - 1), dtype=complex) if i < n_A - 1 else np.eye(1, dtype=complex)
        op = np.kron(np.kron(prefix, XXd_single), suffix)
        ops[f'XXd_{i}'] = np.real(op)

    return ops

def compute_clock_locality(H_E, n_A):
    """
    Decompose H_E into clock-type terms and compute locality fraction.
    Uses least-squares fit (proper projection onto non-orthogonal basis).
    """
    d = 3
    dim = d**n_A
    H_E_tl = H_E - np.trace(H_E) / dim * np.eye(dim)
    total_norm_sq = la.norm(H_E_tl, 'fro')**2
    if total_norm_sq < 1e-15:
        return 0.0, {}, 0.0

    ops = clock_subsystem_operators(n_A)
    labels = list(ops.keys())
    op_matrices = []
    for label in labels:
        O = ops[label]
        O_tl = O - np.trace(O) / dim * np.eye(dim)
        op_matrices.append(O_tl.ravel())

    A_mat = np.array(op_matrices).T
    b_vec = H_E_tl.ravel()

    x_opt, _, _, _ = la.lstsq(A_mat, b_vec)

    clock_fit = A_mat @ x_opt
    fit_norm_sq = np.dot(clock_fit, clock_fit)
    locality = float(np.real(fit_norm_sq / np.dot(b_vec, b_vec)))

    coeffs = {labels[i]: float(np.real(x_opt[i])) for i in range(len(labels))}

    return locality, coeffs, total_norm_sq

def matrix_overlap(A, B):
    nA = la.norm(A, 'fro')
    nB = la.norm(B, 'fro')
    if nA < 1e-15 or nB < 1e-15:
        return 0.0
    return np.real(np.trace(A.T @ B)) / (nA * nB)

def bw_hamiltonian_clock(n_A, h_over_J, envelope_type='sin_inv'):
    """BW prediction: H_BW = Σ β(i) h_i with position-dependent temperature."""
    d = 3
    dim = d**n_A
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

    ops = clock_subsystem_operators(n_A)
    H_BW = np.zeros((dim, dim))

    for i in range(n_A - 1):
        beta = get_beta(i, is_bond=True)
        label = f'ZZd_{i}_{i+1}'
        O = ops[label] - np.trace(ops[label]) / dim * np.eye(dim)
        H_BW -= beta * O

    for i in range(n_A):
        beta = get_beta(i, is_bond=False)
        label = f'XXd_{i}'
        O = ops[label] - np.trace(ops[label]) / dim * np.eye(dim)
        H_BW -= h_over_J * beta * O

    return H_BW

def compute_bw_metrics(H_E, H_BW):
    dim = H_E.shape[0]
    H_E_tl = H_E - np.trace(H_E) / dim * np.eye(dim)
    H_BW_tl = H_BW - np.trace(H_BW) / dim * np.eye(dim)

    fidelity = matrix_overlap(H_E_tl, H_BW_tl)

    denom = np.trace(H_BW_tl.T @ H_BW_tl)
    if denom < 1e-15:
        return 0.0, 0.0, 0.0
    alpha = np.trace(H_E_tl.T @ H_BW_tl) / denom
    residual = H_E_tl - float(alpha) * H_BW_tl
    resid_frac = la.norm(residual, 'fro') / la.norm(H_E_tl, 'fro')
    var_captured = 1 - resid_frac**2

    return float(fidelity), float(var_captured), float(alpha)

# ---- Main sweep ----
n = 8
d = 3
dims = [d] * n
n_A = n // 2
dim_A = d**n_A  # 81

# Sweep h/J through the transition
h_values = np.concatenate([
    [0.05, 0.10],
    np.linspace(0.15, 0.40, 11),
    [0.50, 0.75, 1.0, 1.5, 2.0]
])

print(f"=== Sprint 034b: BW Locality for Z₃ Clock Model (n={n}, n_A={n_A}) ===")
print(f"dim = {d**n}, dim_A = {dim_A}, sweeping {len(h_values)} h/J values")

results_list = []
envelope_types = ['sin_inv', 'linear', 'uniform']

for h_J in h_values:
    t1 = time.time()

    H = sparse_clock_hamiltonian(n, h_J, 'open')
    # Make sure H is Hermitian
    H = (H + H.conj().T) / 2

    eigvals_H, eigvecs_H = eigsh(H, k=4, which='SA')
    sort_idx = np.argsort(eigvals_H)
    psi_gs = eigvecs_H[:, sort_idx[0]]

    rho = np.outer(psi_gs, psi_gs.conj())
    rho = np.real(rho)  # Ground state should be real for this Hermitian H
    rho_A = partial_trace_general(rho, dims, list(range(n_A)))
    rho_A = (rho_A + rho_A.T) / 2

    # Entropy
    ev = la.eigvalsh(rho_A)
    ev_sig = ev[ev > 1e-15]
    S = -np.sum(ev_sig * np.log2(ev_sig))

    # Entanglement Hamiltonian
    H_E = entanglement_hamiltonian(rho_A)

    # Clock locality
    locality, coeffs, total_norm_sq = compute_clock_locality(H_E, n_A)

    # BW envelope comparison
    env_results = {}
    for etype in envelope_types:
        H_BW = bw_hamiltonian_clock(n_A, h_J, etype)
        fid, var_cap, alpha = compute_bw_metrics(H_E, H_BW)
        env_results[etype] = {'fidelity': fid, 'variance_captured': var_cap, 'alpha': alpha}

    best_env = max(envelope_types, key=lambda e: env_results[e]['variance_captured'])

    dt = time.time() - t1

    result = {
        'h_J': float(h_J),
        'entropy_bits': float(S),
        'clock_locality': float(locality),
        'coefficients': coeffs,
        'envelope_results': env_results,
        'best_envelope': best_env,
        'time_s': float(dt)
    }
    results_list.append(result)

    loc_pct = locality * 100
    best_var = env_results[best_env]['variance_captured'] * 100
    print(f"  h/J={h_J:.3f}: S={S:.3f}, locality={loc_pct:.1f}%, "
          f"BW[{best_env}]={best_var:.1f}% [{dt:.1f}s]")

# ---- Summary ----
print("\n=== BW Locality Summary ===")
print("  h/J    Locality%  BW_sin_inv%  BW_linear%  BW_uniform%  best")
for r in results_list:
    loc = r['clock_locality'] * 100
    si = r['envelope_results']['sin_inv']['variance_captured'] * 100
    li = r['envelope_results']['linear']['variance_captured'] * 100
    un = r['envelope_results']['uniform']['variance_captured'] * 100
    print(f"  {r['h_J']:5.3f}  {loc:7.1f}     {si:7.1f}       {li:7.1f}       {un:7.1f}      {r['best_envelope']}")

# ---- Key comparison ----
print("\n=== SYMMETRY COMPARISON: ALL MODELS ===")
print(f"{'Model':<12} {'Group':<8} {'|G|':<5} {'d':<4} {'Peak BW Locality'}")
print(f"{'XXZ':<12} {'U(1)':<8} {'∞':<5} {'2':<4} 100.0%")
print(f"{'TFIM':<12} {'Z₂':<8} {'2':<5} {'2':<4} 91.0%")

peak = max(results_list, key=lambda r: r['clock_locality'])
print(f"{'Potts':<12} {'S₃':<8} {'6':<5} {'3':<4} 76.5%")
print(f"{'Clock':<12} {'Z₃':<8} {'3':<5} {'3':<4} {peak['clock_locality']*100:.1f}% at h/J={peak['h_J']:.3f}")

if peak['clock_locality'] < 0.765:
    print("\n✓ PREDICTION CONFIRMED: Z₃ < S₃ — smaller group → more invariant ops → lower BW locality")
else:
    print("\n✗ PREDICTION OVERTURNED: Z₃ ≥ S₃ — group size alone doesn't explain BW within fixed d")

# ---- Coupling profiles ----
print("\n=== Coupling Profiles at Selected Points ===")
for r in results_list:
    if abs(r['h_J'] - 0.10) < 0.01 or abs(r['h_J'] - peak['h_J']) < 0.01 or abs(r['h_J'] - 0.50) < 0.01 or abs(r['h_J'] - 1.0) < 0.01:
        print(f"\n  h/J = {r['h_J']:.3f} (locality={r['clock_locality']*100:.1f}%):")
        for label, c in sorted(r['coefficients'].items()):
            if abs(c) > 1e-6:
                print(f"    {label}: {c:+.6f}")

# ---- Save ----
results = {
    'experiment': '034b',
    'description': 'BW locality for Z₃ clock model',
    'parameters': {'n': n, 'local_dim': d, 'n_A': n_A},
    'sweep': results_list,
    'peak_locality': peak['clock_locality'],
    'peak_h_J': peak['h_J'],
    'symmetry_comparison': {
        'TFIM_Z2': {'group_order': 2, 'd': 2, 'peak_locality': 0.91},
        'Potts_S3': {'group_order': 6, 'd': 3, 'peak_locality': 0.765},
        'Clock_Z3': {'group_order': 3, 'd': 3, 'peak_locality': peak['clock_locality']},
        'XXZ_U1': {'group_order': float('inf'), 'd': 2, 'peak_locality': 1.0},
    },
    'runtime_seconds': time.time() - t0
}

with open('results/sprint_034b_clock_bw.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nSaved. Runtime: {time.time()-t0:.1f}s")
