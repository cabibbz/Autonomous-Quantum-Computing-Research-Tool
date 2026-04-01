"""
Sprint 033c: BW Locality for 3-State Potts Model — S₃ vs Z₂ vs U(1)

Test Sprint 032's key prediction: symmetry group controls BW accuracy.
- Z₂ (TFIM): 91% BW locality
- S₃ (Potts q=3): ??? (prediction: intermediate, 91-100%)
- U(1) (XXZ): 100% BW locality

The BW prediction for Potts: H_E = Σ β(i) h_i where h_i are local Potts terms
(δ interaction + P+P† transverse field) with position-dependent β(i).

For the Potts model with 3-state local Hilbert space, we decompose H_E into
Potts-type terms (δ bonds and P+P† sites) vs non-Potts corrections.

System: n=8, subsystem A = left 4 sites (dim_A = 3^4 = 81).
"""

import numpy as np
from scipy import linalg as la
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import eigsh
import json, time

t0 = time.time()

# ---- Hamiltonian construction ----
def sparse_potts_hamiltonian(n, h_over_J, boundary='open'):
    d = 3
    dim = d**n
    P = csr_matrix(np.array([[0,0,1],[1,0,0],[0,1,0]], dtype=float))
    Pd = P.T.tocsr()
    delta = csr_matrix((d*d, d*d))
    for a in range(d):
        idx = a * d + a
        delta = delta + csr_matrix(([1.0], ([idx], [idx])), shape=(d*d, d*d))
    H = csr_matrix((dim, dim))
    n_bonds = n - 1 if boundary == 'open' else n
    for bond in range(n_bonds):
        j = (bond + 1) % n
        if j == bond + 1:
            prefix = sp_eye(d**bond, format='csr') if bond > 0 else sp_eye(1, format='csr')
            suffix = sp_eye(d**(n - bond - 2), format='csr') if bond + 2 < n else sp_eye(1, format='csr')
            H = H - sp_kron(sp_kron(prefix, delta), suffix, format='csr')
        else:
            for a in range(d):
                proj = csr_matrix(([1.0], ([a], [a])), shape=(d, d))
                middle = sp_eye(d**(n-2), format='csr')
                H = H - sp_kron(sp_kron(proj, middle), proj, format='csr')
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

def entanglement_hamiltonian(rho_A):
    eigvals, eigvecs = la.eigh(rho_A)
    eigvals = np.maximum(eigvals, 1e-15)
    H_E = -eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.T
    return np.real(H_E)

# ---- Potts operator basis for subsystem A ----
def potts_subsystem_operators(n_A):
    """
    Build all Potts-type operators on subsystem A (n_A sites, local dim 3):
    - δ(σ_i, σ_{i+1}) for adjacent bonds (interaction terms)
    - P_i + P_i† for each site (transverse field terms)

    Returns dict of {label: operator_matrix}.
    """
    d = 3
    dim = d**n_A

    P_single = np.array([[0,0,1],[1,0,0],[0,1,0]], dtype=float)
    Pd_single = P_single.T
    PPd_single = P_single + Pd_single
    I_d = np.eye(d)

    # Delta for a pair
    delta_pair = np.zeros((d*d, d*d))
    for a in range(d):
        idx = a * d + a
        delta_pair[idx, idx] = 1.0

    ops = {}

    # Bond operators: δ(σ_i, σ_{i+1})
    for i in range(n_A - 1):
        prefix = np.eye(d**i) if i > 0 else np.eye(1)
        suffix = np.eye(d**(n_A - i - 2)) if i + 2 < n_A else np.eye(1)
        op = np.kron(np.kron(prefix, delta_pair), suffix)
        ops[f'delta_{i}_{i+1}'] = op

    # Site operators: P_i + P_i†
    for i in range(n_A):
        prefix = np.eye(d**i) if i > 0 else np.eye(1)
        suffix = np.eye(d**(n_A - i - 1)) if i < n_A - 1 else np.eye(1)
        op = np.kron(np.kron(prefix, PPd_single), suffix)
        ops[f'PPd_{i}'] = op

    return ops

def compute_potts_locality(H_E, n_A):
    """
    Decompose H_E into Potts-type terms and compute locality fraction.

    Potts-type terms: δ(σ_i,σ_{i+1}) bonds and (P_i + P_i†) sites.
    Locality = fraction of H_E norm² in Potts-type terms.
    """
    d = 3
    dim = d**n_A
    H_E_tl = H_E - np.trace(H_E) / dim * np.eye(dim)
    total_norm_sq = la.norm(H_E_tl, 'fro')**2
    if total_norm_sq < 1e-15:
        return 0.0, {}, 0.0

    ops = potts_subsystem_operators(n_A)

    # Project H_E onto each operator: c_i = Tr(H_E_tl · O_i) / Tr(O_i · O_i)
    coeffs = {}
    potts_reconstruction = np.zeros_like(H_E_tl)

    for label, O in ops.items():
        O_tl = O - np.trace(O) / dim * np.eye(dim)
        norm_sq_O = la.norm(O_tl, 'fro')**2
        if norm_sq_O < 1e-15:
            coeffs[label] = 0.0
            continue
        c = np.trace(H_E_tl @ O_tl) / norm_sq_O
        coeffs[label] = float(np.real(c))
        potts_reconstruction += float(np.real(c)) * O_tl

    potts_norm_sq = la.norm(potts_reconstruction, 'fro')**2
    locality = potts_norm_sq / total_norm_sq

    # Also compute what fraction using best-fit (Gram-Schmidt-like orthogonal projection)
    # The operators might not be orthogonal, so do a proper least-squares fit
    n_ops = len(ops)
    labels = list(ops.keys())
    op_matrices = []
    for label in labels:
        O = ops[label]
        O_tl = O - np.trace(O) / dim * np.eye(dim)
        op_matrices.append(O_tl.ravel())

    A_mat = np.array(op_matrices).T  # dim² × n_ops
    b_vec = H_E_tl.ravel()

    # Least squares: minimize ||b - A·x||²
    x_opt, residual_arr, _, _ = la.lstsq(A_mat, b_vec)

    potts_fit = A_mat @ x_opt
    fit_norm_sq = np.dot(potts_fit, potts_fit)
    locality_lstsq = float(np.real(fit_norm_sq / np.dot(b_vec, b_vec)))

    coeffs_lstsq = {labels[i]: float(np.real(x_opt[i])) for i in range(n_ops)}

    return float(locality_lstsq), coeffs_lstsq, total_norm_sq

def bw_hamiltonian_potts(n_A, h_over_J, envelope_type='sin_inv'):
    """BW prediction: H_BW = Σ β(i) h_i with position-dependent temperature."""
    d = 3
    dim = d**n_A
    L_total = 2 * n_A

    def get_beta(pos, is_bond=True):
        """β as function of distance from entanglement cut."""
        if is_bond:
            dist = n_A - 1 - pos  # distance of bond from cut
        else:
            dist = n_A - 0.5 - pos  # distance of site from cut

        if envelope_type == 'sin_inv':
            return np.sin(np.pi * (dist + 0.5) / L_total)
        elif envelope_type == 'linear':
            return dist + 0.5
        elif envelope_type == 'uniform':
            return 1.0

    ops = potts_subsystem_operators(n_A)
    H_BW = np.zeros((dim, dim))

    # Add Potts terms with BW envelope
    for i in range(n_A - 1):
        beta = get_beta(i, is_bond=True)
        label = f'delta_{i}_{i+1}'
        O = ops[label] - np.trace(ops[label]) / dim * np.eye(dim)
        H_BW -= beta * O  # -J δ → coefficient is negative

    for i in range(n_A):
        beta = get_beta(i, is_bond=False)
        label = f'PPd_{i}'
        O = ops[label] - np.trace(ops[label]) / dim * np.eye(dim)
        H_BW -= h_over_J * beta * O  # -h (P+P†)

    return H_BW

def matrix_overlap(A, B):
    nA = la.norm(A, 'fro')
    nB = la.norm(B, 'fro')
    if nA < 1e-15 or nB < 1e-15:
        return 0.0
    return np.real(np.trace(A.T @ B)) / (nA * nB)

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
    [0.05, 0.10],                 # deep ordered
    np.linspace(0.15, 0.35, 9),   # transition region
    [0.40, 0.50, 0.75, 1.0, 2.0]  # disordered
])

print(f"=== Sprint 033c: BW Locality for 3-State Potts (n={n}, n_A={n_A}) ===")
print(f"dim_A = {dim_A}, sweeping {len(h_values)} h/J values")

results_list = []
envelope_types = ['sin_inv', 'linear', 'uniform']

for h_J in h_values:
    t1 = time.time()

    H = sparse_potts_hamiltonian(n, h_J, 'open')
    eigvals_H, eigvecs_H = eigsh(H, k=4, which='SA')
    sort_idx = np.argsort(eigvals_H)
    psi_gs = eigvecs_H[:, sort_idx[0]]

    rho = np.outer(psi_gs, psi_gs.conj())
    rho_A = partial_trace_general(rho, dims, list(range(n_A)))
    rho_A = (rho_A + rho_A.T) / 2  # symmetrize

    # Entropy
    ev = la.eigvalsh(rho_A)
    ev_sig = ev[ev > 1e-15]
    S = -np.sum(ev_sig * np.log2(ev_sig))

    # Entanglement Hamiltonian
    H_E = entanglement_hamiltonian(rho_A)

    # Potts locality (fraction of H_E in Potts-type terms)
    locality, coeffs, total_norm_sq = compute_potts_locality(H_E, n_A)

    # BW envelope comparison
    env_results = {}
    for etype in envelope_types:
        H_BW = bw_hamiltonian_potts(n_A, h_J, etype)
        fid, var_cap, alpha = compute_bw_metrics(H_E, H_BW)
        env_results[etype] = {'fidelity': fid, 'variance_captured': var_cap, 'alpha': alpha}

    best_env = max(envelope_types, key=lambda e: env_results[e]['variance_captured'])

    dt = time.time() - t1

    result = {
        'h_J': float(h_J),
        'entropy_bits': float(S),
        'potts_locality': float(locality),
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
    loc = r['potts_locality'] * 100
    si = r['envelope_results']['sin_inv']['variance_captured'] * 100
    li = r['envelope_results']['linear']['variance_captured'] * 100
    un = r['envelope_results']['uniform']['variance_captured'] * 100
    print(f"  {r['h_J']:5.3f}  {loc:7.1f}     {si:7.1f}       {li:7.1f}       {un:7.1f}      {r['best_envelope']}")

# ---- Key comparison: S₃ vs Z₂ vs U(1) ----
print("\n=== SYMMETRY COMPARISON ===")
print("Model        Symmetry  |G|  BW Locality at Peak")
print("TFIM         Z₂        2    91% (Sprint 032a)")
print("XXZ          U(1)      ∞    100% (Sprint 032c)")

# Find peak locality for Potts
peak = max(results_list, key=lambda r: r['potts_locality'])
print(f"Potts (q=3)  S₃        6    {peak['potts_locality']*100:.1f}% at h/J={peak['h_J']:.3f}")

# Ordered phase locality
ordered = [r for r in results_list if r['h_J'] <= 0.10]
if ordered:
    ord_loc = np.mean([r['potts_locality'] for r in ordered])
    print(f"\nOrdered phase average locality: {ord_loc*100:.1f}%")

# Disordered phase locality
disordered = [r for r in results_list if r['h_J'] >= 0.75]
if disordered:
    dis_loc = np.mean([r['potts_locality'] for r in disordered])
    print(f"Disordered phase average locality: {dis_loc*100:.1f}%")

# ---- Extract coupling profiles ----
print("\n=== Coupling Profiles at Selected Points ===")
for r in results_list:
    if r['h_J'] in [0.10, 0.225, 0.50, 1.0]:
        print(f"\n  h/J = {r['h_J']:.3f} (locality={r['potts_locality']*100:.1f}%):")
        for label, c in sorted(r['coefficients'].items()):
            if abs(c) > 1e-6:
                print(f"    {label}: {c:+.6f}")

# ---- Save ----
results = {
    'experiment': '033c',
    'description': 'BW locality for 3-state Potts model',
    'parameters': {'n': n, 'local_dim': d, 'n_A': n_A},
    'sweep': results_list,
    'symmetry_comparison': {
        'TFIM_Z2': {'group_order': 2, 'peak_locality': 0.91},
        'Potts_S3': {'group_order': 6, 'peak_locality': peak['potts_locality']},
        'XXZ_U1': {'group_order': float('inf'), 'peak_locality': 1.0},
    },
    'runtime_seconds': time.time() - t0
}

with open('results/sprint_033c_potts_bw.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nSaved. Runtime: {time.time()-t0:.1f}s")
