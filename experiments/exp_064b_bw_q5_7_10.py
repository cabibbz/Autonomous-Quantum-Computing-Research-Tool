"""
Sprint 064b: BW Locality for q=5,7,10 Potts at True Critical Points

Extends 064a to the novel CFT regime (q>4).
- q=5: n=8, dim=390625, dim_A=5^4=625. g_c=0.441
- q=7: n=6, dim=117649, dim_A=7^3=343. g_c=0.535
- q=10: n=6, dim=1000000, dim_A=10^3=1000. g_c=0.684

For each q, sweep g through transition and measure BW locality + coupling profiles.
"""

import numpy as np
from scipy import linalg as la
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import eigsh
import json, time

t0 = time.time()

def sparse_potts_hamiltonian(n, q, g, boundary='open'):
    dim = q**n
    X = np.zeros((q, q))
    for a in range(q):
        X[(a+1) % q, a] = 1.0
    X = csr_matrix(X)
    XXd = csr_matrix(X + X.T)
    delta = csr_matrix((q*q, q*q))
    for a in range(q):
        idx = a * q + a
        delta = delta + csr_matrix(([1.0], ([idx], [idx])), shape=(q*q, q*q))
    H = csr_matrix((dim, dim))
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
    for site in range(n):
        prefix = sp_eye(q**site) if site > 0 else sp_eye(1)
        suffix = sp_eye(q**(n - site - 1)) if site < n - 1 else sp_eye(1)
        H = H - g * sp_kron(sp_kron(prefix, XXd, format='csr'), suffix, format='csr')
    return H

def partial_trace_via_reshape(psi, n, q, keep_sites):
    dims = [q] * n
    n_A = len(keep_sites)
    trace_sites = [i for i in range(n) if i not in keep_sites]
    psi_tensor = psi.reshape(dims)
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
    for i in range(n_A - 1):
        prefix = np.eye(q**i) if i > 0 else np.eye(1)
        suffix = np.eye(q**(n_A - i - 2)) if i + 2 < n_A else np.eye(1)
        ops[f'delta_{i}_{i+1}'] = np.kron(np.kron(prefix, delta_pair), suffix)
    for i in range(n_A):
        prefix = np.eye(q**i) if i > 0 else np.eye(1)
        suffix = np.eye(q**(n_A - i - 1)) if i < n_A - 1 else np.eye(1)
        ops[f'XXd_{i}'] = np.kron(np.kron(prefix, XXd_single), suffix)
    return ops

def compute_potts_locality(H_E, n_A, q):
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

def bw_hamiltonian(n_A, q, g, envelope_type='linear'):
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

def run_bw_sweep(q, n, g_c, g_values):
    """Run BW locality sweep for given q, n, g_c."""
    n_A = n // 2
    print(f"\n{'='*60}")
    print(f"q={q}, n={n}, n_A={n_A}, dim_full={q**n}, dim_A={q**n_A}, g_c={g_c}")
    print(f"{'='*60}")

    results_list = []
    envelope_types = ['sin_inv', 'linear', 'uniform']

    for g in g_values:
        t1 = time.time()

        H = sparse_potts_hamiltonian(n, q, g, 'open')
        eigvals_H, eigvecs_H = eigsh(H, k=2, which='SA')
        sort_idx = np.argsort(eigvals_H)
        psi_gs = eigvecs_H[:, sort_idx[0]]
        gap = eigvals_H[sort_idx[1]] - eigvals_H[sort_idx[0]]

        rho_A = partial_trace_via_reshape(psi_gs, n, q, list(range(n_A)))
        rho_A = (rho_A + rho_A.T) / 2

        ev = la.eigvalsh(rho_A)
        ev_sig = ev[ev > 1e-15]
        S = -np.sum(ev_sig * np.log2(ev_sig))

        H_E = entanglement_hamiltonian(rho_A)
        locality, coeffs, total_norm_sq = compute_potts_locality(H_E, n_A, q)

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
        print(f"  g={g:.3f}: S={S:.3f}, locality={loc_pct:.1f}%, "
              f"BW[{best_env}]={best_var:.1f}% [{dt:.1f}s]")

    peak = max(results_list, key=lambda r: r['potts_locality'])
    print(f"  Peak locality: {peak['potts_locality']*100:.1f}% at g={peak['g']:.3f}")

    # Coupling profile at g_c
    gc_result = min(results_list, key=lambda r: abs(r['g'] - g_c))
    print(f"  Coupling profile at g~{gc_result['g']:.3f}:")
    for label, c in sorted(gc_result['coefficients'].items()):
        if abs(c) > 1e-6:
            print(f"    {label}: {c:+.6f}")

    return results_list, peak

# ---- Run for q=5, 7, 10 ----
all_results = {}

# q=5: n=8 (dim=390625, ~10s per point)
g_vals_5 = np.array([0.05, 0.15, 0.30, 0.40, 0.441, 0.50, 0.60, 0.80, 1.0])
res5, peak5 = run_bw_sweep(5, 8, 0.441, g_vals_5)
all_results['q5'] = {'results': res5, 'peak': peak5, 'q': 5, 'n': 8, 'g_c': 0.441}

# q=7: n=6 (dim=117649, ~5s per point)
g_vals_7 = np.array([0.05, 0.20, 0.40, 0.50, 0.535, 0.60, 0.70, 0.90])
res7, peak7 = run_bw_sweep(7, 6, 0.535, g_vals_7)
all_results['q7'] = {'results': res7, 'peak': peak7, 'q': 7, 'n': 6, 'g_c': 0.535}

# q=10: n=6 (dim=1000000, ~60s per point)
g_vals_10 = np.array([0.05, 0.30, 0.55, 0.684, 0.80, 1.0])
res10, peak10 = run_bw_sweep(10, 6, 0.684, g_vals_10)
all_results['q10'] = {'results': res10, 'peak': peak10, 'q': 10, 'n': 6, 'g_c': 0.684}

# ---- Grand comparison ----
print(f"\n{'='*60}")
print("GRAND BW LOCALITY COMPARISON")
print(f"{'='*60}")
print(f"q=2 TFIM (Z_2, d=2):  peak = 91.0% (Sprint 032)")
print(f"q=2 XXZ  (U(1), d=2): peak = 100.0% (Sprint 032)")
print(f"q=3 Potts (S_3, d=3): peak = 76.5% (Sprint 033)")
print(f"q=4 Potts (Z_4, d=4): peak = 56.0% (Sprint 064a)")
print(f"q=5 Potts (Z_5, d=5): peak = {peak5['potts_locality']*100:.1f}%")
print(f"q=7 Potts (Z_7, d=7): peak = {peak7['potts_locality']*100:.1f}%")
print(f"q=10 Potts(Z_10,d=10):peak = {peak10['potts_locality']*100:.1f}%")

# ---- Save ----
save_data = {
    'experiment': '064b',
    'description': 'BW locality for q=5,7,10 Potts at true critical points',
    'q5': {'n': 8, 'g_c': 0.441, 'sweep': res5,
           'peak': {'g': peak5['g'], 'locality': peak5['potts_locality']}},
    'q7': {'n': 6, 'g_c': 0.535, 'sweep': res7,
           'peak': {'g': peak7['g'], 'locality': peak7['potts_locality']}},
    'q10': {'n': 6, 'g_c': 0.684, 'sweep': res10,
            'peak': {'g': peak10['g'], 'locality': peak10['potts_locality']}},
    'runtime_seconds': time.time() - t0
}

with open('results/sprint_064b_bw_q5_7_10.json', 'w') as f:
    json.dump(save_data, f, indent=2)

print(f"\nSaved. Total runtime: {time.time()-t0:.1f}s")
