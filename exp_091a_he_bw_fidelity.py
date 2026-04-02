#!/usr/bin/env python3
"""Sprint 091a: Entanglement Hamiltonian BW fidelity for S_q Potts across walking boundary.

Compute H_E = -log(rho_A) for periodic S_q Potts at g_c=1/q.
Project onto BW form (Potts NN operators with sin envelope).
Measure fidelity, variance captured, and coupling profiles.

q=2 n=12, q=3 n=10, q=4 n=8, q=5 n=8, q=7 n=6.
"""
import numpy as np
from scipy import linalg as la
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
import json, time

t0_global = time.time()

results = {
    'experiment': '091a_he_bw_fidelity',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_091a_he_bw_fidelity.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_periodic(n, q, g):
    """Build S_q Potts Hamiltonian on periodic chain.
    H = -sum_{<ij>} delta(s_i,s_j) - g * sum_i sum_{k=1}^{q-1} X_i^k
    """
    dim = q**n
    H = lil_matrix((dim, dim), dtype=complex)
    for idx in range(dim):
        c = []
        tmp = idx
        for _ in range(n):
            c.append(tmp % q)
            tmp //= q
        # Diagonal: Potts coupling (periodic)
        diag = 0.0
        for site in range(n):
            nxt = (site + 1) % n
            if c[site] == c[nxt]:
                diag -= 1.0
        H[idx, idx] += diag
        # Off-diagonal: S_q transverse field (all cyclic shifts)
        for site in range(n):
            for k in range(1, q):
                c2 = c.copy()
                c2[site] = (c[site] + k) % q
                idx2 = 0
                for i in range(n - 1, -1, -1):
                    idx2 = idx2 * q + c2[i]
                H[idx, idx2] += -g
    return csr_matrix(H)

def get_rho_A(psi, n, q, nA):
    """Get reduced density matrix for left nA sites."""
    dimA = q**nA
    dimB = q**(n - nA)
    psi_mat = psi.reshape(dimA, dimB)
    rho_A = psi_mat @ psi_mat.conj().T
    return np.real(rho_A)

def entanglement_hamiltonian(rho_A):
    """H_E = -log(rho_A), regularized."""
    eigvals, eigvecs = la.eigh(rho_A)
    eigvals = np.maximum(eigvals, 1e-15)
    H_E = -eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.T
    return np.real(H_E)

def build_potts_subsystem_ops(nA, q):
    """Build Potts-type operators on subsystem A (nA sites, q states each).
    Returns dict of operators: delta bonds, field terms."""
    dimA = q**nA
    ops = {}

    # delta(s_i, s_j) projector for bond between sites i and i+1
    delta_pair = np.zeros((q*q, q*q))
    for a in range(q):
        idx = a * q + a
        delta_pair[idx, idx] = 1.0

    for i in range(nA - 1):
        prefix = np.eye(q**i) if i > 0 else np.eye(1)
        suffix = np.eye(q**(nA - i - 2)) if i + 2 < nA else np.eye(1)
        ops[f'delta_{i}_{i+1}'] = np.kron(np.kron(prefix, delta_pair), suffix)

    # S_q field: sum_{k=1}^{q-1} X^k at each site
    X_single = np.zeros((q, q))
    for a in range(q):
        X_single[(a + 1) % q, a] = 1.0
    field_single = np.zeros((q, q))
    for k in range(1, q):
        field_single += np.linalg.matrix_power(X_single, k)

    for i in range(nA):
        prefix = np.eye(q**i) if i > 0 else np.eye(1)
        suffix = np.eye(q**(nA - i - 1)) if i < nA - 1 else np.eye(1)
        ops[f'field_{i}'] = np.kron(np.kron(prefix, field_single), suffix)

    return ops

def bw_envelope_periodic(nA, j, is_bond=False):
    """BW entanglement temperature envelope for periodic chain half-cut.
    For periodic chain of N=2*nA sites, half-chain A = sites 0..nA-1.
    Two cuts: between sites nA-1|nA and sites N-1|0.
    beta(x) = pi * sin(2*pi*x / N) where x is distance from left cut.
    """
    if is_bond:
        x = j + 1.0  # bond between j and j+1, centered at j+1
    else:
        x = j + 0.5  # site j centered at j+0.5
    return np.pi * np.sin(np.pi * x / nA)  # sin(2*pi*x/(2*nA)) = sin(pi*x/nA)

def build_bw_hamiltonian(nA, q, g, ops):
    """Build BW entanglement Hamiltonian with sin envelope."""
    dimA = q**nA
    H_BW = np.zeros((dimA, dimA))

    # Coupling terms with BW envelope
    for i in range(nA - 1):
        beta = bw_envelope_periodic(nA, i, is_bond=True)
        O = ops[f'delta_{i}_{i+1}']
        O_tl = O - np.trace(O) / dimA * np.eye(dimA)
        H_BW -= beta * O_tl

    # Field terms with BW envelope
    for i in range(nA):
        beta = bw_envelope_periodic(nA, i, is_bond=False)
        O = ops[f'field_{i}']
        O_tl = O - np.trace(O) / dimA * np.eye(dimA)
        H_BW -= g * beta * O_tl

    return H_BW

def compute_bw_metrics(H_E, H_BW):
    """Compute BW fidelity and variance captured."""
    dim = H_E.shape[0]
    H_E_tl = H_E - np.trace(H_E) / dim * np.eye(dim)
    H_BW_tl = H_BW - np.trace(H_BW) / dim * np.eye(dim)

    norm_E = la.norm(H_E_tl, 'fro')
    norm_BW = la.norm(H_BW_tl, 'fro')

    if norm_E < 1e-15 or norm_BW < 1e-15:
        return 0.0, 0.0, 0.0

    # Fidelity = normalized overlap
    fidelity = np.trace(H_E_tl.T @ H_BW_tl) / (norm_E * norm_BW)

    # Optimal rescaling alpha
    alpha = np.trace(H_E_tl.T @ H_BW_tl) / np.trace(H_BW_tl.T @ H_BW_tl)
    residual = H_E_tl - float(alpha) * H_BW_tl
    resid_frac = la.norm(residual, 'fro') / norm_E
    var_captured = 1 - resid_frac**2

    return float(np.real(fidelity)), float(np.real(var_captured)), float(np.real(alpha))

def compute_potts_locality(H_E, nA, q, ops):
    """Fraction of H_E captured by ALL Potts NN operators (without BW envelope).
    Uses least-squares fit: best linear combination of operators."""
    dimA = q**nA
    H_E_tl = H_E - np.trace(H_E) / dimA * np.eye(dimA)
    total_norm_sq = la.norm(H_E_tl, 'fro')**2
    if total_norm_sq < 1e-15:
        return 0.0, {}

    labels = list(ops.keys())
    op_vecs = []
    for label in labels:
        O = ops[label]
        O_tl = O - np.trace(O) / dimA * np.eye(dimA)
        op_vecs.append(O_tl.ravel())
    A_mat = np.array(op_vecs).T
    b_vec = H_E_tl.ravel()
    x_opt, _, _, _ = la.lstsq(A_mat, b_vec)
    fit = A_mat @ x_opt
    fit_norm_sq = np.dot(fit, fit)
    locality = float(np.real(fit_norm_sq / np.dot(b_vec, b_vec)))
    coeffs = {labels[i]: float(np.real(x_opt[i])) for i in range(len(labels))}
    return locality, coeffs

# ---- Main loop ----
configs = [
    (2, 12),  # dim=4096, dimA=64
    (3, 10),  # dim=59049, dimA=243
    (4, 8),   # dim=65536, dimA=256
    (5, 8),   # dim=390625, dimA=625
    (7, 6),   # dim=117649, dimA=343
]

# Known values
Re_c = {2: 0.500, 3: 0.800, 4: 1.000, 5: 1.138, 7: 1.351}
M_ratio = {2: 1.352, 3: 1.086, 4: 1.003, 5: 0.964, 7: 0.930}  # M/[(q-1)/q] at n=16
c_eff_ratio = {2: 1.00, 3: 1.12, 4: 1.00, 5: 1.01, 7: 0.78}

print("Sprint 091a: Entanglement Hamiltonian BW fidelity — S_q Potts")
print("=" * 70, flush=True)

for q, n in configs:
    gc = 1.0 / q
    nA = n // 2
    dimA = q**nA
    dim = q**n

    print(f"\nq={q}, n={n}, nA={nA}, dim={dim:,}, dimA={dimA:,}, g_c={gc:.4f}")
    print("-" * 60, flush=True)

    # Build Hamiltonian
    t1 = time.time()
    H = build_sq_potts_periodic(n, q, gc)
    t_build = time.time() - t1
    print(f"  H built: {t_build:.1f}s", flush=True)

    # Ground state
    t1 = time.time()
    evals, evecs = eigsh(H, k=2, which='SA')
    idx_sort = np.argsort(evals)
    psi = evecs[:, idx_sort[0]]
    gap = evals[idx_sort[1]] - evals[idx_sort[0]]
    t_eig = time.time() - t1
    print(f"  Eigsolve: {t_eig:.1f}s, E0={evals[idx_sort[0]]:.6f}, gap={gap:.6f}", flush=True)

    # Reduced density matrix
    rho_A = get_rho_A(psi, n, q, nA)
    rho_A = (rho_A + rho_A.T) / 2  # symmetrize

    # Entanglement entropy
    ev = la.eigvalsh(rho_A)
    ev_pos = ev[ev > 1e-15]
    S_vn = -np.sum(ev_pos * np.log(ev_pos))
    print(f"  S_vN = {S_vn:.6f} (nats)", flush=True)

    # Entanglement Hamiltonian
    H_E = entanglement_hamiltonian(rho_A)
    norm_HE = la.norm(H_E, 'fro')
    print(f"  ||H_E|| = {norm_HE:.4f}", flush=True)

    # Build subsystem operators
    ops = build_potts_subsystem_ops(nA, q)

    # Potts locality (best linear combination of NN operators)
    locality, coeffs = compute_potts_locality(H_E, nA, q, ops)
    print(f"  Potts NN locality: {locality*100:.2f}%", flush=True)

    # BW fidelity (with sin envelope)
    H_BW = build_bw_hamiltonian(nA, q, gc, ops)
    bw_fid, bw_var, bw_alpha = compute_bw_metrics(H_E, H_BW)
    print(f"  BW fidelity: {bw_fid:.4f}", flush=True)
    print(f"  BW variance captured: {bw_var*100:.2f}%", flush=True)
    print(f"  BW optimal alpha: {bw_alpha:.4f}", flush=True)

    # Also try uniform and linear envelopes for comparison
    dimA_loc = dimA
    env_results = {}
    for env_name, env_func in [('sin', lambda nA, j, ib: bw_envelope_periodic(nA, j, ib)),
                                 ('uniform', lambda nA, j, ib: 1.0),
                                 ('linear', lambda nA, j, ib: (j + 1.0 if ib else j + 0.5))]:
        # Build H_BW with this envelope
        H_test = np.zeros((dimA_loc, dimA_loc))
        for i in range(nA - 1):
            beta = env_func(nA, i, True)
            O = ops[f'delta_{i}_{i+1}']
            O_tl = O - np.trace(O) / dimA_loc * np.eye(dimA_loc)
            H_test -= beta * O_tl
        for i in range(nA):
            beta = env_func(nA, i, False)
            O = ops[f'field_{i}']
            O_tl = O - np.trace(O) / dimA_loc * np.eye(dimA_loc)
            H_test -= gc * beta * O_tl
        fid, var, alpha = compute_bw_metrics(H_E, H_test)
        env_results[env_name] = {'fidelity': float(fid), 'variance': float(var), 'alpha': float(alpha)}
    best_env = max(env_results, key=lambda e: env_results[e]['variance'])
    print(f"  Best envelope: {best_env} (var={env_results[best_env]['variance']*100:.2f}%)")

    # Coupling profile from fit
    print(f"  Fit coefficients:")
    for label in sorted(coeffs.keys()):
        if abs(coeffs[label]) > 1e-6:
            print(f"    {label}: {coeffs[label]:+.6f}")

    # BW envelope prediction for comparison
    print(f"  BW sin envelope prediction (scaled by alpha={bw_alpha:.3f}):")
    for i in range(nA - 1):
        beta = bw_envelope_periodic(nA, i, is_bond=True)
        print(f"    delta_{i}_{i+1}: fit={coeffs.get(f'delta_{i}_{i+1}', 0):+.6f}, "
              f"BW={-bw_alpha * beta:+.6f}")
    for i in range(nA):
        beta = bw_envelope_periodic(nA, i, is_bond=False)
        print(f"    field_{i}: fit={coeffs.get(f'field_{i}', 0):+.6f}, "
              f"BW={-gc * bw_alpha * beta:+.6f}")

    # Save data for this q
    results['data'][f'q{q}'] = {
        'q': q, 'n': n, 'nA': nA, 'dim': dim, 'dimA': dimA,
        'gc': gc, 'Re_c': Re_c.get(q, None),
        'S_vn': float(S_vn), 'gap': float(gap),
        'HE_norm': float(norm_HE),
        'potts_locality': float(locality),
        'bw_fidelity': float(bw_fid),
        'bw_variance': float(bw_var),
        'bw_alpha': float(bw_alpha),
        'envelope_comparison': env_results,
        'best_envelope': best_env,
        'fit_coefficients': coeffs,
        'time_s': float(time.time() - t1),
    }
    save()  # incremental save

# ---- Summary ----
print(f"\n{'='*70}")
print("SUMMARY: BW fidelity across walking boundary")
print(f"{'='*70}")
print(f"{'q':>3} {'n':>3} {'nA':>3} {'locality%':>10} {'BW_fid':>8} {'BW_var%':>8} "
      f"{'M/(q-1)/q':>10} {'c_eff/Rec':>10}")
for q, n in configs:
    d = results['data'][f'q{q}']
    mr = M_ratio.get(q, '?')
    cr = c_eff_ratio.get(q, '?')
    print(f"{q:>3} {n:>3} {d['nA']:>3} {d['potts_locality']*100:>10.2f} "
          f"{d['bw_fidelity']:>8.4f} {d['bw_variance']*100:>8.2f} "
          f"{mr:>10} {cr:>10}")

results['runtime_seconds'] = time.time() - t0_global
save()
print(f"\nTotal runtime: {time.time()-t0_global:.1f}s")
