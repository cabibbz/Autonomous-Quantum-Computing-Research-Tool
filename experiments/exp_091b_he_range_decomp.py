#!/usr/bin/env python3
"""Sprint 091b: H_E operator range decomposition at FIXED nA=4 for q=2,3,4,5.

091a showed BW fidelity varies with nA (confound). Fix nA=4 for all q to get
a clean comparison. Also add NNN delta operators to measure how much non-local
content exists at each q.

Additionally: nA scaling for q=2 (nA=3,4,5,6) to see size dependence.
"""
import numpy as np
from scipy import linalg as la
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
import json, time

t0_global = time.time()

results = {
    'experiment': '091b_he_range_decomp',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_091b_he_range_decomp.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_periodic(n, q, g):
    dim = q**n
    H = lil_matrix((dim, dim), dtype=complex)
    for idx in range(dim):
        c = []
        tmp = idx
        for _ in range(n):
            c.append(tmp % q)
            tmp //= q
        diag = 0.0
        for site in range(n):
            nxt = (site + 1) % n
            if c[site] == c[nxt]:
                diag -= 1.0
        H[idx, idx] += diag
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
    dimA = q**nA
    dimB = q**(n - nA)
    psi_mat = psi.reshape(dimA, dimB)
    rho_A = psi_mat @ psi_mat.conj().T
    return np.real(rho_A)

def entanglement_hamiltonian(rho_A):
    eigvals, eigvecs = la.eigh(rho_A)
    eigvals = np.maximum(eigvals, 1e-15)
    H_E = -eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.T
    return np.real(H_E)

def build_operators(nA, q, max_range=3):
    """Build Potts operators of different ranges on subsystem A.
    Range 1: single-site field operators
    Range 2 (NN): nearest-neighbor delta operators
    Range 2 (NNN): next-nearest-neighbor delta operators
    Range 3: 3-body delta chain delta(s_i,s_{i+1})*delta(s_{i+1},s_{i+2})
    """
    dimA = q**nA
    ops = {'range1': {}, 'range2_nn': {}, 'range2_nnn': {}, 'range3': {}}

    # delta(s_i, s_j) projector
    delta_pair = np.zeros((q*q, q*q))
    for a in range(q):
        idx = a * q + a
        delta_pair[idx, idx] = 1.0

    # Field: sum_{k=1}^{q-1} X^k
    X_single = np.zeros((q, q))
    for a in range(q):
        X_single[(a + 1) % q, a] = 1.0
    field_single = np.zeros((q, q))
    for k in range(1, q):
        field_single += np.linalg.matrix_power(X_single, k)

    # Range 1: field terms
    for i in range(nA):
        prefix = np.eye(q**i) if i > 0 else np.eye(1)
        suffix = np.eye(q**(nA - i - 1)) if i < nA - 1 else np.eye(1)
        ops['range1'][f'field_{i}'] = np.kron(np.kron(prefix, field_single), suffix)

    # Range 2 NN: delta nearest-neighbor
    for i in range(nA - 1):
        prefix = np.eye(q**i) if i > 0 else np.eye(1)
        suffix = np.eye(q**(nA - i - 2)) if i + 2 < nA else np.eye(1)
        ops['range2_nn'][f'delta_{i}_{i+1}'] = np.kron(np.kron(prefix, delta_pair), suffix)

    # Range 2 NNN: delta next-nearest-neighbor
    if max_range >= 2 and nA >= 3:
        # delta(s_i, s_{i+2}): need to build as q^3 operator then embed
        delta_nnn = np.zeros((q**3, q**3))
        for a in range(q):
            for b in range(q):  # middle site free
                row = a * q*q + b * q + a  # s0=a, s1=b, s2=a
                delta_nnn[row, row] = 1.0  # actually we want delta(s0, s2)
        # Correct: delta(s_i, s_{i+2}) = sum_a |a><a|_i tensor I_{i+1} tensor |a><a|_{i+2}
        delta_nnn = np.zeros((q**3, q**3))
        for a in range(q):
            for b in range(q):
                row = a * q**2 + b * q + a
                col = row
                delta_nnn[row, col] = 1.0

        for i in range(nA - 2):
            prefix = np.eye(q**i) if i > 0 else np.eye(1)
            suffix = np.eye(q**(nA - i - 3)) if i + 3 < nA else np.eye(1)
            ops['range2_nnn'][f'delta_{i}_{i+2}'] = np.kron(np.kron(prefix, delta_nnn), suffix)

    # Range 3: 3-body delta chain = delta(s_i,s_{i+1}) * delta(s_{i+1},s_{i+2})
    if max_range >= 3 and nA >= 3:
        # = sum_a |a,a,a><a,a,a| : all three sites equal
        delta_3body = np.zeros((q**3, q**3))
        for a in range(q):
            idx = a * q**2 + a * q + a
            delta_3body[idx, idx] = 1.0

        for i in range(nA - 2):
            prefix = np.eye(q**i) if i > 0 else np.eye(1)
            suffix = np.eye(q**(nA - i - 3)) if i + 3 < nA else np.eye(1)
            ops['range3'][f'delta3_{i}_{i+1}_{i+2}'] = np.kron(np.kron(prefix, delta_3body), suffix)

    return ops

def project_onto_operators(H_E, ops_dict, dimA):
    """Project H_E onto a set of operators via least-squares. Return variance captured."""
    H_E_tl = H_E - np.trace(H_E) / dimA * np.eye(dimA)
    total_norm_sq = la.norm(H_E_tl, 'fro')**2
    if total_norm_sq < 1e-15:
        return 0.0, {}

    # Flatten all operators from all groups
    labels = []
    for group_name in sorted(ops_dict.keys()):
        for label in sorted(ops_dict[group_name].keys()):
            labels.append((group_name, label))

    op_vecs = []
    for group_name, label in labels:
        O = ops_dict[group_name][label]
        O_tl = O - np.trace(O) / dimA * np.eye(dimA)
        op_vecs.append(O_tl.ravel())

    A_mat = np.array(op_vecs).T
    b_vec = H_E_tl.ravel()
    x_opt, _, _, _ = la.lstsq(A_mat, b_vec)
    fit = A_mat @ x_opt
    fit_norm_sq = np.dot(fit, fit)
    total_var = float(np.real(fit_norm_sq / np.dot(b_vec, b_vec)))

    coeffs = {f"{gn}/{lb}": float(np.real(x_opt[i])) for i, (gn, lb) in enumerate(labels)}
    return total_var, coeffs

def bw_envelope_periodic(nA, j, is_bond=False):
    if is_bond:
        x = j + 1.0
    else:
        x = j + 0.5
    return np.pi * np.sin(np.pi * x / nA)

def bw_metrics(H_E, nA, q, gc, ops):
    """Compute BW fidelity using sin envelope on NN operators."""
    dimA = q**nA
    H_E_tl = H_E - np.trace(H_E) / dimA * np.eye(dimA)
    norm_E = la.norm(H_E_tl, 'fro')
    if norm_E < 1e-15:
        return 0.0, 0.0

    H_BW = np.zeros((dimA, dimA))
    for label, O in ops['range2_nn'].items():
        i = int(label.split('_')[1])
        beta = bw_envelope_periodic(nA, i, is_bond=True)
        O_tl = O - np.trace(O) / dimA * np.eye(dimA)
        H_BW -= beta * O_tl
    for label, O in ops['range1'].items():
        i = int(label.split('_')[1])
        beta = bw_envelope_periodic(nA, i, is_bond=False)
        O_tl = O - np.trace(O) / dimA * np.eye(dimA)
        H_BW -= gc * beta * O_tl

    H_BW_tl = H_BW - np.trace(H_BW) / dimA * np.eye(dimA)
    norm_BW = la.norm(H_BW_tl, 'fro')
    if norm_BW < 1e-15:
        return 0.0, 0.0

    fid = float(np.real(np.trace(H_E_tl.T @ H_BW_tl) / (norm_E * norm_BW)))
    alpha = float(np.real(np.trace(H_E_tl.T @ H_BW_tl) / np.trace(H_BW_tl.T @ H_BW_tl)))
    resid = H_E_tl - alpha * H_BW_tl
    var = 1 - (la.norm(resid, 'fro') / norm_E)**2
    return fid, float(var)


# ---- Part 1: Fixed nA=4, varying q ----
print("=" * 70)
print("PART 1: Fixed nA=4, varying q (controlled comparison)")
print("=" * 70, flush=True)

fixed_nA_configs = [
    (2, 8),   # dim=256, dimA=16
    (3, 8),   # dim=6561, dimA=81
    (4, 8),   # dim=65536, dimA=256
    (5, 8),   # dim=390625, dimA=625
]

for q, n in fixed_nA_configs:
    gc = 1.0 / q
    nA = n // 2  # = 4 for all
    dimA = q**nA
    dim = q**n
    print(f"\nq={q}, n={n}, nA={nA}, dim={dim:,}, dimA={dimA:,}")
    print("-" * 50, flush=True)

    t1 = time.time()
    H = build_sq_potts_periodic(n, q, gc)
    print(f"  H built: {time.time()-t1:.1f}s", flush=True)

    t2 = time.time()
    evals, evecs = eigsh(H, k=2, which='SA')
    idx_sort = np.argsort(evals)
    psi = evecs[:, idx_sort[0]]
    print(f"  Eigsolve: {time.time()-t2:.1f}s", flush=True)

    rho_A = get_rho_A(psi, n, q, nA)
    rho_A = (rho_A + rho_A.T) / 2
    H_E = entanglement_hamiltonian(rho_A)

    # Build operators at all ranges
    ops = build_operators(nA, q, max_range=3)
    n_ops = {k: len(v) for k, v in ops.items()}
    print(f"  Operators: {n_ops}", flush=True)

    # Progressive projection: NN only, then +NNN, then +3body
    # NN only (range1 + range2_nn)
    ops_nn = {'range1': ops['range1'], 'range2_nn': ops['range2_nn']}
    var_nn, _ = project_onto_operators(H_E, ops_nn, dimA)

    # NN + NNN
    ops_nnn = {**ops_nn, 'range2_nnn': ops['range2_nnn']}
    var_nnn, _ = project_onto_operators(H_E, ops_nnn, dimA)

    # NN + NNN + 3body
    var_all, coeffs_all = project_onto_operators(H_E, ops, dimA)

    # BW fidelity
    bw_fid, bw_var = bw_metrics(H_E, nA, q, gc, ops)

    print(f"  NN locality:       {var_nn*100:.3f}%")
    print(f"  +NNN locality:     {var_nnn*100:.3f}%")
    print(f"  +3body locality:   {var_all*100:.3f}%")
    print(f"  BW fidelity:       {bw_fid:.4f}")
    print(f"  BW variance:       {bw_var*100:.3f}%")
    print(f"  NNN contribution:  {(var_nnn-var_nn)*100:.3f}%")
    print(f"  3body contribution:{(var_all-var_nnn)*100:.3f}%")
    print(f"  Residual (non-local): {(1-var_all)*100:.3f}%")

    results['data'][f'fixed_nA4_q{q}'] = {
        'q': q, 'n': n, 'nA': nA, 'dimA': dimA,
        'var_nn': float(var_nn),
        'var_nnn': float(var_nnn),
        'var_all': float(var_all),
        'bw_fidelity': float(bw_fid),
        'bw_variance': float(bw_var),
        'nnn_contribution': float(var_nnn - var_nn),
        'threebody_contribution': float(var_all - var_nnn),
        'residual': float(1 - var_all),
        'coefficients': coeffs_all,
        'time_s': float(time.time() - t1),
    }
    save()

# ---- Part 2: nA scaling for q=2 ----
print(f"\n{'='*70}")
print("PART 2: nA scaling for q=2 (BW fidelity vs subsystem size)")
print("=" * 70, flush=True)

q2_configs = [
    (2, 6),   # nA=3, dimA=8
    (2, 8),   # nA=4, dimA=16
    (2, 10),  # nA=5, dimA=32
    (2, 12),  # nA=6, dimA=64
    (2, 14),  # nA=7, dimA=128
]

for q, n in q2_configs:
    gc = 1.0 / q
    nA = n // 2
    dimA = q**nA
    dim = q**n
    print(f"\n  q={q}, n={n}, nA={nA}, dim={dim:,}, dimA={dimA:,} ... ", end='', flush=True)

    t1 = time.time()
    H = build_sq_potts_periodic(n, q, gc)
    evals, evecs = eigsh(H, k=2, which='SA')
    idx_sort = np.argsort(evals)
    psi = evecs[:, idx_sort[0]]
    rho_A = get_rho_A(psi, n, q, nA)
    rho_A = (rho_A + rho_A.T) / 2
    H_E = entanglement_hamiltonian(rho_A)

    ops = build_operators(nA, q, max_range=3)
    ops_nn = {'range1': ops['range1'], 'range2_nn': ops['range2_nn']}
    var_nn, _ = project_onto_operators(H_E, ops_nn, dimA)
    ops_nnn = {**ops_nn, 'range2_nnn': ops['range2_nnn']}
    var_nnn, _ = project_onto_operators(H_E, ops_nnn, dimA)
    var_all, _ = project_onto_operators(H_E, ops, dimA)
    bw_fid, bw_var = bw_metrics(H_E, nA, q, gc, ops)

    dt = time.time() - t1
    print(f"NN={var_nn*100:.2f}%, +NNN={var_nnn*100:.2f}%, BW_fid={bw_fid:.4f} [{dt:.1f}s]")

    results['data'][f'q2_nA{nA}'] = {
        'q': q, 'n': n, 'nA': nA, 'dimA': dimA,
        'var_nn': float(var_nn),
        'var_nnn': float(var_nnn),
        'var_all': float(var_all),
        'bw_fidelity': float(bw_fid),
        'bw_variance': float(bw_var),
        'time_s': float(dt),
    }
    save()

# ---- Part 3: nA scaling for q=5 ----
print(f"\n{'='*70}")
print("PART 3: nA scaling for q=5 (walking regime)")
print("=" * 70, flush=True)

q5_configs = [
    (5, 6),   # nA=3, dimA=125
    (5, 8),   # nA=4, dimA=625
]

for q, n in q5_configs:
    gc = 1.0 / q
    nA = n // 2
    dimA = q**nA
    dim = q**n
    print(f"\n  q={q}, n={n}, nA={nA}, dim={dim:,}, dimA={dimA:,} ... ", end='', flush=True)

    t1 = time.time()
    H = build_sq_potts_periodic(n, q, gc)
    evals, evecs = eigsh(H, k=2, which='SA')
    idx_sort = np.argsort(evals)
    psi = evecs[:, idx_sort[0]]
    rho_A = get_rho_A(psi, n, q, nA)
    rho_A = (rho_A + rho_A.T) / 2
    H_E = entanglement_hamiltonian(rho_A)

    ops = build_operators(nA, q, max_range=3)
    ops_nn = {'range1': ops['range1'], 'range2_nn': ops['range2_nn']}
    var_nn, _ = project_onto_operators(H_E, ops_nn, dimA)
    ops_nnn = {**ops_nn, 'range2_nnn': ops['range2_nnn']}
    var_nnn, _ = project_onto_operators(H_E, ops_nnn, dimA)
    var_all, _ = project_onto_operators(H_E, ops, dimA)
    bw_fid, bw_var = bw_metrics(H_E, nA, q, gc, ops)

    dt = time.time() - t1
    print(f"NN={var_nn*100:.2f}%, +NNN={var_nnn*100:.2f}%, BW_fid={bw_fid:.4f} [{dt:.1f}s]")

    results['data'][f'q5_nA{nA}'] = {
        'q': q, 'n': n, 'nA': nA, 'dimA': dimA,
        'var_nn': float(var_nn),
        'var_nnn': float(var_nnn),
        'var_all': float(var_all),
        'bw_fidelity': float(bw_fid),
        'bw_variance': float(bw_var),
        'time_s': float(dt),
    }
    save()

# ---- Summary ----
print(f"\n{'='*70}")
print("SUMMARY: Fixed nA=4")
print(f"{'='*70}")
print(f"{'q':>3} {'NN%':>8} {'+NNN%':>8} {'+3body%':>8} {'BW_fid':>8} {'BW_var%':>8} {'resid%':>8}")
for q in [2, 3, 4, 5]:
    d = results['data'][f'fixed_nA4_q{q}']
    print(f"{q:>3} {d['var_nn']*100:>8.3f} {d['var_nnn']*100:>8.3f} {d['var_all']*100:>8.3f} "
          f"{d['bw_fidelity']:>8.4f} {d['bw_variance']*100:>8.3f} {d['residual']*100:>8.3f}")

print(f"\nSUMMARY: q=2 nA scaling")
print(f"{'nA':>3} {'NN%':>8} {'+NNN%':>8} {'BW_fid':>8}")
for nA in [3, 4, 5, 6, 7]:
    d = results['data'].get(f'q2_nA{nA}')
    if d:
        print(f"{nA:>3} {d['var_nn']*100:>8.3f} {d['var_nnn']*100:>8.3f} {d['bw_fidelity']:>8.4f}")

print(f"\nSUMMARY: q=5 nA scaling")
print(f"{'nA':>3} {'NN%':>8} {'+NNN%':>8} {'BW_fid':>8}")
for nA in [3, 4]:
    d = results['data'].get(f'q5_nA{nA}')
    if d:
        print(f"{nA:>3} {d['var_nn']*100:>8.3f} {d['var_nnn']*100:>8.3f} {d['bw_fidelity']:>8.4f}")

results['runtime_seconds'] = time.time() - t0_global
save()
print(f"\nTotal runtime: {time.time()-t0_global:.1f}s")
