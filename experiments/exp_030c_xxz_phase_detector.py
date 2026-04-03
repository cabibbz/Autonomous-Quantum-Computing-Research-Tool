"""
Sprint 030c: MI Uniformity CV as Universal Phase Detector
High-resolution sweep testing whether MI-CV and its derivatives detect:
- The first-order FM transition at Δ=-1
- The BKT transition at Δ=1 (notoriously hard to detect)
Compare XXZ to TFIM (Sprint 029), test finite-size scaling n=4,6,8
"""

import numpy as np
from scipy.sparse import kron, eye, csr_matrix
from scipy.sparse.linalg import eigsh
import json, time

# Pauli matrices
I2 = eye(2, format='csr')
Sp = csr_matrix(np.array([[0, 1], [0, 0]], dtype=complex))
Sm = csr_matrix(np.array([[0, 0], [1, 0]], dtype=complex))
Sz_op = csr_matrix(np.array([[0.5, 0], [0, -0.5]], dtype=complex))
X = csr_matrix(np.array([[0, 1], [1, 0]], dtype=complex))
Z = csr_matrix(np.array([[1, 0], [0, -1]], dtype=complex))

def chain_op(op, qubit, n):
    ops = [I2] * n
    ops[qubit] = op
    result = ops[0]
    for o in ops[1:]:
        result = kron(result, o, format='csr')
    return result

def xxz_hamiltonian(n, delta):
    dim = 2**n
    H = csr_matrix((dim, dim), dtype=complex)
    for i in range(n - 1):
        H += 0.5 * (chain_op(Sp, i, n) @ chain_op(Sm, i+1, n) +
                     chain_op(Sm, i, n) @ chain_op(Sp, i+1, n))
        H += delta * chain_op(Sz_op, i, n) @ chain_op(Sz_op, i+1, n)
    return H

def tfim_hamiltonian(n, h_over_J):
    dim = 2**n
    H = csr_matrix((dim, dim), dtype=complex)
    for i in range(n - 1):
        H -= chain_op(Z, i, n) @ chain_op(Z, i + 1, n)
    for i in range(n):
        H -= h_over_J * chain_op(X, i, n)
    return H

def project_to_sz_sector(H, n, sz_target):
    dim = 2**n
    basis_indices = []
    for i in range(dim):
        bits = format(i, f'0{n}b')
        sz = sum(0.5 if b == '0' else -0.5 for b in bits)
        if abs(sz - sz_target) < 1e-10:
            basis_indices.append(i)
    P = np.zeros((dim, len(basis_indices)))
    for j, idx in enumerate(basis_indices):
        P[idx, j] = 1.0
    H_dense = H.toarray()
    H_sector = P.T @ H_dense @ P
    return H_sector, P

def partial_trace_simple(state_vec, keep, n):
    keep = sorted(keep)
    trace_out = sorted(set(range(n)) - set(keep))
    rho = np.outer(state_vec, state_vec.conj()).reshape([2]*n + [2]*n)
    offset = 0
    for q in sorted(trace_out, reverse=True):
        rho = np.trace(rho, axis1=q, axis2=q + n - offset)
        offset += 1
    k = len(keep)
    return rho.reshape(2**k, 2**k)

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return -np.sum(evals * np.log2(evals))

def compute_mi_cv(psi, n):
    """Compute MI uniformity CV (coefficient of variation of all pairwise MI)."""
    mi_values = []
    for i in range(n):
        for j in range(i+1, n):
            S_i = von_neumann_entropy(partial_trace_simple(psi, [i], n))
            S_j = von_neumann_entropy(partial_trace_simple(psi, [j], n))
            S_ij = von_neumann_entropy(partial_trace_simple(psi, [i, j], n))
            mi = S_i + S_j - S_ij
            mi_values.append(mi)
    mi_arr = np.array(mi_values)
    if np.mean(mi_arr) > 1e-10:
        return float(np.std(mi_arr) / np.mean(mi_arr))
    return float('inf')

def compute_concurrence_pair(psi, i, j, n):
    """Concurrence between qubits i and j."""
    rho_ij = partial_trace_simple(psi, [i, j], n)
    # Spin-flip matrix
    sy = np.array([[0, -1j], [1j, 0]])
    sigma_yy = np.kron(sy, sy)
    rho_tilde = sigma_yy @ rho_ij.conj() @ sigma_yy
    product = rho_ij @ rho_tilde
    evals = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(product), 0))))[::-1]
    C = max(0, evals[0] - evals[1] - evals[2] - evals[3])
    return float(C)

def ground_state_xxz(n, delta):
    """Get XXZ ground state, using Sz=0 projection for FM phase."""
    H = xxz_hamiltonian(n, delta)
    if delta < -0.9:  # Use Sz=0 sector near/below FM transition
        H_sector, P = project_to_sz_sector(H, n, sz_target=0.0)
        evals, evecs = np.linalg.eigh(H_sector)
        psi = P @ evecs[:, 0]
        gap = evals[1] - evals[0] if len(evals) > 1 else 0
    else:
        evals, evecs = eigsh(H, k=2, which='SA')
        psi = evecs[:, 0]
        gap = evals[1] - evals[0]
    psi = psi / np.linalg.norm(psi)
    return psi, float(evals[0]), float(gap)

def ground_state_tfim(n, h):
    """Get TFIM ground state."""
    H = tfim_hamiltonian(n, h)
    evals, evecs = eigsh(H, k=2, which='SA')
    psi = evecs[:, 0]
    gap = evals[1] - evals[0]
    return psi / np.linalg.norm(psi), float(evals[0]), float(gap)

# ---- Part 1: High-resolution XXZ MI-CV sweep ----
print("=== Part 1: High-Resolution XXZ MI-CV Sweep ===\n")

results = {'xxz': {}, 'tfim': {}}

for n in [4, 6, 8]:
    print(f"--- n={n} ---")
    t0 = time.time()

    # Dense sweep near transitions
    delta_vals = np.unique(np.sort(np.concatenate([
        np.linspace(-2.0, -1.3, 8),
        np.linspace(-1.3, -0.7, 20),  # near FM transition
        np.linspace(-0.7, 0.5, 10),
        np.linspace(0.5, 1.5, 20),    # near BKT transition
        np.linspace(1.5, 3.0, 8),
    ])))

    n_results = {
        'delta': [], 'mi_cv': [], 'half_cut_S': [], 'gap': [],
        'avg_concurrence': [], 'energy': [],
    }

    for delta in delta_vals:
        if n == 8 and delta < -1.1:
            continue  # Skip deep FM for n=8 (Sz=0 sector is large)

        psi, E, gap = ground_state_xxz(n, delta)
        cv = compute_mi_cv(psi, n)
        S_half = von_neumann_entropy(partial_trace_simple(psi, list(range(n//2)), n))

        # Nearest-neighbor concurrence (fast)
        conc_vals = []
        for i in range(n-1):
            conc_vals.append(compute_concurrence_pair(psi, i, i+1, n))
        avg_conc = np.mean(conc_vals)

        n_results['delta'].append(float(delta))
        n_results['mi_cv'].append(float(cv) if np.isfinite(cv) else 999.0)
        n_results['half_cut_S'].append(float(S_half))
        n_results['gap'].append(float(gap))
        n_results['avg_concurrence'].append(float(avg_conc))
        n_results['energy'].append(float(E))

        if abs(delta - round(delta, 1)) < 0.02 or abs(abs(delta) - 1.0) < 0.05:
            print(f"  Δ={delta:+.3f}: CV={cv:.4f}, S={S_half:.4f}, gap={gap:.4f}, C={avg_conc:.4f}")

    elapsed = time.time() - t0
    print(f"  n={n} done in {elapsed:.1f}s")
    results['xxz'][f'n={n}'] = n_results

# ---- Part 2: TFIM comparison sweep ----
print("\n=== Part 2: TFIM Comparison Sweep ===\n")

for n in [4, 6, 8]:
    print(f"--- TFIM n={n} ---")
    t0 = time.time()

    h_vals = np.unique(np.sort(np.concatenate([
        np.linspace(0.05, 0.5, 5),
        np.linspace(0.5, 1.5, 20),  # near critical point
        np.linspace(1.5, 3.0, 5),
    ])))

    n_results = {'h': [], 'mi_cv': [], 'half_cut_S': [], 'gap': []}

    for h in h_vals:
        psi, E, gap = ground_state_tfim(n, h)
        cv = compute_mi_cv(psi, n)
        S_half = von_neumann_entropy(partial_trace_simple(psi, list(range(n//2)), n))

        n_results['h'].append(float(h))
        n_results['mi_cv'].append(float(cv) if np.isfinite(cv) else 999.0)
        n_results['half_cut_S'].append(float(S_half))
        n_results['gap'].append(float(gap))

        if abs(h - round(h, 1)) < 0.02 or abs(h - 1.0) < 0.05:
            print(f"  h={h:.3f}: CV={cv:.4f}, S={S_half:.4f}, gap={gap:.4f}")

    elapsed = time.time() - t0
    print(f"  n={n} done in {elapsed:.1f}s")
    results['tfim'][f'n={n}'] = n_results

# Save
with open('results/sprint_030c_xxz_phase_detector.json', 'w') as f:
    json.dump(results, f, indent=2)

# ---- Analysis ----
print("\n=== ANALYSIS: MI-CV as Phase Detector ===\n")

# XXZ: compute dCV/dΔ to look for transitions
for n_key in sorted(results['xxz'].keys()):
    data = results['xxz'][n_key]
    d = np.array(data['delta'])
    cv = np.array(data['mi_cv'])
    cv[cv > 100] = np.nan  # clean infinities

    # Numerical derivative
    dcv = np.gradient(cv, d)

    # Find peaks in |dCV/dΔ|
    valid = ~np.isnan(dcv)
    if np.any(valid):
        d_valid = d[valid]
        dcv_valid = dcv[valid]

        # Look near Δ=-1 and Δ=1
        near_fm = np.abs(d_valid - (-1.0)) < 0.3
        near_bkt = np.abs(d_valid - 1.0) < 0.3

        if np.any(near_fm):
            fm_peak = d_valid[near_fm][np.argmax(np.abs(dcv_valid[near_fm]))]
            fm_mag = np.max(np.abs(dcv_valid[near_fm]))
            print(f"{n_key} FM transition: peak |dCV/dΔ| = {fm_mag:.3f} at Δ = {fm_peak:.3f}")

        if np.any(near_bkt):
            bkt_peak = d_valid[near_bkt][np.argmax(np.abs(dcv_valid[near_bkt]))]
            bkt_mag = np.max(np.abs(dcv_valid[near_bkt]))
            print(f"{n_key} BKT transition: peak |dCV/dΔ| = {bkt_mag:.3f} at Δ = {bkt_peak:.3f}")

# TFIM: compute dCV/dh
print()
for n_key in sorted(results['tfim'].keys()):
    data = results['tfim'][n_key]
    h = np.array(data['h'])
    cv = np.array(data['mi_cv'])
    cv[cv > 100] = np.nan

    dcv = np.gradient(cv, h)
    valid = ~np.isnan(dcv)
    if np.any(valid):
        peak_idx = np.argmax(np.abs(dcv[valid]))
        print(f"TFIM {n_key}: peak |dCV/dh| = {np.abs(dcv[valid][peak_idx]):.3f} at h = {h[valid][peak_idx]:.3f}")

# Comparison: detector power
print("\n=== DETECTOR COMPARISON ===")
print("Metric           | FM (Δ=-1) | BKT (Δ=1) | TFIM (h=1)")
print("-" * 60)
for n_key in ['n=6']:
    xxz_data = results['xxz'][n_key]
    d = np.array(xxz_data['delta'])
    cv = np.array(xxz_data['mi_cv'])
    S = np.array(xxz_data['half_cut_S'])
    gap = np.array(xxz_data['gap'])
    conc = np.array(xxz_data['avg_concurrence'])

    # dS/dΔ
    dS = np.gradient(S, d)
    # dC/dΔ
    dC = np.gradient(conc, d)
    # dCV/dΔ
    cv_clean = cv.copy()
    cv_clean[cv_clean > 100] = np.nan
    dCV = np.gradient(cv_clean, d)

    for name, deriv in [('dS/dΔ', dS), ('dC/dΔ', dC), ('dCV/dΔ', dCV)]:
        near_fm = np.abs(d - (-1.0)) < 0.2
        near_bkt = np.abs(d - 1.0) < 0.2
        fm_val = np.nanmax(np.abs(deriv[near_fm])) if np.any(near_fm) else 0
        bkt_val = np.nanmax(np.abs(deriv[near_bkt])) if np.any(near_bkt) else 0
        print(f"{name:<17} | {fm_val:>9.3f} | {bkt_val:>9.3f} |")

print("\nResults saved to results/sprint_030c_xxz_phase_detector.json")
