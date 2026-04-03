"""
Sprint 034c: Chiral Clock Model — Genuine Z₃ BW Locality + Prediction Verification

34b showed the standard clock model is secretly S₃-symmetric (ZZ†+h.c. = 3δ-I).
The CHIRAL clock model breaks S₃ → Z₃:
  H = -J Σ (e^{iφ} Z_i Z_{i+1}† + e^{-iφ} Z_i† Z_{i+1}) - h Σ (X + X†)

When φ=0: S₃ symmetry (= Potts), peak BW = 76.5%
When φ≠0: only Z₃ symmetry (cyclic permutation X survives, but transposition breaks)

From 034a operator counting:
- Z₃ at d=3: 125 G-invariant operators (n_A=3)
- S₃ at d=3: 69 G-invariant operators (n_A=3)

PREDICTION: Chiral (Z₃-only) should have LOWER BW locality than φ=0 (S₃) = 76.5%

Also compile final prediction table comparing H/G-inv ratio to measured BW.
"""

import numpy as np
from scipy import linalg as la
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import eigsh
import json, time

t0 = time.time()

omega = np.exp(2j * np.pi / 3)

def clock_Z():
    return np.diag([1.0+0j, omega, omega**2])

def clock_X():
    return np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=complex)

def sparse_chiral_clock_hamiltonian(n, h_over_J, phi, boundary='open'):
    """
    H = -J Σ (e^{iφ} Z_i Z_{i+1}† + e^{-iφ} Z_i† Z_{i+1}) - h Σ (X + X†)
    """
    d = 3
    dim = d**n

    Z = clock_Z()
    Zd = Z.conj().T
    X = clock_X()
    Xd = X.conj().T

    # Chiral bond operator (Hermitian): e^{iφ} Z⊗Z† + e^{-iφ} Z†⊗Z
    ZZd_chiral = np.exp(1j*phi) * np.kron(Z, Zd) + np.exp(-1j*phi) * np.kron(Zd, Z)
    ZZd_sparse = csr_matrix(ZZd_chiral)

    XXd = X + Xd
    XXd_sparse = csr_matrix(XXd)

    H = csr_matrix((dim, dim), dtype=complex)

    n_bonds = n - 1 if boundary == 'open' else n
    for bond in range(n_bonds):
        j = (bond + 1) % n
        if j == bond + 1:
            prefix = sp_eye(d**bond, format='csr') if bond > 0 else sp_eye(1, format='csr')
            suffix = sp_eye(d**(n - bond - 2), format='csr') if bond + 2 < n else sp_eye(1, format='csr')
            H = H - sp_kron(sp_kron(prefix, ZZd_sparse), suffix, format='csr')

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
    H_E = -eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.conj().T  # V D V† for Hermitian
    return H_E

def chiral_subsystem_operators(n_A, phi):
    """
    Build chiral clock-type operators on subsystem A:
    - e^{iφ} Z_i Z_{i+1}† + e^{-iφ} Z_i† Z_{i+1} for adjacent bonds (Hermitian)
    - X_i + X_i† for each site (Hermitian)
    """
    d = 3
    dim = d**n_A
    Z = clock_Z()
    Zd = Z.conj().T
    X = clock_X()
    Xd = X.conj().T

    ZZd_chiral = np.exp(1j*phi) * np.kron(Z, Zd) + np.exp(-1j*phi) * np.kron(Zd, Z)
    XXd_single = X + Xd

    ops = {}

    for i in range(n_A - 1):
        prefix = np.eye(d**i, dtype=complex) if i > 0 else np.eye(1, dtype=complex)
        suffix = np.eye(d**(n_A - i - 2), dtype=complex) if i + 2 < n_A else np.eye(1, dtype=complex)
        op = np.kron(np.kron(prefix, ZZd_chiral), suffix)
        ops[f'ZZd_{i}_{i+1}'] = op

    for i in range(n_A):
        prefix = np.eye(d**i, dtype=complex) if i > 0 else np.eye(1, dtype=complex)
        suffix = np.eye(d**(n_A - i - 1), dtype=complex) if i < n_A - 1 else np.eye(1, dtype=complex)
        op = np.kron(np.kron(prefix, XXd_single), suffix)
        ops[f'XXd_{i}'] = op

    return ops

def compute_chiral_locality(H_E, n_A, phi):
    """Decompose H_E into chiral clock-type terms via least squares."""
    d = 3
    dim = d**n_A
    H_E_tl = H_E - np.trace(H_E) / dim * np.eye(dim)
    total_norm_sq = float(np.real(np.trace(H_E_tl.conj().T @ H_E_tl)))
    if total_norm_sq < 1e-15:
        return 0.0, {}, 0.0

    ops = chiral_subsystem_operators(n_A, phi)
    labels = list(ops.keys())

    # Need to work with Hermitian operators. The operators are Hermitian,
    # and we want real coefficients. Flatten as real vectors.
    op_vecs = []
    for label in labels:
        O = ops[label]
        O_tl = O - np.trace(O) / dim * np.eye(dim, dtype=complex)
        # For Hermitian operators, we can use the real part of the flattened vector
        # but need to be careful with complex entries
        op_vecs.append(O_tl.ravel())

    A_mat = np.array(op_vecs).T  # dim² × n_ops (complex)
    b_vec = H_E_tl.ravel()  # dim² (complex)

    # Complex least squares
    x_opt, _, _, _ = la.lstsq(A_mat, b_vec)

    fit = A_mat @ x_opt
    fit_norm_sq = float(np.real(np.vdot(fit, fit)))
    b_norm_sq = float(np.real(np.vdot(b_vec, b_vec)))
    locality = fit_norm_sq / b_norm_sq if b_norm_sq > 0 else 0

    coeffs = {labels[i]: float(np.real(x_opt[i])) for i in range(len(labels))}

    return float(locality), coeffs, float(b_norm_sq)

# ---- Sweep over chirality φ at fixed h/J near peak ----
n = 8
d = 3
dims = [d] * n
n_A = n // 2

# First: sweep φ at h/J = 0.75 (peak from 34b)
h_J_peak = 0.75
phi_values = np.linspace(0, np.pi/3, 13)  # 0 to π/3 (Z₃ period is 2π/3)

print(f"=== Sprint 034c: Chiral Clock BW Locality (n={n}, n_A={n_A}) ===")
print(f"\n--- Part 1: Chirality sweep at h/J={h_J_peak} ---")

phi_results = []

for phi in phi_values:
    t1 = time.time()

    H = sparse_chiral_clock_hamiltonian(n, h_J_peak, phi, 'open')
    H = (H + H.conj().T) / 2  # Ensure Hermitian

    eigvals_H, eigvecs_H = eigsh(H, k=4, which='SA')
    sort_idx = np.argsort(eigvals_H)
    psi_gs = eigvecs_H[:, sort_idx[0]]

    rho = np.outer(psi_gs, psi_gs.conj())
    rho_A = partial_trace_general(rho, dims, list(range(n_A)))
    rho_A = (rho_A + rho_A.conj().T) / 2

    ev = la.eigvalsh(rho_A)
    ev_sig = ev[ev > 1e-15]
    S = -np.sum(ev_sig * np.log2(ev_sig))

    H_E = entanglement_hamiltonian(rho_A)
    locality, coeffs, norm_sq = compute_chiral_locality(H_E, n_A, phi)

    dt = time.time() - t1

    result = {
        'phi': float(phi),
        'phi_over_pi': float(phi/np.pi),
        'h_J': h_J_peak,
        'entropy_bits': float(S),
        'locality': float(locality),
        'time_s': float(dt)
    }
    phi_results.append(result)

    sym = "S₃" if abs(phi) < 1e-10 else "Z₃"
    print(f"  φ/π={phi/np.pi:.4f} ({sym}): S={S:.3f}, locality={locality*100:.1f}% [{dt:.1f}s]")

# ---- Part 2: Full h/J sweep at selected chirality ----
print(f"\n--- Part 2: h/J sweep at φ=π/6 (genuine Z₃) ---")
phi_test = np.pi / 6
h_values = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 1.0, 1.5, 2.0]

h_results = []
for h_J in h_values:
    t1 = time.time()

    H = sparse_chiral_clock_hamiltonian(n, h_J, phi_test, 'open')
    H = (H + H.conj().T) / 2

    eigvals_H, eigvecs_H = eigsh(H, k=4, which='SA')
    sort_idx = np.argsort(eigvals_H)
    psi_gs = eigvecs_H[:, sort_idx[0]]

    rho = np.outer(psi_gs, psi_gs.conj())
    rho_A = partial_trace_general(rho, dims, list(range(n_A)))
    rho_A = (rho_A + rho_A.conj().T) / 2

    ev = la.eigvalsh(rho_A)
    ev_sig = ev[ev > 1e-15]
    S = -np.sum(ev_sig * np.log2(ev_sig))

    H_E = entanglement_hamiltonian(rho_A)
    locality, coeffs, norm_sq = compute_chiral_locality(H_E, n_A, phi_test)

    dt = time.time() - t1

    result = {
        'h_J': float(h_J),
        'phi': float(phi_test),
        'entropy_bits': float(S),
        'locality': float(locality),
        'time_s': float(dt)
    }
    h_results.append(result)

    print(f"  h/J={h_J:.3f}: S={S:.3f}, locality={locality*100:.1f}% [{dt:.1f}s]")

peak_chiral = max(h_results, key=lambda r: r['locality'])
peak_phi_sweep = max(phi_results, key=lambda r: r['locality'])
phi0_result = phi_results[0]

# ---- Part 3: Compile prediction table ----
print("\n=== FINAL PREDICTION TABLE ===")
print(f"{'Model':<16} {'Sym':<6} {'|G|':<5} {'d':<4} {'G-inv(n=3)':<12} {'H/G-inv':<10} {'BW Locality'}")
print("-" * 75)
print(f"{'XXZ':<16} {'U(1)':<6} {'∞':<5} {'2':<4} {'42(n=4)':<12} {'0.143':<10} {'100.0%'}")
print(f"{'TFIM':<16} {'Z₂':<6} {'2':<5} {'2':<4} {'71(n=4)':<12} {'0.099':<10} {'91.0%'}")
print(f"{'Potts/Clock':<16} {'S₃':<6} {'6':<5} {'3':<4} {'69':<12} {'0.072':<10} {'76.5%'}")
chiral_loc_pct = f"{peak_chiral['locality']*100:.1f}%"
print(f"{'Chiral Clock':<16} {'Z₃':<6} {'3':<5} {'3':<4} {'125':<12} {'0.040':<10} {chiral_loc_pct}")

# ---- Key comparison ----
print(f"\n=== KEY RESULTS ===")
print(f"φ=0 (S₃ Potts) peak BW: {phi0_result['locality']*100:.1f}%")
print(f"φ=π/6 (Z₃ chiral) peak BW: {peak_chiral['locality']*100:.1f}%")
delta = phi0_result['locality'] - peak_chiral['locality']
print(f"S₃ → Z₃ locality change: {delta*100:+.1f} percentage points")

if peak_chiral['locality'] < phi0_result['locality'] - 0.01:
    print("\n✓ CONFIRMED: Breaking S₃ → Z₃ REDUCES BW locality")
    print("  More symmetry → fewer invariant operators → more constrained H_E → better BW")
elif peak_chiral['locality'] > phi0_result['locality'] + 0.01:
    print("\n✗ OVERTURNED: Z₃ has BETTER BW despite less symmetry")
else:
    print("\n~ INCONCLUSIVE: Z₃ and S₃ BW within 1% — chirality doesn't affect BW strongly")

# ---- Chirality trend ----
print(f"\n=== CHIRALITY TREND (h/J={h_J_peak}) ===")
print(f"{'φ/π':<10} {'Locality%':<12} {'ΔS₃'}")
for r in phi_results:
    delta_loc = (r['locality'] - phi0_result['locality']) * 100
    print(f"{r['phi_over_pi']:<10.4f} {r['locality']*100:<12.1f} {delta_loc:+.1f}%")

# ---- Save ----
results = {
    'experiment': '034c',
    'description': 'Chiral clock model BW locality and prediction verification',
    'parameters': {'n': n, 'local_dim': d, 'n_A': n_A},
    'chirality_sweep': phi_results,
    'h_sweep_chiral': h_results,
    'peak_chiral': {
        'locality': peak_chiral['locality'],
        'h_J': peak_chiral['h_J'],
        'phi': float(phi_test)
    },
    'phi0_locality': phi0_result['locality'],
    'prediction_table': {
        'XXZ_U1': {'group_order': 'inf', 'd': 2, 'bw_locality': 1.0},
        'TFIM_Z2': {'group_order': 2, 'd': 2, 'bw_locality': 0.91},
        'Potts_S3': {'group_order': 6, 'd': 3, 'bw_locality': 0.765},
        'Chiral_Z3': {'group_order': 3, 'd': 3, 'bw_locality': peak_chiral['locality']},
    },
    'runtime_seconds': time.time() - t0
}

with open('results/sprint_034c_chiral_bw.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nSaved. Runtime: {time.time()-t0:.1f}s")
