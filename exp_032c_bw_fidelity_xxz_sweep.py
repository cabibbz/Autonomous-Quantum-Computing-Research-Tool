"""
Sprint 032c: BW Fidelity for XXZ Model Across Phases

Test BW locality and envelope across the XXZ phase diagram:
- FM phase (Δ < -1): first-order transition
- XY phase (-1 < Δ < 1): gapless critical
- Néel phase (Δ > 1): BKT transition

Compare to TFIM results. Does BW locality depend on symmetry (Z2 vs U(1))?
Does gaplessness enhance or suppress BW accuracy?

System: n=8, subsystem A = left 4 qubits.
"""

import numpy as np
from scipy import linalg as la
import json, time

t0 = time.time()

# ---- Reusable functions ----
def pauli_matrices():
    I = np.eye(2)
    X = np.array([[0,1],[1,0]])
    Y = np.array([[0,-1j],[1j,0]])
    Z = np.array([[1,0],[0,-1]])
    return I, X, Y, Z

def kron_list(ops):
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result

def xxz_hamiltonian(n, delta):
    """H = Σ (X_i X_{i+1} + Y_i Y_{i+1} + Δ Z_i Z_{i+1}) [open boundary]"""
    I, X, Y, Z = pauli_matrices()
    dim = 2**n
    H = np.zeros((dim, dim), dtype=complex)
    for i in range(n-1):
        for P in [X, Y]:
            ops = [I]*n; ops[i] = P; ops[i+1] = P
            H += kron_list(ops)
        ops = [I]*n; ops[i] = Z; ops[i+1] = Z
        H += delta * kron_list(ops)
    return H

def partial_trace(rho, n_total, keep_qubits):
    n_keep = len(keep_qubits)
    trace_qubits = [q for q in range(n_total) if q not in keep_qubits]
    rho_tensor = rho.reshape([2]*n_total + [2]*n_total)
    for q in sorted(trace_qubits, reverse=True):
        n_remaining = rho_tensor.ndim // 2
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q+n_remaining)
    dim_keep = 2**n_keep
    return rho_tensor.reshape(dim_keep, dim_keep)

def entanglement_hamiltonian(rho_A):
    eigvals, eigvecs = la.eigh(rho_A)
    eigvals = np.maximum(eigvals, 1e-15)
    H_E = -eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.T
    return H_E

def matrix_overlap(A, B):
    nA = la.norm(A, 'fro')
    nB = la.norm(B, 'fro')
    if nA < 1e-15 or nB < 1e-15:
        return 0.0
    return np.real(np.trace(A.conj().T @ B)) / (nA * nB)

def bw_hamiltonian_xxz(n_A, delta, envelope_type='sin_inv'):
    """BW prediction for XXZ subsystem."""
    I, X, Y, Z = pauli_matrices()
    dim = 2**n_A
    L = 2 * n_A

    def get_beta_bond(i):
        d = n_A - 1 - i
        if envelope_type == 'linear':
            return d + 0.5
        elif envelope_type == 'sin_inv':
            return np.sin(np.pi * (n_A - i - 1) / L)
        elif envelope_type == 'uniform':
            return 1.0

    H_BW = np.zeros((dim, dim), dtype=complex)
    for i in range(n_A - 1):
        beta = get_beta_bond(i)
        for P in [X, Y]:
            ops = [I]*n_A; ops[i] = P; ops[i+1] = P
            H_BW += beta * kron_list(ops)
        ops = [I]*n_A; ops[i] = Z; ops[i+1] = Z
        H_BW += delta * beta * kron_list(ops)
    return H_BW

def compute_bw_metrics(H_E, H_BW):
    H_E_tl = H_E - np.trace(H_E)/H_E.shape[0] * np.eye(H_E.shape[0])
    H_BW_tl = H_BW - np.trace(H_BW)/H_BW.shape[0] * np.eye(H_BW.shape[0])

    # Handle complex matrices
    H_E_tl = np.real(H_E_tl) if np.max(np.abs(np.imag(H_E_tl))) < 1e-10 else H_E_tl
    H_BW_tl = np.real(H_BW_tl) if np.max(np.abs(np.imag(H_BW_tl))) < 1e-10 else H_BW_tl

    fidelity = matrix_overlap(H_E_tl, H_BW_tl)

    denom = np.real(np.trace(H_BW_tl.conj().T @ H_BW_tl))
    if denom < 1e-15:
        return 0.0, 0.0, 0.0
    alpha = np.real(np.trace(H_E_tl.conj().T @ H_BW_tl)) / denom
    residual = H_E_tl - alpha * H_BW_tl
    resid_frac = la.norm(residual, 'fro') / la.norm(H_E_tl, 'fro')
    var_captured = 1 - resid_frac**2

    return float(fidelity), float(var_captured), float(alpha)

def compute_locality_xxz(H_E, n_A, delta):
    """Fraction of H_E in XXZ-type terms."""
    I, X, Y, Z = pauli_matrices()
    dim = 2**n_A
    H_E_r = np.real(H_E) if np.max(np.abs(np.imag(H_E))) < 1e-10 else H_E
    H_E_tl = H_E_r - np.trace(H_E_r)/dim * np.eye(dim)
    total_norm_sq = la.norm(H_E_tl, 'fro')**2

    if total_norm_sq < 1e-15:
        return 0.0

    xxz_norm_sq = 0
    for i in range(n_A - 1):
        for P in [X, Y, Z]:
            ops = [I]*n_A; ops[i] = P; ops[i+1] = P
            op = kron_list(ops)
            coeff = np.real(np.trace(H_E_tl @ op)) / dim
            xxz_norm_sq += coeff**2

    return float(xxz_norm_sq * dim / total_norm_sq)

# ---- Main sweep ----
n = 8
n_A = 4

# Dense sampling around transitions
delta_values = np.concatenate([
    np.linspace(-2.0, -1.2, 5),   # deep FM
    np.linspace(-1.1, -0.9, 5),   # FM transition
    np.linspace(-0.5, 0.5, 5),    # XY phase
    np.linspace(0.7, 1.3, 7),     # BKT region
    np.linspace(1.5, 2.5, 4),     # Néel phase
])

envelope_types = ['linear', 'sin_inv', 'uniform']

print(f"=== Sprint 032c: BW Fidelity Across XXZ Phase Diagram ===")
print(f"n={n}, n_A={n_A}, {len(delta_values)} Δ values")

results_list = []

for delta in delta_values:
    H = xxz_hamiltonian(n, delta)
    eigvals, eigvecs = la.eigh(H)
    psi_gs = eigvecs[:, 0]
    gap = float(eigvals[1] - eigvals[0])

    rho = np.outer(psi_gs, psi_gs.conj())
    rho_A = partial_trace(rho, n, list(range(n_A)))

    # Make Hermitian (numerical cleanup)
    rho_A = (rho_A + rho_A.conj().T) / 2
    rho_A = np.real(rho_A) if np.max(np.abs(np.imag(rho_A))) < 1e-10 else rho_A

    H_E = entanglement_hamiltonian(rho_A)

    ev = la.eigvalsh(rho_A)
    ev = ev[ev > 1e-15]
    S = float(-np.sum(ev * np.log2(ev)))

    loc_frac = compute_locality_xxz(H_E, n_A, delta)

    env_results = {}
    for etype in envelope_types:
        H_BW = bw_hamiltonian_xxz(n_A, delta, etype)
        # Make real if possible
        H_BW = np.real(H_BW) if np.max(np.abs(np.imag(H_BW))) < 1e-10 else H_BW
        fid, var_cap, alpha = compute_bw_metrics(H_E, H_BW)
        env_results[etype] = {'fidelity': fid, 'variance_captured': var_cap, 'alpha': alpha}

    best_env = max(envelope_types, key=lambda e: env_results[e]['variance_captured'])

    # Extract effective coupling ratios for bonds
    I2, X2, Y2, Z2 = pauli_matrices()
    bond_coeffs = []
    for i in range(n_A - 1):
        xx_coeff = yy_coeff = zz_coeff = 0.0
        for P, label in [(X2, 'XX'), (Y2, 'YY'), (Z2, 'ZZ')]:
            ops = [I2]*n_A; ops[i] = P; ops[i+1] = P
            op = kron_list(ops)
            c = float(np.real(np.trace(H_E @ op)) / (2**n_A))
            if label == 'XX': xx_coeff = c
            elif label == 'YY': yy_coeff = c
            elif label == 'ZZ': zz_coeff = c
        bond_coeffs.append({'bond': f'{i}-{i+1}', 'XX': xx_coeff, 'YY': yy_coeff, 'ZZ': zz_coeff})

    result = {
        'delta': float(delta),
        'energy_gap': gap,
        'entropy_bits': S,
        'locality_fraction': loc_frac,
        'envelope_results': env_results,
        'best_envelope': best_env,
        'best_variance_captured': env_results[best_env]['variance_captured'],
        'bond_coefficients': bond_coeffs
    }
    results_list.append(result)

    print(f"Δ={delta:+.2f}: gap={gap:.3f}, S={S:.3f}, locality={loc_frac:.3f}, "
          f"best={best_env}({env_results[best_env]['variance_captured']:.3f})")

# ---- Analysis ----
print("\n=== Analysis ===")

# Phase classification
fm = [r for r in results_list if r['delta'] < -1.05]
fm_trans = [r for r in results_list if -1.05 <= r['delta'] <= -0.95]
xy = [r for r in results_list if -0.95 < r['delta'] < 0.95]
bkt = [r for r in results_list if 0.95 <= r['delta'] <= 1.1]
neel = [r for r in results_list if r['delta'] > 1.1]

for phase_name, phase_data in [('FM', fm), ('FM transition', fm_trans),
                                ('XY', xy), ('BKT', bkt), ('Néel', neel)]:
    if not phase_data:
        continue
    avg_loc = np.mean([r['locality_fraction'] for r in phase_data])
    avg_bw = np.mean([r['best_variance_captured'] for r in phase_data])
    print(f"\n{phase_name} (n={len(phase_data)}):")
    print(f"  Avg locality: {avg_loc:.4f}")
    print(f"  Avg best BW variance: {avg_bw:.4f}")

# Compare to TFIM
print("\n--- XXZ vs TFIM BW comparison ---")
max_loc_xxz = max(results_list, key=lambda r: r['locality_fraction'])
min_loc_xxz = min(results_list, key=lambda r: r['locality_fraction'])
print(f"XXZ max locality: {max_loc_xxz['locality_fraction']:.4f} at Δ={max_loc_xxz['delta']:.2f}")
print(f"XXZ min locality: {min_loc_xxz['locality_fraction']:.4f} at Δ={min_loc_xxz['delta']:.2f}")
print(f"TFIM max (from 32b): ~0.921 at h/J=0.80")
print(f"TFIM critical: ~0.909 at h/J=1.00")

# XX=YY check (U(1) symmetry)
print("\n--- U(1) symmetry check: XX vs YY coefficients ---")
for r in results_list[::5]:  # sample every 5th
    bc = r['bond_coefficients'][0]
    xy_ratio = bc['XX'] / bc['YY'] if abs(bc['YY']) > 1e-10 else float('inf')
    print(f"  Δ={r['delta']:+.2f}: XX/YY ratio at bond 0-1 = {xy_ratio:.4f}")

# Position dependence
print("\n--- Bond coefficient ratios (testing BW envelope shape) ---")
for r in results_list[::5]:
    bcs = r['bond_coefficients']
    if abs(bcs[2]['XX']) > 1e-10:
        ratio_01_23 = bcs[0]['XX'] / bcs[2]['XX']
        ratio_12_23 = bcs[1]['XX'] / bcs[2]['XX']
        print(f"  Δ={r['delta']:+.2f}: XX ratios (0-1)/(2-3) = {ratio_01_23:.2f}, (1-2)/(2-3) = {ratio_12_23:.2f}")

# ---- Save ----
output = {
    'experiment': '032c',
    'description': 'BW fidelity across XXZ phase diagram',
    'parameters': {'n': n, 'n_A': n_A, 'n_points': len(delta_values)},
    'sweep_results': results_list,
    'max_locality': {'delta': float(max_loc_xxz['delta']), 'value': float(max_loc_xxz['locality_fraction'])},
    'min_locality': {'delta': float(min_loc_xxz['delta']), 'value': float(min_loc_xxz['locality_fraction'])},
    'runtime_seconds': time.time() - t0
}

with open('results/sprint_032c_bw_fidelity_xxz_sweep.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved. Runtime: {time.time()-t0:.1f}s")
