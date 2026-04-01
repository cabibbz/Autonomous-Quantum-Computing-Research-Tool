"""
Sprint 032b: BW Fidelity Across TFIM Phase Diagram

Sweep h/J from 0.1 to 3.0, compute BW fidelity at each point.
Test multiple BW envelope forms to find best fit.
Prediction: BW should work best at criticality (CFT applies).

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
    Z = np.array([[1,0],[0,-1]])
    return I, X, Z

def kron_list(ops):
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result

def tfim_hamiltonian(n, h_over_J):
    I, X, Z = pauli_matrices()
    dim = 2**n
    H = np.zeros((dim, dim))
    for i in range(n-1):
        ops = [I]*n; ops[i] = Z; ops[i+1] = Z
        H -= kron_list(ops)
    for i in range(n):
        ops = [I]*n; ops[i] = X
        H -= h_over_J * kron_list(ops)
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

# ---- BW Hamiltonians with different envelopes ----
def bw_hamiltonian(n_A, h_over_J, envelope_type='linear'):
    """
    Build BW prediction with different envelope choices.
    Subsystem A = left n_A qubits of a chain of length 2*n_A.
    Site 0 = far from cut, site n_A-1 = adjacent to cut.

    Envelopes (distance from cut = n_A - 1 - i for site i):
    - 'linear': β(i) ∝ (n_A - i)  [half-infinite BW]
    - 'sin_L': β(i) ∝ sin(π(i+0.5)/L)  [finite periodic CFT]
    - 'sin_inv': β(i) ∝ sin(π(n_A-i-0.5)/L)  [finite, from cut side]
    - 'uniform': β(i) = 1  [no position dependence]
    - 'log': β(i) ∝ log(n_A - i + 1)  [logarithmic]
    """
    I2, X, Z = pauli_matrices()
    dim = 2**n_A
    L = 2 * n_A
    H_BW = np.zeros((dim, dim))

    def get_beta_bond(i, etype):
        """β for bond between site i and i+1. Distance from cut = n_A - 1.5 - i."""
        d = n_A - 1 - i  # distance of bond midpoint from cut
        if etype == 'linear':
            return d + 0.5
        elif etype == 'sin_L':
            return np.sin(np.pi * (i + 1) / L)
        elif etype == 'sin_inv':
            return np.sin(np.pi * (n_A - i - 1) / L)
        elif etype == 'uniform':
            return 1.0
        elif etype == 'log':
            return np.log(d + 1.5)

    def get_beta_site(i, etype):
        """β for site i. Distance from cut = n_A - 1 - i."""
        d = n_A - 1 - i
        if etype == 'linear':
            return d + 0.5
        elif etype == 'sin_L':
            return np.sin(np.pi * (i + 0.5) / L)
        elif etype == 'sin_inv':
            return np.sin(np.pi * (n_A - i - 0.5) / L)
        elif etype == 'uniform':
            return 1.0
        elif etype == 'log':
            return np.log(d + 1.5)

    # ZZ terms
    for i in range(n_A - 1):
        beta = get_beta_bond(i, envelope_type)
        ops = [I2]*n_A; ops[i] = Z; ops[i+1] = Z
        H_BW -= beta * kron_list(ops)

    # X terms
    for i in range(n_A):
        beta = get_beta_site(i, envelope_type)
        ops = [I2]*n_A; ops[i] = X
        H_BW -= h_over_J * beta * kron_list(ops)

    return H_BW

def compute_bw_metrics(H_E, H_BW):
    """Compute BW fidelity and variance captured."""
    H_E_tl = H_E - np.trace(H_E)/H_E.shape[0] * np.eye(H_E.shape[0])
    H_BW_tl = H_BW - np.trace(H_BW)/H_BW.shape[0] * np.eye(H_BW.shape[0])

    fidelity = matrix_overlap(H_E_tl, H_BW_tl)

    alpha = np.real(np.trace(H_E_tl @ H_BW_tl)) / np.real(np.trace(H_BW_tl @ H_BW_tl))
    residual = H_E_tl - alpha * H_BW_tl
    resid_frac = la.norm(residual, 'fro') / la.norm(H_E_tl, 'fro')
    var_captured = 1 - resid_frac**2

    return fidelity, var_captured, alpha

def compute_locality_fraction(H_E, n_A, h_over_J):
    """Fraction of H_E norm² in TFIM-type terms (excluding identity)."""
    I2 = np.eye(2)
    X = np.array([[0,1],[1,0]])
    Z = np.array([[1,0],[0,-1]])

    # Build individual TFIM operators in subsystem A
    tfim_ops = []
    # ZZ bonds
    for i in range(n_A - 1):
        ops = [I2]*n_A; ops[i] = Z; ops[i+1] = Z
        tfim_ops.append(kron_list(ops))
    # X sites
    for i in range(n_A):
        ops = [I2]*n_A; ops[i] = X
        tfim_ops.append(kron_list(ops))

    dim = 2**n_A
    H_E_tl = H_E - np.trace(H_E)/dim * np.eye(dim)
    total_norm_sq = la.norm(H_E_tl, 'fro')**2

    # Project H_E onto TFIM subspace
    tfim_norm_sq = 0
    for op in tfim_ops:
        coeff = np.real(np.trace(H_E_tl @ op)) / dim
        tfim_norm_sq += coeff**2

    # Each Pauli has norm² = dim, so coefficient² * dim = contribution to Frobenius norm²
    # Actually: ||Σ c_i P_i||² = Σ c_i² * Tr(P_i²) = Σ c_i² * dim
    # And total_norm_sq = ||H_E_tl||² = Σ c_i² * dim (sum over all non-identity Paulis)
    # So fraction = (Σ c_TFIM² * dim) / (Σ c_all² * dim) = Σ c_TFIM² / Σ c_all²

    return tfim_norm_sq * dim / total_norm_sq if total_norm_sq > 0 else 0

# ---- Main sweep ----
n = 8
n_A = 4
h_J_values = np.concatenate([
    np.linspace(0.1, 0.8, 8),    # ordered phase
    np.linspace(0.85, 1.15, 7),  # critical region
    np.linspace(1.2, 3.0, 10),   # disordered phase
])

envelope_types = ['linear', 'sin_inv', 'uniform', 'log']

print(f"=== Sprint 032b: BW Fidelity Across TFIM Phase Diagram ===")
print(f"n={n}, n_A={n_A}, {len(h_J_values)} h/J values, {len(envelope_types)} envelope types")

results_list = []

for h_J in h_J_values:
    # Ground state
    H = tfim_hamiltonian(n, h_J)
    eigvals, eigvecs = la.eigh(H)
    psi_gs = eigvecs[:, 0]
    gap = eigvals[1] - eigvals[0]

    # Reduced density matrix and H_E
    rho = np.outer(psi_gs, psi_gs.conj())
    rho_A = partial_trace(rho, n, list(range(n_A)))
    H_E = entanglement_hamiltonian(rho_A)

    # Entanglement entropy
    ev = la.eigvalsh(rho_A)
    ev = ev[ev > 1e-15]
    S = -np.sum(ev * np.log2(ev))

    # Locality fraction
    loc_frac = compute_locality_fraction(H_E, n_A, h_J)

    # BW fidelity for each envelope
    env_results = {}
    for etype in envelope_types:
        H_BW = bw_hamiltonian(n_A, h_J, etype)
        fid, var_cap, alpha = compute_bw_metrics(H_E, H_BW)
        env_results[etype] = {'fidelity': float(fid), 'variance_captured': float(var_cap), 'alpha': float(alpha)}

    best_env = max(envelope_types, key=lambda e: env_results[e]['variance_captured'])

    result = {
        'h_over_J': float(h_J),
        'energy_gap': float(gap),
        'entropy_bits': float(S),
        'locality_fraction': float(loc_frac),
        'envelope_results': env_results,
        'best_envelope': best_env,
        'best_variance_captured': env_results[best_env]['variance_captured']
    }
    results_list.append(result)

    print(f"h/J={h_J:.2f}: gap={gap:.3f}, S={S:.3f}, locality={loc_frac:.3f}, "
          f"best={best_env}({env_results[best_env]['variance_captured']:.3f}), "
          f"linear={env_results['linear']['variance_captured']:.3f}")

# ---- Analysis ----
print("\n=== Analysis ===")

# Find where locality is highest
max_loc = max(results_list, key=lambda r: r['locality_fraction'])
print(f"\nHighest locality: {max_loc['locality_fraction']:.4f} at h/J={max_loc['h_over_J']:.2f}")

# Find where BW fidelity is highest
max_bw = max(results_list, key=lambda r: r['best_variance_captured'])
print(f"Highest BW variance: {max_bw['best_variance_captured']:.4f} at h/J={max_bw['h_over_J']:.2f} ({max_bw['best_envelope']})")

# Phase-averaged results
ordered = [r for r in results_list if r['h_over_J'] < 0.9]
critical = [r for r in results_list if 0.9 <= r['h_over_J'] <= 1.1]
disordered = [r for r in results_list if r['h_over_J'] > 1.1]

for phase_name, phase_data in [('Ordered', ordered), ('Critical', critical), ('Disordered', disordered)]:
    if not phase_data:
        continue
    avg_loc = np.mean([r['locality_fraction'] for r in phase_data])
    avg_bw = np.mean([r['best_variance_captured'] for r in phase_data])
    print(f"\n{phase_name} phase (n={len(phase_data)}):")
    print(f"  Avg locality: {avg_loc:.4f}")
    print(f"  Avg best BW variance: {avg_bw:.4f}")
    for etype in envelope_types:
        avg_var = np.mean([r['envelope_results'][etype]['variance_captured'] for r in phase_data])
        print(f"  Avg {etype} variance: {avg_var:.4f}")

# Envelope comparison
print("\n--- Envelope Comparison (wins across all points) ---")
for etype in envelope_types:
    wins = sum(1 for r in results_list if r['best_envelope'] == etype)
    avg_var = np.mean([r['envelope_results'][etype]['variance_captured'] for r in results_list])
    print(f"  {etype}: {wins}/{len(results_list)} wins, avg variance captured = {avg_var:.4f}")

# ---- Save ----
output = {
    'experiment': '032b',
    'description': 'BW fidelity across TFIM phase diagram',
    'parameters': {'n': n, 'n_A': n_A, 'n_points': len(h_J_values)},
    'envelope_types': envelope_types,
    'sweep_results': results_list,
    'max_locality': {'h_over_J': float(max_loc['h_over_J']), 'value': float(max_loc['locality_fraction'])},
    'max_bw_variance': {'h_over_J': float(max_bw['h_over_J']), 'value': float(max_bw['best_variance_captured']),
                        'envelope': max_bw['best_envelope']},
    'phase_averages': {
        'ordered': {'avg_locality': float(np.mean([r['locality_fraction'] for r in ordered])),
                    'avg_best_bw': float(np.mean([r['best_variance_captured'] for r in ordered]))},
        'critical': {'avg_locality': float(np.mean([r['locality_fraction'] for r in critical])),
                     'avg_best_bw': float(np.mean([r['best_variance_captured'] for r in critical]))},
        'disordered': {'avg_locality': float(np.mean([r['locality_fraction'] for r in disordered])),
                       'avg_best_bw': float(np.mean([r['best_variance_captured'] for r in disordered]))}
    },
    'runtime_seconds': time.time() - t0
}

with open('results/sprint_032b_bw_fidelity_tfim_sweep.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved. Runtime: {time.time()-t0:.1f}s")
