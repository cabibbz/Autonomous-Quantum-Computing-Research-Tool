"""
Sprint 020c: Toric code Page curve and noise fingerprint
1. Page curve: MI(R:subset) for logical information recovery vs subset size
2. Noise fingerprint: depolarizing, amplitude damping, phase damping response
Compare with [[5,1,3]] (Sprint 015) and code archetypes.
"""

import numpy as np
from itertools import combinations
import json, time

start = time.time()

def density_matrix(sv):
    sv = sv.reshape(-1, 1)
    return sv @ sv.conj().T

def partial_trace(rho, keep, n):
    d = 2
    rho = rho.reshape([d] * (2 * n))
    trace_out = sorted(set(range(n)) - set(keep))
    for i, q in enumerate(sorted(trace_out)):
        q_adj = q - i
        n_cur = n - i
        rho = np.trace(rho, axis1=q_adj, axis2=q_adj + n_cur)
    return rho.reshape(d**len(keep), d**len(keep))

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return -np.sum(evals * np.log2(evals))

def fidelity(rho, sigma):
    """Quantum fidelity between two density matrices."""
    sqrt_rho = np.linalg.cholesky(rho + 1e-14 * np.eye(len(rho)))
    M = sqrt_rho @ sigma @ sqrt_rho.conj().T
    evals = np.linalg.eigvalsh(M)
    evals = np.maximum(evals, 0)
    return (np.sum(np.sqrt(evals)))**2

n = 8
dim = 2**n
I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Z = np.array([[1, 0], [0, -1]])

def pauli_string(qubits, pauli, n):
    result = np.eye(1)
    for q in range(n):
        result = np.kron(result, pauli if q in qubits else I2)
    return result

# Build toric code states
vertex_stabs = [[0, 1, 4, 6], [0, 1, 5, 7], [2, 3, 4, 6]]
plaquette_stabs = [[0, 2, 4, 5], [1, 3, 4, 5], [0, 2, 6, 7]]

proj = np.eye(dim)
for qubits in vertex_stabs:
    S = pauli_string(qubits, X, n)
    proj = proj @ (np.eye(dim) + S) / 2
for qubits in plaquette_stabs:
    S = pauli_string(qubits, Z, n)
    proj = proj @ (np.eye(dim) + S) / 2

# Get all 4 ground states
eigenvalues, eigenvectors = np.linalg.eigh(proj)
code_basis_idx = np.where(eigenvalues > 0.5)[0]
code_basis = eigenvectors[:, code_basis_idx]

Z_L1 = pauli_string([0, 1], Z, n)
Z_L2 = pauli_string([4, 6], Z, n)
Z_L1_code = code_basis.conj().T @ Z_L1 @ code_basis
Z_L2_code = code_basis.conj().T @ Z_L2 @ code_basis

evals1, evecs1 = np.linalg.eigh(Z_L1_code)
ground_states = {}
for eval_target in [1.0, -1.0]:
    idx = np.where(np.abs(evals1 - eval_target) < 0.1)[0]
    sub_evecs = evecs1[:, idx]
    Z_L2_sub = sub_evecs.conj().T @ Z_L2_code @ sub_evecs
    evals2, evecs2 = np.linalg.eigh(Z_L2_sub)
    for j, ev2 in enumerate(evals2):
        vec_code = sub_evecs @ evecs2[:, j]
        vec_full = code_basis @ vec_code
        z1_sign = '+' if eval_target > 0 else '-'
        z2_sign = '+' if ev2 > 0 else '-'
        ground_states[f"{z1_sign}{z2_sign}"] = vec_full

# ============================================================
# Page Curve: Holevo information vs subset size
# ============================================================
print("=== Page Curve: Logical Information Recovery ===\n")

# Encode a logical state and see how much information a subset can recover
# Use the reference state approach: encode |0>_L1 and |1>_L1 (first logical qubit)
# with second logical qubit in |+>_L2

# States for logical qubit 1: |0>_L1 = (|++> + |+->)/sqrt(2), |1>_L1 = (|-+> + |-->)/sqrt(2)
# (averaging over logical qubit 2)
sv_0L = (ground_states['++'] + ground_states['+-']) / np.sqrt(2)
sv_1L = (ground_states['-+'] + ground_states['--']) / np.sqrt(2)
rho_0L = density_matrix(sv_0L)
rho_1L = density_matrix(sv_1L)

# Holevo information: chi = S(rho_avg) - (S(rho_0) + S(rho_1))/2
# where rho_avg = (rho_0 + rho_1)/2
# Computed for each subset
print("Toric code [[8,2,2]] Page curve (logical qubit 1):")
page_curve_toric = {}
for size in range(1, n):
    holevo_values = []
    for subset in combinations(range(n), size):
        sub = list(subset)
        rdm_0 = partial_trace(rho_0L, sub, n)
        rdm_1 = partial_trace(rho_1L, sub, n)
        rdm_avg = (rdm_0 + rdm_1) / 2
        chi = von_neumann_entropy(rdm_avg) - (von_neumann_entropy(rdm_0) + von_neumann_entropy(rdm_1)) / 2
        holevo_values.append(chi)
    page_curve_toric[size] = {
        'mean': round(np.mean(holevo_values), 4),
        'min': round(np.min(holevo_values), 4),
        'max': round(np.max(holevo_values), 4),
        'std': round(np.std(holevo_values), 4),
        'count': len(holevo_values)
    }
    print(f"  |A|={size}: Holevo={np.mean(holevo_values):.4f} ± {np.std(holevo_values):.4f} [{np.min(holevo_values):.4f}, {np.max(holevo_values):.4f}] (n={len(holevo_values)})")

# Also do [[5,1,3]] for comparison
print("\n[[5,1,3]] Page curve (from Sprint 015 setup):")
n_513 = 5
dim_513 = 2**n_513
Y = np.array([[0, -1j], [1j, 0]])

def multi_pauli_513(ops_dict, n):
    result = np.eye(1)
    for q in range(n):
        result = np.kron(result, ops_dict.get(q, I2))
    return result

stabs_513 = [
    {0: X, 1: Z, 2: Z, 3: X},
    {1: X, 2: Z, 3: Z, 4: X},
    {0: X, 2: X, 3: Z, 4: Z},
    {0: Z, 1: X, 3: X, 4: Z},
]
proj_513 = np.eye(dim_513)
for s in stabs_513:
    S = multi_pauli_513(s, n_513)
    proj_513 = proj_513 @ (np.eye(dim_513) + S) / 2

# Logical Z for [[5,1,3]]: ZZZZZ
Z_L_513 = multi_pauli_513({q: Z for q in range(5)}, 5)

# Get the two logical states
state_0 = np.zeros(dim_513); state_0[0] = 1.0
sv_0_513 = proj_513 @ state_0
sv_0_513 = sv_0_513 / np.linalg.norm(sv_0_513)

# |1>_L by applying logical X (= XXXXX)
X_L_513 = multi_pauli_513({q: X for q in range(5)}, 5)
sv_1_513 = X_L_513 @ sv_0_513
sv_1_513 = sv_1_513 / np.linalg.norm(sv_1_513)

rho_0_513 = density_matrix(sv_0_513)
rho_1_513 = density_matrix(sv_1_513)

page_curve_513 = {}
for size in range(1, n_513):
    holevo_values = []
    for subset in combinations(range(n_513), size):
        sub = list(subset)
        rdm_0 = partial_trace(rho_0_513, sub, n_513)
        rdm_1 = partial_trace(rho_1_513, sub, n_513)
        rdm_avg = (rdm_0 + rdm_1) / 2
        chi = von_neumann_entropy(rdm_avg) - (von_neumann_entropy(rdm_0) + von_neumann_entropy(rdm_1)) / 2
        holevo_values.append(chi)
    page_curve_513[size] = {
        'mean': round(np.mean(holevo_values), 4),
        'min': round(np.min(holevo_values), 4),
        'max': round(np.max(holevo_values), 4),
        'std': round(np.std(holevo_values), 4)
    }
    print(f"  |A|={size}: Holevo={np.mean(holevo_values):.4f} ± {np.std(holevo_values):.4f} [{np.min(holevo_values):.4f}, {np.max(holevo_values):.4f}]")

# ============================================================
# Noise Fingerprint
# ============================================================
print("\n=== Noise Fingerprint ===\n")

def apply_depolarizing(rho, p, n):
    """Apply single-qubit depolarizing noise to all qubits."""
    Y_mat = np.array([[0, -1j], [1j, 0]])
    paulis = [I2, X, Y_mat, Z]
    result = rho.astype(complex)
    for q in range(n):
        rho_new = np.zeros_like(result)
        for P in paulis:
            P_full = np.eye(1, dtype=complex)
            for i in range(n):
                P_full = np.kron(P_full, P if i == q else I2)
            if P is I2:
                rho_new += (1 - 3*p/4) * result
            else:
                rho_new += (p/4) * P_full @ result @ P_full.conj().T
        result = rho_new
    return result

def apply_amplitude_damping(rho, gamma, n):
    """Apply amplitude damping to all qubits."""
    K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]])
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]])
    result = rho.copy()
    for q in range(n):
        rho_new = np.zeros_like(result)
        for K in [K0, K1]:
            K_full = np.eye(1)
            for i in range(n):
                K_full = np.kron(K_full, K if i == q else I2)
            rho_new += K_full @ result @ K_full.conj().T
        result = rho_new
    return result

def apply_phase_damping(rho, lam, n):
    """Apply phase damping to all qubits."""
    K0 = np.array([[1, 0], [0, np.sqrt(1 - lam)]])
    K1 = np.array([[0, 0], [0, np.sqrt(lam)]])
    result = rho.copy()
    for q in range(n):
        rho_new = np.zeros_like(result)
        for K in [K0, K1]:
            K_full = np.eye(1)
            for i in range(n):
                K_full = np.kron(K_full, K if i == q else I2)
            rho_new += K_full @ result @ K_full.conj().T
        result = rho_new
    return result

# Compute Holevo information under noise for the toric code
# Using logical |0>_L and |1>_L (for first logical qubit, second in |0>_L)
sv_0L_clean = ground_states['++']
sv_1L_clean = ground_states['-+']
rho_0L_clean = density_matrix(sv_0L_clean)
rho_1L_clean = density_matrix(sv_1L_clean)

noise_levels = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]

noise_types = {
    'depolarizing': apply_depolarizing,
    'amplitude_damping': apply_amplitude_damping,
    'phase_damping': apply_phase_damping
}

noise_results = {}
for noise_name, noise_fn in noise_types.items():
    print(f"\n{noise_name}:")
    noise_results[noise_name] = {}
    for p in noise_levels:
        rho_0_noisy = noise_fn(rho_0L_clean, p, n)
        rho_1_noisy = noise_fn(rho_1L_clean, p, n)
        rho_avg_noisy = (rho_0_noisy + rho_1_noisy) / 2

        # Holevo information (full system)
        chi = von_neumann_entropy(rho_avg_noisy) - (von_neumann_entropy(rho_0_noisy) + von_neumann_entropy(rho_1_noisy)) / 2

        # Fidelity with clean state
        f0 = np.real(np.trace(rho_0L_clean @ rho_0_noisy))
        f1 = np.real(np.trace(rho_1L_clean @ rho_1_noisy))
        avg_fidelity = (f0 + f1) / 2

        # Codespace retention
        retention = np.real(np.trace(proj @ rho_0_noisy))

        noise_results[noise_name][str(p)] = {
            'holevo': round(chi, 4),
            'avg_fidelity': round(avg_fidelity, 4),
            'codespace_retention': round(retention, 4)
        }
        print(f"  p={p:.2f}: Holevo={chi:.4f}, Fidelity={avg_fidelity:.4f}, Retention={retention:.4f}")

# ============================================================
# Compare noise response: toric vs [[5,1,3]]
# ============================================================
print("\n=== [[5,1,3]] Noise Comparison ===")
noise_results_513 = {}
for noise_name, noise_fn in noise_types.items():
    print(f"\n{noise_name}:")
    noise_results_513[noise_name] = {}
    for p in noise_levels:
        rho_0_noisy = noise_fn(rho_0_513, p, n_513)
        rho_1_noisy = noise_fn(rho_1_513, p, n_513)
        rho_avg_noisy = (rho_0_noisy + rho_1_noisy) / 2

        chi = von_neumann_entropy(rho_avg_noisy) - (von_neumann_entropy(rho_0_noisy) + von_neumann_entropy(rho_1_noisy)) / 2
        f0 = np.real(np.trace(rho_0_513 @ rho_0_noisy))
        f1 = np.real(np.trace(rho_1_513 @ rho_1_noisy))
        avg_fidelity = (f0 + f1) / 2
        retention = np.real(np.trace(proj_513 @ rho_0_noisy))

        noise_results_513[noise_name][str(p)] = {
            'holevo': round(chi, 4),
            'avg_fidelity': round(avg_fidelity, 4),
            'codespace_retention': round(retention, 4)
        }
        print(f"  p={p:.2f}: Holevo={chi:.4f}, Fidelity={avg_fidelity:.4f}, Retention={retention:.4f}")

# Summary comparison
print("\n=== Summary: Noise Break-Even Points ===")
print("(Where Holevo < 0.5, i.e., less than half the logical information survives)")
for noise_name in noise_types:
    for code_name, nr in [("Toric [[8,2,2]]", noise_results), ("[[5,1,3]]", noise_results_513)]:
        for p in noise_levels:
            if nr[noise_name][str(p)]['holevo'] < 0.5:
                print(f"  {code_name} under {noise_name}: break-even at p≈{p:.2f}")
                break

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save results
results = {
    'sprint': '020c',
    'description': 'Toric code Page curve and noise fingerprint',
    'page_curve': {
        'toric': page_curve_toric,
        'five_qubit': page_curve_513
    },
    'noise_fingerprint': {
        'toric': noise_results,
        'five_qubit': noise_results_513
    },
    'elapsed_seconds': round(elapsed, 1)
}

with open('results/sprint_020c_toric_page_noise.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_020c_toric_page_noise.json")
