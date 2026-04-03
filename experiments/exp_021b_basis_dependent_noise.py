"""
Sprint 021b: Basis-dependent code performance under combined T1+T2 noise
The 3-qubit code appeared to dominate in 21a, but that used Z-basis logical states
which are phase-damping eigenstates. Test with X-basis and Y-basis logical states.

Key question: Does the 3-qubit advantage survive when logical information is NOT
aligned with the noise channel's preferred basis?
"""

import numpy as np
import json, time

start = time.time()

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])
H_gate = np.array([[1, 1], [1, -1]]) / np.sqrt(2)

def density_matrix(sv):
    sv = sv.reshape(-1, 1)
    return sv @ sv.conj().T

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return -np.sum(evals * np.log2(evals))

def multi_pauli(ops_dict, n):
    result = np.eye(1)
    for q in range(n):
        result = np.kron(result, ops_dict.get(q, I2))
    return result

def apply_amplitude_damping(rho, gamma, n):
    K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]])
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]])
    result = rho.copy().astype(complex)
    for q in range(n):
        rho_new = np.zeros_like(result)
        for K in [K0, K1]:
            K_full = np.eye(1, dtype=complex)
            for i in range(n):
                K_full = np.kron(K_full, K if i == q else I2)
            rho_new += K_full @ result @ K_full.conj().T
        result = rho_new
    return result

def apply_phase_damping(rho, lam, n):
    K0 = np.array([[1, 0], [0, np.sqrt(1 - lam)]])
    K1 = np.array([[0, 0], [0, np.sqrt(lam)]])
    result = rho.copy().astype(complex)
    for q in range(n):
        rho_new = np.zeros_like(result)
        for K in [K0, K1]:
            K_full = np.eye(1, dtype=complex)
            for i in range(n):
                K_full = np.kron(K_full, K if i == q else I2)
            rho_new += K_full @ result @ K_full.conj().T
        result = rho_new
    return result

def apply_combined_noise(rho, gamma, lam, n):
    rho = apply_amplitude_damping(rho, gamma, n)
    rho = apply_phase_damping(rho, lam, n)
    return rho

def holevo_info(rho_0, rho_1):
    rho_avg = (rho_0 + rho_1) / 2
    return von_neumann_entropy(rho_avg) - (von_neumann_entropy(rho_0) + von_neumann_entropy(rho_1)) / 2

# ============================================================
# Build code states in multiple bases
# ============================================================

# --- 3-qubit bit-flip code ---
n3 = 3
dim3 = 2**n3
# Z-basis: |0>_L = |000>, |1>_L = |111>
sv_0z_3q = np.zeros(dim3); sv_0z_3q[0] = 1.0
sv_1z_3q = np.zeros(dim3); sv_1z_3q[7] = 1.0
# X-basis: |+>_L = (|000>+|111>)/sqrt(2), |->_L = (|000>-|111>)/sqrt(2)
sv_0x_3q = (sv_0z_3q + sv_1z_3q) / np.sqrt(2)
sv_1x_3q = (sv_0z_3q - sv_1z_3q) / np.sqrt(2)
# Y-basis: |+i>_L = (|000>+i|111>)/sqrt(2), |-i>_L = (|000>-i|111>)/sqrt(2)
sv_0y_3q = (sv_0z_3q + 1j * sv_1z_3q) / np.sqrt(2)
sv_1y_3q = (sv_0z_3q - 1j * sv_1z_3q) / np.sqrt(2)

# --- [[5,1,3]] code ---
n5 = 5
dim5 = 2**n5
stabs_513 = [
    {0: X, 1: Z, 2: Z, 3: X},
    {1: X, 2: Z, 3: Z, 4: X},
    {0: X, 2: X, 3: Z, 4: Z},
    {0: Z, 1: X, 3: X, 4: Z},
]
proj_513 = np.eye(dim5)
for s in stabs_513:
    S = multi_pauli(s, n5)
    proj_513 = proj_513 @ (np.eye(dim5) + S) / 2

state_0 = np.zeros(dim5); state_0[0] = 1.0
sv_0z_513 = proj_513 @ state_0
sv_0z_513 = sv_0z_513 / np.linalg.norm(sv_0z_513)
X_L_513 = multi_pauli({q: X for q in range(5)}, 5)
sv_1z_513 = X_L_513 @ sv_0z_513
sv_1z_513 = sv_1z_513 / np.linalg.norm(sv_1z_513)
# X-basis logical states
sv_0x_513 = (sv_0z_513 + sv_1z_513) / np.sqrt(2)
sv_1x_513 = (sv_0z_513 - sv_1z_513) / np.sqrt(2)
# Y-basis logical states
sv_0y_513 = (sv_0z_513 + 1j * sv_1z_513) / np.sqrt(2)
sv_1y_513 = (sv_0z_513 - 1j * sv_1z_513) / np.sqrt(2)

# --- Toric [[8,2,2]] code ---
n8 = 8
dim8 = 2**n8

def pauli_string(qubits, pauli, n):
    result = np.eye(1)
    for q in range(n):
        result = np.kron(result, pauli if q in qubits else I2)
    return result

vertex_stabs = [[0, 1, 4, 6], [0, 1, 5, 7], [2, 3, 4, 6]]
plaquette_stabs = [[0, 2, 4, 5], [1, 3, 4, 5], [0, 2, 6, 7]]
proj_toric = np.eye(dim8)
for qubits in vertex_stabs:
    S = pauli_string(qubits, X, n8)
    proj_toric = proj_toric @ (np.eye(dim8) + S) / 2
for qubits in plaquette_stabs:
    S = pauli_string(qubits, Z, n8)
    proj_toric = proj_toric @ (np.eye(dim8) + S) / 2

eigenvalues, eigenvectors = np.linalg.eigh(proj_toric)
code_basis_idx = np.where(eigenvalues > 0.5)[0]
code_basis = eigenvectors[:, code_basis_idx]
Z_L1 = pauli_string([0, 1], Z, n8)
Z_L2 = pauli_string([4, 6], Z, n8)
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

# Toric Z-basis (logical qubit 1)
sv_0z_toric = ground_states['++']
sv_1z_toric = ground_states['-+']
# Toric X-basis
sv_0x_toric = (sv_0z_toric + sv_1z_toric) / np.sqrt(2)
sv_1x_toric = (sv_0z_toric - sv_1z_toric) / np.sqrt(2)
# Toric Y-basis
sv_0y_toric = (sv_0z_toric + 1j * sv_1z_toric) / np.sqrt(2)
sv_1y_toric = (sv_0z_toric - 1j * sv_1z_toric) / np.sqrt(2)

# --- Uncoded single qubit ---
sv_0z_unc = np.array([1.0, 0.0])
sv_1z_unc = np.array([0.0, 1.0])
sv_0x_unc = np.array([1.0, 1.0]) / np.sqrt(2)
sv_1x_unc = np.array([1.0, -1.0]) / np.sqrt(2)
sv_0y_unc = np.array([1.0, 1j]) / np.sqrt(2)
sv_1y_unc = np.array([1.0, -1j]) / np.sqrt(2)

print(f"Code states built in {time.time()-start:.1f}s")

# ============================================================
# Compare Z, X, Y bases under combined noise
# ============================================================
print("\n=== Holevo Information: Z vs X vs Y Basis ===\n")

# Sample a range of (gamma, lambda) points
test_points = [
    (0.00, 0.10), (0.00, 0.30), (0.00, 0.50),  # pure dephasing
    (0.10, 0.00), (0.30, 0.00), (0.50, 0.00),  # pure relaxation
    (0.10, 0.10), (0.20, 0.20), (0.30, 0.30),  # balanced
    (0.05, 0.30), (0.10, 0.40), (0.30, 0.10),  # asymmetric
]

codes = {
    'uncoded': {
        'n': 1,
        'Z': (sv_0z_unc, sv_1z_unc),
        'X': (sv_0x_unc, sv_1x_unc),
        'Y': (sv_0y_unc, sv_1y_unc),
    },
    '3-qubit': {
        'n': n3,
        'Z': (sv_0z_3q, sv_1z_3q),
        'X': (sv_0x_3q, sv_1x_3q),
        'Y': (sv_0y_3q, sv_1y_3q),
    },
    '[[5,1,3]]': {
        'n': n5,
        'Z': (sv_0z_513, sv_1z_513),
        'X': (sv_0x_513, sv_1x_513),
        'Y': (sv_0y_513, sv_1y_513),
    },
    'toric': {
        'n': n8,
        'Z': (sv_0z_toric, sv_1z_toric),
        'X': (sv_0x_toric, sv_1x_toric),
        'Y': (sv_0y_toric, sv_1y_toric),
    },
}

all_results = {}
for gamma, lam in test_points:
    key = f"({gamma:.2f},{lam:.2f})"
    all_results[key] = {}
    print(f"\nγ={gamma:.2f}, λ={lam:.2f}:")
    for code_name, code_data in codes.items():
        n = code_data['n']
        all_results[key][code_name] = {}
        for basis in ['Z', 'X', 'Y']:
            sv0, sv1 = code_data[basis]
            rho_0 = density_matrix(sv0)
            rho_1 = density_matrix(sv1)
            rho_0n = apply_combined_noise(rho_0, gamma, lam, n)
            rho_1n = apply_combined_noise(rho_1, gamma, lam, n)
            h = holevo_info(rho_0n, rho_1n)
            all_results[key][code_name][basis] = round(float(h), 4)
        z, x, y = all_results[key][code_name]['Z'], all_results[key][code_name]['X'], all_results[key][code_name]['Y']
        # Basis asymmetry: max - min
        asym = max(z, x, y) - min(z, x, y)
        all_results[key][code_name]['asymmetry'] = round(float(asym), 4)
        all_results[key][code_name]['average'] = round(float((z + x + y) / 3), 4)
        print(f"  {code_name:10s}: Z={z:.3f}  X={x:.3f}  Y={y:.3f}  avg={all_results[key][code_name]['average']:.3f}  asym={asym:.3f}")

# ============================================================
# Analysis: average-over-basis Holevo — fair comparison
# ============================================================
print("\n\n=== Average-Basis Holevo: Fair Code Comparison ===\n")
print(f"{'(gamma,lambda)':>16s} | {'uncoded':>8s} | {'3-qubit':>8s} | {'[[5,1,3]]':>10s} | {'toric':>8s} | WINNER")
print("-" * 80)
for gamma, lam in test_points:
    key = f"({gamma:.2f},{lam:.2f})"
    avgs = {name: all_results[key][name]['average'] for name in codes}
    winner = max(avgs, key=avgs.get)
    print(f"  ({gamma:.2f},{lam:.2f})     | {avgs['uncoded']:8.3f} | {avgs['3-qubit']:8.3f} | {avgs['[[5,1,3]]']:10.3f} | {avgs['toric']:8.3f} | {winner}")

# ============================================================
# Basis asymmetry analysis
# ============================================================
print("\n\n=== Basis Asymmetry (max Holevo - min Holevo across Z,X,Y) ===\n")
print(f"{'(gamma,lambda)':>16s} | {'uncoded':>8s} | {'3-qubit':>8s} | {'[[5,1,3]]':>10s} | {'toric':>8s}")
print("-" * 70)
for gamma, lam in test_points:
    key = f"({gamma:.2f},{lam:.2f})"
    asyms = {name: all_results[key][name]['asymmetry'] for name in codes}
    print(f"  ({gamma:.2f},{lam:.2f})     | {asyms['uncoded']:8.3f} | {asyms['3-qubit']:8.3f} | {asyms['[[5,1,3]]']:10.3f} | {asyms['toric']:8.3f}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save
results = {
    'sprint': '021b',
    'description': 'Basis-dependent code performance under combined T1+T2 noise',
    'test_points': [(g, l) for g, l in test_points],
    'results': all_results,
    'elapsed_seconds': round(elapsed, 1)
}

with open('results/sprint_021b_basis_dependent_noise.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_021b_basis_dependent_noise.json")
