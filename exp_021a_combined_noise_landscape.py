"""
Sprint 021a: Combined T1+T2 noise fidelity landscape
Sweep (gamma, lambda) for 3-qubit bit-flip, [[5,1,3]], and toric [[8,2,2]] codes.
Combined channel: amplitude damping(gamma) composed with phase damping(lambda).
Physical constraint: lambda <= 1 - sqrt(1-gamma) approximately (T2 <= 2T1 bound).
"""

import numpy as np
from itertools import combinations
import json, time

start = time.time()

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])

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
    """Apply amplitude damping then phase damping (order: T1 then T2)."""
    rho = apply_amplitude_damping(rho, gamma, n)
    rho = apply_phase_damping(rho, lam, n)
    return rho

def holevo_info(rho_0, rho_1):
    rho_avg = (rho_0 + rho_1) / 2
    return von_neumann_entropy(rho_avg) - (von_neumann_entropy(rho_0) + von_neumann_entropy(rho_1)) / 2

# ============================================================
# Build code states
# ============================================================

# --- 3-qubit bit-flip code ---
n3 = 3
dim3 = 2**n3
sv_0_3q = np.zeros(dim3); sv_0_3q[0] = 1.0  # |000>
sv_1_3q = np.zeros(dim3); sv_1_3q[7] = 1.0  # |111>
rho_0_3q = density_matrix(sv_0_3q)
rho_1_3q = density_matrix(sv_1_3q)

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
sv_0_513 = proj_513 @ state_0
sv_0_513 = sv_0_513 / np.linalg.norm(sv_0_513)
X_L_513 = multi_pauli({q: X for q in range(5)}, 5)
sv_1_513 = X_L_513 @ sv_0_513
sv_1_513 = sv_1_513 / np.linalg.norm(sv_1_513)
rho_0_513 = density_matrix(sv_0_513)
rho_1_513 = density_matrix(sv_1_513)

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

sv_0_toric = ground_states['++']
sv_1_toric = ground_states['-+']
rho_0_toric = density_matrix(sv_0_toric)
rho_1_toric = density_matrix(sv_1_toric)

print(f"Code states built in {time.time()-start:.1f}s")

# ============================================================
# Sweep combined noise landscape
# ============================================================
print("\n=== Combined Noise (gamma, lambda) Landscape ===\n")

# Use a grid: gamma from 0 to 0.5, lambda from 0 to 0.5
gammas = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]
lambdas = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]

results_3q = {}
results_513 = {}
results_toric = {}

for gamma in gammas:
    for lam in lambdas:
        key = f"({gamma:.2f},{lam:.2f})"

        # 3-qubit
        r0 = apply_combined_noise(rho_0_3q, gamma, lam, n3)
        r1 = apply_combined_noise(rho_1_3q, gamma, lam, n3)
        f0 = np.real(np.trace(rho_0_3q @ r0))
        f1 = np.real(np.trace(rho_1_3q @ r1))
        h = holevo_info(r0, r1)
        results_3q[key] = {'fidelity': round(float((f0+f1)/2), 4), 'holevo': round(float(h), 4)}

        # [[5,1,3]]
        r0 = apply_combined_noise(rho_0_513, gamma, lam, n5)
        r1 = apply_combined_noise(rho_1_513, gamma, lam, n5)
        f0 = np.real(np.trace(rho_0_513 @ r0))
        f1 = np.real(np.trace(rho_1_513 @ r1))
        h = holevo_info(r0, r1)
        results_513[key] = {'fidelity': round(float((f0+f1)/2), 4), 'holevo': round(float(h), 4)}

        # Toric
        r0 = apply_combined_noise(rho_0_toric, gamma, lam, n8)
        r1 = apply_combined_noise(rho_1_toric, gamma, lam, n8)
        f0 = np.real(np.trace(rho_0_toric @ r0))
        f1 = np.real(np.trace(rho_1_toric @ r1))
        h = holevo_info(r0, r1)
        results_toric[key] = {'fidelity': round(float((f0+f1)/2), 4), 'holevo': round(float(h), 4)}

        print(f"  γ={gamma:.2f}, λ={lam:.2f} | 3q: F={results_3q[key]['fidelity']:.3f} H={results_3q[key]['holevo']:.3f} | [[5,1,3]]: F={results_513[key]['fidelity']:.3f} H={results_513[key]['holevo']:.3f} | Toric: F={results_toric[key]['fidelity']:.3f} H={results_toric[key]['holevo']:.3f}")

# ============================================================
# Also compute uncoded single-qubit baseline
# ============================================================
print("\n=== Uncoded (single qubit) baseline ===")
rho_0_uncoded = np.array([[1,0],[0,0]], dtype=complex)
rho_1_uncoded = np.array([[0,0],[0,1]], dtype=complex)
results_uncoded = {}
for gamma in gammas:
    for lam in lambdas:
        key = f"({gamma:.2f},{lam:.2f})"
        r0 = apply_combined_noise(rho_0_uncoded, gamma, lam, 1)
        r1 = apply_combined_noise(rho_1_uncoded, gamma, lam, 1)
        f0 = np.real(np.trace(rho_0_uncoded @ r0))
        f1 = np.real(np.trace(rho_1_uncoded @ r1))
        h = holevo_info(r0, r1)
        results_uncoded[key] = {'fidelity': round(float((f0+f1)/2), 4), 'holevo': round(float(h), 4)}

# ============================================================
# Analysis: where does each code beat uncoded?
# ============================================================
print("\n=== Code vs Uncoded: Holevo advantage ===")
print("(Positive = code beats uncoded)\n")
for gamma in gammas:
    for lam in lambdas:
        key = f"({gamma:.2f},{lam:.2f})"
        h_unc = results_uncoded[key]['holevo']
        h_3q = results_3q[key]['holevo']
        h_513 = results_513[key]['holevo']
        h_tor = results_toric[key]['holevo']
        winner = 'uncoded'
        best_h = h_unc
        for name, h in [('3-qubit', h_3q), ('[[5,1,3]]', h_513), ('toric', h_tor)]:
            if h > best_h:
                best_h = h
                winner = name
        if gamma > 0 or lam > 0:
            print(f"  γ={gamma:.2f}, λ={lam:.2f} | uncoded={h_unc:.3f} | 3q={h_3q:.3f} | [[5,1,3]]={h_513:.3f} | toric={h_tor:.3f} | WINNER: {winner}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save
results = {
    'sprint': '021a',
    'description': 'Combined T1+T2 noise fidelity landscape',
    'gammas': gammas,
    'lambdas': lambdas,
    'codes': {
        '3-qubit': results_3q,
        '[[5,1,3]]': results_513,
        'toric': results_toric,
        'uncoded': results_uncoded
    },
    'elapsed_seconds': round(elapsed, 1)
}

with open('results/sprint_021a_combined_noise_landscape.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_021a_combined_noise_landscape.json")
