"""
Sprint 021c: Holevo information vs T2/T1 ratio at fixed total noise
Physical model: T1 = 1/gamma, T2 = 1/(gamma/2 + lambda_pure)
  => lambda_pure = 1/T2 - gamma/2
  => for fixed total "noise time" t, gamma_eff = 1 - exp(-t/T1), lambda_eff = 1 - exp(-t/T_phi)
We parametrize by: noise_strength p (controls gamma), and ratio r = T2/T1 in [0, 2]
  gamma = p, lambda_pure = p*(1/r - 1/2) approximately (small p limit)
For simplicity: sweep lambda from 0 to p (matching physical constraint T2 >= T1/2 roughly)

Key question: where in the T2/T1 space does each code perform best, averaged over bases?
"""

import numpy as np
import json, time

start = time.time()

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])

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

def avg_basis_holevo(sv_0z, sv_1z, gamma, lam, n):
    """Compute average Holevo over Z, X, Y bases."""
    sv_0x = (sv_0z + sv_1z) / np.sqrt(2)
    sv_1x = (sv_0z - sv_1z) / np.sqrt(2)
    sv_0y = (sv_0z + 1j * sv_1z) / np.sqrt(2)
    sv_1y = (sv_0z - 1j * sv_1z) / np.sqrt(2)

    total = 0.0
    basis_vals = {}
    for basis_name, (s0, s1) in [('Z', (sv_0z, sv_1z)), ('X', (sv_0x, sv_1x)), ('Y', (sv_0y, sv_1y))]:
        rho_0 = density_matrix(s0)
        rho_1 = density_matrix(s1)
        rho_0n = apply_combined_noise(rho_0, gamma, lam, n)
        rho_1n = apply_combined_noise(rho_1, gamma, lam, n)
        h = holevo_info(rho_0n, rho_1n)
        total += h
        basis_vals[basis_name] = float(h)
    return total / 3, basis_vals

# ============================================================
# Build code states
# ============================================================

# 3-qubit
n3 = 3; dim3 = 2**n3
sv_0z_3q = np.zeros(dim3); sv_0z_3q[0] = 1.0
sv_1z_3q = np.zeros(dim3); sv_1z_3q[7] = 1.0

# [[5,1,3]]
n5 = 5; dim5 = 2**n5
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

# Toric [[8,2,2]]
n8 = 8; dim8 = 2**n8
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
sv_0z_toric = ground_states['++']
sv_1z_toric = ground_states['-+']

# Uncoded
sv_0z_unc = np.array([1.0, 0.0])
sv_1z_unc = np.array([0.0, 1.0])

print(f"Code states built in {time.time()-start:.1f}s")

# ============================================================
# Sweep T2/T1 ratio at fixed noise strengths
# ============================================================
print("\n=== Average-Basis Holevo vs Dephasing Fraction ===")
print("Parametrize: gamma fixed, lambda swept from 0 to gamma")
print("(lambda=0 means T2=2*T1, lambda=gamma means strong additional dephasing)\n")

noise_strengths = [0.05, 0.10, 0.20, 0.30]
# For each gamma, sweep lambda from 0 to 2*gamma (beyond physical bound for completeness)
lambda_fractions = [0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0]

codes_info = {
    'uncoded': (sv_0z_unc, sv_1z_unc, 1),
    '3-qubit': (sv_0z_3q, sv_1z_3q, n3),
    '[[5,1,3]]': (sv_0z_513, sv_1z_513, n5),
    'toric': (sv_0z_toric, sv_1z_toric, n8),
}

all_results = {}
for gamma in noise_strengths:
    print(f"\n--- γ = {gamma:.2f} ---")
    print(f"{'lambda/gamma':>12s} {'lambda':>8s} | {'uncoded':>8s} {'3-qubit':>8s} {'[[5,1,3]]':>10s} {'toric':>8s} | WINNER")
    all_results[str(gamma)] = {}

    for frac in lambda_fractions:
        lam = gamma * frac
        if lam > 0.95:  # cap at physical limit
            continue
        key = f"{frac:.2f}"
        all_results[str(gamma)][key] = {}

        winner_name = 'uncoded'
        winner_val = 0
        for code_name, (sv0, sv1, n) in codes_info.items():
            avg_h, basis_vals = avg_basis_holevo(sv0, sv1, gamma, lam, n)
            all_results[str(gamma)][key][code_name] = {
                'avg_holevo': round(avg_h, 4),
                'Z': round(basis_vals['Z'], 4),
                'X': round(basis_vals['X'], 4),
                'Y': round(basis_vals['Y'], 4),
                'asymmetry': round(max(basis_vals.values()) - min(basis_vals.values()), 4),
            }
            if avg_h > winner_val:
                winner_val = avg_h
                winner_name = code_name

        all_results[str(gamma)][key]['winner'] = winner_name
        vals = [all_results[str(gamma)][key][c]['avg_holevo'] for c in ['uncoded', '3-qubit', '[[5,1,3]]', 'toric']]
        print(f"  {frac:>10.2f}  {lam:>8.3f} | {vals[0]:8.3f} {vals[1]:8.3f} {vals[2]:10.3f} {vals[3]:8.3f} | {winner_name}")

# ============================================================
# Key analysis: crossover detection
# ============================================================
print("\n\n=== Crossover Analysis ===\n")
print("Finding lambda/gamma ratios where code rankings change:\n")

for gamma in noise_strengths:
    print(f"γ = {gamma:.2f}:")
    prev_winner = None
    for frac in lambda_fractions:
        lam = gamma * frac
        if lam > 0.95:
            continue
        key = f"{frac:.2f}"
        w = all_results[str(gamma)][key]['winner']
        if w != prev_winner and prev_winner is not None:
            print(f"  Crossover at λ/γ ≈ {frac:.2f} (λ={lam:.3f}): {prev_winner} → {w}")
        prev_winner = w

# ============================================================
# "Perfect code" metric: average Holevo * (1 - asymmetry)
# ============================================================
print("\n\n=== Code Quality Score: avg_Holevo * (1 - asymmetry) ===")
print("(Rewards both high info AND basis-isotropy)\n")

for gamma in [0.10, 0.20, 0.30]:
    print(f"γ = {gamma:.2f}:")
    for frac in [0.0, 0.50, 1.0, 2.0]:
        lam = gamma * frac
        if lam > 0.95:
            continue
        key = f"{frac:.2f}"
        if key not in all_results[str(gamma)]:
            continue
        print(f"  λ/γ={frac:.1f} (λ={lam:.3f}):", end="")
        for code_name in ['uncoded', '3-qubit', '[[5,1,3]]', 'toric']:
            d = all_results[str(gamma)][key][code_name]
            score = d['avg_holevo'] * (1 - d['asymmetry'])
            print(f"  {code_name}={score:.3f}", end="")
        print()

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save
results = {
    'sprint': '021c',
    'description': 'Holevo vs T2/T1 ratio at fixed noise strengths',
    'noise_strengths': noise_strengths,
    'lambda_fractions': lambda_fractions,
    'results': all_results,
    'elapsed_seconds': round(elapsed, 1)
}

with open('results/sprint_021c_t2_t1_ratio.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_021c_t2_t1_ratio.json")
