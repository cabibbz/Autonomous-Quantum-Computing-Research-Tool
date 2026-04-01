"""
Sprint 018b: Holevo information for [[5,1,3]] code under depolarizing noise.

Sprint 017 showed the threshold is a phase transition for the bit-flip
repetition code (Holevo info → step function with concatenation).
Does the [[5,1,3]] code, which corrects ALL single-qubit errors, show
the same phase transition? Compare its Holevo curve shape to the
repetition code.

Also: compare the [[5,1,3]] code's Holevo information to uncoded and
3-qubit code under depolarizing noise (apples-to-apples comparison).
"""

import numpy as np
import json

# Pauli matrices
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

def kron_list(mats):
    result = mats[0]
    for m in mats[1:]:
        result = np.kron(result, m)
    return result

def pauli_string(paulis, n):
    pmap = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
    return kron_list([pmap[c] for c in paulis])

def von_neumann_entropy(rho):
    """Von Neumann entropy in bits from density matrix (numpy array)."""
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-15]
    return float(-np.sum(evals * np.log2(evals)))

def depolarizing_channel(rho, p, n):
    """Apply depolarizing noise with probability p to each qubit independently."""
    for q in range(n):
        ops = [I2] * n
        ops_x, ops_y, ops_z = [I2]*n, [I2]*n, [I2]*n
        ops_x[q] = X; ops_y[q] = Y; ops_z[q] = Z
        XI = kron_list(ops_x); YI = kron_list(ops_y); ZI = kron_list(ops_z)
        rho = (1 - p) * rho + (p / 3) * (XI @ rho @ XI + YI @ rho @ YI + ZI @ rho @ ZI)
    return rho

def build_five_qubit_code():
    """Build encoding isometry for [[5,1,3]] code."""
    n = 5
    dim = 2**n
    stabilizers = ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']
    state_0 = np.zeros(dim, dtype=complex)
    state_0[0] = 1.0
    for stab in stabilizers:
        S = pauli_string(stab, n)
        proj = (np.eye(dim) + S) / 2.0
        state_0 = proj @ state_0
        norm = np.linalg.norm(state_0)
        if norm > 1e-15:
            state_0 = state_0 / norm
    X_all = pauli_string('XXXXX', n)
    state_1 = X_all @ state_0
    V = np.zeros((dim, 2), dtype=complex)
    V[:, 0] = state_0
    V[:, 1] = state_1
    return V

def five_qubit_syndrome_correct(rho, V):
    """Syndrome measurement and correction for [[5,1,3]] code."""
    n = 5
    dim = 2**n
    stabilizers = ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']
    stab_ops = [pauli_string(s, n) for s in stabilizers]
    rho_corrected = np.zeros_like(rho)

    for syndrome in range(16):
        proj = np.eye(dim, dtype=complex)
        for k, S in enumerate(stab_ops):
            sign = 1 if (syndrome >> k) & 1 == 0 else -1
            proj = (np.eye(dim) + sign * S) / 2.0 @ proj

        if syndrome == 0:
            correction = np.eye(dim)
        else:
            correction = np.eye(dim)
            for q in range(n):
                found = False
                for P in [X, Y, Z]:
                    error = np.eye(dim, dtype=complex)
                    ops = [I2] * n
                    ops[q] = P
                    error = kron_list(ops)
                    s = 0
                    for k, S in enumerate(stab_ops):
                        diff = S @ error - error @ S
                        if np.linalg.norm(diff) > 1e-10:
                            s |= (1 << k)
                    if s == syndrome:
                        correction = error
                        found = True
                        break
                if found:
                    break

        projected = proj @ rho @ proj.conj().T
        rho_corrected += correction @ projected @ correction.conj().T

    return rho_corrected

def three_qubit_encode(logical_rho):
    V = np.zeros((8, 2), dtype=complex)
    V[0b000, 0] = 1.0
    V[0b111, 1] = 1.0
    return V @ logical_rho @ V.conj().T

def three_qubit_syndrome(rho):
    n = 3
    Z0Z1 = kron_list([Z, Z, I2])
    Z1Z2 = kron_list([I2, Z, Z])
    P_pp = (np.eye(8) + Z0Z1) / 2 @ (np.eye(8) + Z1Z2) / 2
    P_pm = (np.eye(8) + Z0Z1) / 2 @ (np.eye(8) - Z1Z2) / 2
    P_mp = (np.eye(8) - Z0Z1) / 2 @ (np.eye(8) + Z1Z2) / 2
    P_mm = (np.eye(8) - Z0Z1) / 2 @ (np.eye(8) - Z1Z2) / 2
    X0 = kron_list([X, I2, I2])
    X1 = kron_list([I2, X, I2])
    X2 = kron_list([I2, I2, X])
    return (P_pp @ rho @ P_pp +
            X2 @ P_pm @ rho @ P_pm @ X2 +
            X0 @ P_mp @ rho @ P_mp @ X0 +
            X1 @ P_mm @ rho @ P_mm @ X1)

def three_qubit_decode(rho):
    V = np.zeros((8, 2), dtype=complex)
    V[0b000, 0] = 1.0
    V[0b111, 1] = 1.0
    return V.conj().T @ rho @ V

def holevo_information(rho0, rho1):
    """Holevo info = S(avg) - avg(S) for equal priors."""
    rho_avg = (rho0 + rho1) / 2
    return von_neumann_entropy(rho_avg) - 0.5 * von_neumann_entropy(rho0) - 0.5 * von_neumann_entropy(rho1)

def fidelity_with_pure(rho, state_vec):
    return float(np.real(state_vec.conj() @ rho @ state_vec))


# ============================================================
# Build [[5,1,3]] code
# ============================================================
V5 = build_five_qubit_code()
print("[[5,1,3]] code built.\n")

# Logical states
psi0 = np.array([1, 0], dtype=complex)
psi1 = np.array([0, 1], dtype=complex)
psi_plus = np.array([1, 1], dtype=complex) / np.sqrt(2)

p_values = [0.0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15,
            0.18, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.75]

# ============================================================
# 1. Uncoded single qubit under depolarizing
# ============================================================
print("=== UNCODED QUBIT (depolarizing) ===")
uncoded_results = {'holevo': [], 'fidelity_0': [], 'fidelity_plus': []}

for p in p_values:
    rho0 = np.outer(psi0, psi0.conj())
    rho1 = np.outer(psi1, psi1.conj())
    rho0_n = (1-p)*rho0 + (p/3)*(X@rho0@X + Y@rho0@Y + Z@rho0@Z)
    rho1_n = (1-p)*rho1 + (p/3)*(X@rho1@X + Y@rho1@Y + Z@rho1@Z)
    hol = holevo_information(rho0_n, rho1_n)
    f0 = fidelity_with_pure(rho0_n, psi0)

    rho_p = np.outer(psi_plus, psi_plus.conj())
    rho_p_n = (1-p)*rho_p + (p/3)*(X@rho_p@X + Y@rho_p@Y + Z@rho_p@Z)
    f_plus = fidelity_with_pure(rho_p_n, psi_plus)

    uncoded_results['holevo'].append(round(hol, 6))
    uncoded_results['fidelity_0'].append(round(f0, 6))
    uncoded_results['fidelity_plus'].append(round(f_plus, 6))

# ============================================================
# 2. 3-qubit bit-flip code under depolarizing (encode, noise, correct, decode)
# ============================================================
print("=== 3-QUBIT CODE (depolarizing) ===")
three_q_results = {'holevo': [], 'fidelity_0': [], 'fidelity_plus': []}

for p in p_values:
    print(f"  p={p:.2f}...")
    holevo_vals = []
    for logical_sv in [psi0, psi1]:
        logical_rho = np.outer(logical_sv, logical_sv.conj())
        rho_coded = three_qubit_encode(logical_rho)
        rho_noisy = depolarizing_channel(rho_coded, p, 3)
        rho_corrected = three_qubit_syndrome(rho_noisy)
        rho_decoded = three_qubit_decode(rho_corrected)
        tr = np.trace(rho_decoded).real
        if tr > 1e-12:
            rho_decoded = rho_decoded / tr
        holevo_vals.append(rho_decoded)

    hol = holevo_information(holevo_vals[0], holevo_vals[1])
    f0 = fidelity_with_pure(holevo_vals[0], psi0)

    # |+> state
    logical_rho = np.outer(psi_plus, psi_plus.conj())
    rho_coded = three_qubit_encode(logical_rho)
    rho_noisy = depolarizing_channel(rho_coded, p, 3)
    rho_corrected = three_qubit_syndrome(rho_noisy)
    rho_decoded = three_qubit_decode(rho_corrected)
    tr = np.trace(rho_decoded).real
    if tr > 1e-12:
        rho_decoded = rho_decoded / tr
    f_plus = fidelity_with_pure(rho_decoded, psi_plus)

    three_q_results['holevo'].append(round(hol, 6))
    three_q_results['fidelity_0'].append(round(f0, 6))
    three_q_results['fidelity_plus'].append(round(f_plus, 6))

# ============================================================
# 3. [[5,1,3]] code under depolarizing
# ============================================================
print("=== [[5,1,3]] CODE (depolarizing) ===")
five_q_results = {'holevo': [], 'fidelity_0': [], 'fidelity_plus': []}

for p in p_values:
    print(f"  p={p:.2f}...")
    holevo_vals = []
    for logical_sv in [psi0, psi1]:
        logical_rho = np.outer(logical_sv, logical_sv.conj())
        rho_coded = V5 @ logical_rho @ V5.conj().T
        rho_noisy = depolarizing_channel(rho_coded, p, 5)
        rho_corrected = five_qubit_syndrome_correct(rho_noisy, V5)
        rho_decoded = V5.conj().T @ rho_corrected @ V5
        tr = np.trace(rho_decoded).real
        if tr > 1e-12:
            rho_decoded = rho_decoded / tr
        holevo_vals.append(rho_decoded)

    hol = holevo_information(holevo_vals[0], holevo_vals[1])
    f0 = fidelity_with_pure(holevo_vals[0], psi0)

    # |+> state
    logical_rho = np.outer(psi_plus, psi_plus.conj())
    rho_coded = V5 @ logical_rho @ V5.conj().T
    rho_noisy = depolarizing_channel(rho_coded, p, 5)
    rho_corrected = five_qubit_syndrome_correct(rho_noisy, V5)
    rho_decoded = V5.conj().T @ rho_corrected @ V5
    tr = np.trace(rho_decoded).real
    if tr > 1e-12:
        rho_decoded = rho_decoded / tr
    f_plus = fidelity_with_pure(rho_decoded, psi_plus)

    five_q_results['holevo'].append(round(hol, 6))
    five_q_results['fidelity_0'].append(round(f0, 6))
    five_q_results['fidelity_plus'].append(round(f_plus, 6))

# ============================================================
# Print comparison table
# ============================================================
print(f"\n{'='*75}")
print(f"HOLEVO INFORMATION COMPARISON (depolarizing noise)")
print(f"{'='*75}")
print(f"{'p':>5} | {'Uncoded':>8} | {'3-qubit':>8} | {'[[5,1,3]]':>9} | {'Best':>9}")
print("-" * 50)
for i, p in enumerate(p_values):
    u = uncoded_results['holevo'][i]
    t = three_q_results['holevo'][i]
    f = five_q_results['holevo'][i]
    best = "5q" if f >= t and f >= u else ("3q" if t >= u else "uncoded")
    print(f"{p:5.2f} | {u:8.4f} | {t:8.4f} | {f:9.4f} | {best:>9}")

print(f"\nFIDELITY COMPARISON (|0> state)")
print("-" * 50)
print(f"{'p':>5} | {'Uncoded':>8} | {'3-qubit':>8} | {'[[5,1,3]]':>9}")
print("-" * 50)
for i, p in enumerate(p_values):
    print(f"{p:5.2f} | {uncoded_results['fidelity_0'][i]:8.4f} | "
          f"{three_q_results['fidelity_0'][i]:8.4f} | {five_q_results['fidelity_0'][i]:9.4f}")

print(f"\nFIDELITY COMPARISON (|+> state)")
print("-" * 50)
print(f"{'p':>5} | {'Uncoded':>8} | {'3-qubit':>8} | {'[[5,1,3]]':>9}")
print("-" * 50)
for i, p in enumerate(p_values):
    print(f"{p:5.2f} | {uncoded_results['fidelity_plus'][i]:8.4f} | "
          f"{three_q_results['fidelity_plus'][i]:8.4f} | {five_q_results['fidelity_plus'][i]:9.4f}")

# Find crossover points
print("\n=== CROSSOVER ANALYSIS ===")
for i in range(1, len(p_values)):
    p = p_values[i]
    u = uncoded_results['holevo'][i]
    f = five_q_results['holevo'][i]
    t = three_q_results['holevo'][i]
    if f < u and five_q_results['holevo'][i-1] >= uncoded_results['holevo'][i-1]:
        print(f"[[5,1,3]] Holevo crosses below uncoded between p={p_values[i-1]:.2f} and p={p:.2f}")
    if t < u and three_q_results['holevo'][i-1] >= uncoded_results['holevo'][i-1]:
        print(f"3-qubit Holevo crosses below uncoded between p={p_values[i-1]:.2f} and p={p:.2f}")

# Save
all_results = {
    'p_values': p_values,
    'uncoded': uncoded_results,
    'three_qubit': three_q_results,
    'five_qubit': five_q_results,
    'description': 'Holevo information for QEC codes under depolarizing noise',
}

with open('results/sprint_018b_five_qubit_holevo.json', 'w') as f:
    json.dump(all_results, f, indent=2)

print("\nResults saved to results/sprint_018b_five_qubit_holevo.json")
