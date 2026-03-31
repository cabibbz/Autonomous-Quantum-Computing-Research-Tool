"""Sprint 014c: QEC meets scrambling — information spreading in code states.

Tests the hypothesis: good QEC codes should scramble information (spread it
uniformly), connecting our scrambling results (Sprint 012-013) to QEC.

Measures:
1. OTOC-like information spreading in code states vs archetypes
2. Mutual information between encoded qubit and each physical qubit
3. Recoverability: can you recover the logical qubit from any k physical qubits?
"""

import numpy as np
import json
from itertools import combinations

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

def kron_list(mats):
    result = mats[0]
    for m in mats[1:]:
        result = np.kron(result, m)
    return result

def partial_trace(rho, keep, n):
    keep = sorted(keep)
    remove = sorted(set(range(n)) - set(keep))
    rho_t = rho.reshape([2]*n + [2]*n)
    for i, q in enumerate(remove):
        ax_row = q - i
        ax_col = ax_row + (n - i)
        rho_t = np.trace(rho_t, axis1=ax_row, axis2=ax_col)
    k = len(keep)
    return rho_t.reshape(2**k, 2**k)

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return float(-np.sum(evals * np.log2(evals)))

def mutual_information_AB(rho_full, A, B, n):
    """MI between subsystems A and B (lists of qubit indices)."""
    rho_A = partial_trace(rho_full, A, n)
    rho_B = partial_trace(rho_full, B, n)
    rho_AB = partial_trace(rho_full, sorted(A + B), n)
    return von_neumann_entropy(rho_A) + von_neumann_entropy(rho_B) - von_neumann_entropy(rho_AB)

def pauli_string(paulis, n):
    pmap = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
    return kron_list([pmap[c] for c in paulis])


# ============================================================
# Build states
# ============================================================

def build_five_qubit_code():
    n, dim = 5, 32
    stabilizers = ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']
    state_0 = np.zeros(dim, dtype=complex)
    state_0[0] = 1.0
    for stab in stabilizers:
        S = pauli_string(stab, n)
        proj = (np.eye(dim) + S) / 2.0
        state_0 = proj @ state_0
        state_0 /= np.linalg.norm(state_0)
    X_all = pauli_string('XXXXX', n)
    state_1 = X_all @ state_0
    V = np.zeros((dim, 2), dtype=complex)
    V[:, 0] = state_0
    V[:, 1] = state_1
    return V, state_0, state_1

# Reference state for MI with "logical qubit":
# Use a 6-qubit system: qubit 0 = reference, qubits 1-5 = code
# Bell pair between reference and logical, then encode
def build_reference_code_state(V):
    """Build |Φ⟩_ref = (|0⟩_R|0⟩_L + |1⟩_R|1⟩_L)/√2 encoded in 6 qubits."""
    n_total = 6
    dim = 2**n_total
    state = np.zeros(dim, dtype=complex)
    # |0⟩_R ⊗ |0⟩_L (encoded)
    state_0L = V[:, 0]  # 5-qubit |0>_L
    state_1L = V[:, 1]  # 5-qubit |1>_L
    # Reference qubit is qubit 0, code is qubits 1-5
    for i, amp in enumerate(state_0L):
        if abs(amp) > 1e-15:
            # |0>_R ⊗ state_0L[i]  → index = 0*32 + i
            state[i] += amp / np.sqrt(2)
    for i, amp in enumerate(state_1L):
        if abs(amp) > 1e-15:
            # |1>_R ⊗ state_1L[i]  → index = 1*32 + i
            state[32 + i] += amp / np.sqrt(2)
    return state / np.linalg.norm(state)

# Same for 3-qubit code
def build_reference_3qubit():
    """Bell pair between reference and 3-qubit logical."""
    n_total = 4  # 1 ref + 3 code
    dim = 2**n_total
    state = np.zeros(dim, dtype=complex)
    # |0>_L = |000>, |1>_L = |111>
    # |0>_R|000> + |1>_R|111>
    state[0b0000] = 1.0 / np.sqrt(2)  # |0>_R |000>
    state[0b1111] = 1.0 / np.sqrt(2)  # |1>_R |111>
    return state

# GHZ-like reference state for comparison
def build_reference_ghz(n_code):
    """Bell pair between reference and GHZ-encoded logical."""
    n_total = 1 + n_code
    dim = 2**n_total
    state = np.zeros(dim, dtype=complex)
    # |0>_R|00...0> + |1>_R|11...1>
    state[0] = 1.0 / np.sqrt(2)  # |0...0>
    state[dim - 1] = 1.0 / np.sqrt(2)  # |1...1>
    return state

# Product (unencoded) reference
def build_reference_product():
    """Bell pair between reference and single physical qubit (no encoding)."""
    n_total = 2
    dim = 2**n_total
    state = np.zeros(dim, dtype=complex)
    state[0b00] = 1.0 / np.sqrt(2)
    state[0b11] = 1.0 / np.sqrt(2)
    return state


# ============================================================
# Experiment 1: MI between reference and subsets of code qubits
# ============================================================
print("=== MI between reference qubit and code qubit subsets ===")
print("(How much of the logical info can you recover from k physical qubits?)")

results = {}

V5, _, _ = build_five_qubit_code()
ref_5q = build_reference_code_state(V5)
rho_5q = np.outer(ref_5q, ref_5q.conj())
n_total_5 = 6  # ref(0) + code(1-5)

# For each subset size k=1..5, compute MI(ref : subset)
print("\n[[5,1,3]] code:")
mi_by_size_5q = {}
for k in range(1, 6):
    mis = []
    for subset in combinations(range(1, 6), k):
        mi = mutual_information_AB(rho_5q, [0], list(subset), n_total_5)
        mis.append(round(float(mi), 4))
    avg_mi = round(float(np.mean(mis)), 4)
    std_mi = round(float(np.std(mis)), 4)
    mi_by_size_5q[k] = {'values': mis, 'mean': avg_mi, 'std': std_mi}
    print(f"  k={k}: MI={avg_mi} ± {std_mi}  (values: {mis})")

results['5qubit_mi_recovery'] = mi_by_size_5q

# 3-qubit code
ref_3q = build_reference_3qubit()
rho_3q = np.outer(ref_3q, ref_3q.conj())
n_total_3 = 4

print("\n3-qubit code:")
mi_by_size_3q = {}
for k in range(1, 4):
    mis = []
    for subset in combinations(range(1, 4), k):
        mi = mutual_information_AB(rho_3q, [0], list(subset), n_total_3)
        mis.append(round(float(mi), 4))
    avg_mi = round(float(np.mean(mis)), 4)
    std_mi = round(float(np.std(mis)), 4)
    mi_by_size_3q[k] = {'values': mis, 'mean': avg_mi, 'std': std_mi}
    print(f"  k={k}: MI={avg_mi} ± {std_mi}  (values: {mis})")

results['3qubit_mi_recovery'] = mi_by_size_3q

# GHZ (same as 3-qubit code for comparison)
ref_ghz5 = build_reference_ghz(5)
rho_ghz5 = np.outer(ref_ghz5, ref_ghz5.conj())

print("\nGHZ (n=5, for comparison):")
mi_by_size_ghz = {}
for k in range(1, 6):
    mis = []
    for subset in combinations(range(1, 6), k):
        mi = mutual_information_AB(rho_ghz5, [0], list(subset), 6)
        mis.append(round(float(mi), 4))
    avg_mi = round(float(np.mean(mis)), 4)
    std_mi = round(float(np.std(mis)), 4)
    mi_by_size_ghz[k] = {'values': mis, 'mean': avg_mi, 'std': std_mi}
    print(f"  k={k}: MI={avg_mi} ± {std_mi}  (values: {mis})")

results['ghz5_mi_recovery'] = mi_by_size_ghz


# ============================================================
# Experiment 2: OTOC-like measure — does encoding spread information?
# ============================================================
print("\n=== Information spreading: operator weight after encoding ===")
print("(Apply Z to logical qubit, measure where it spreads in code)")

# For [[5,1,3]], the logical Z operator is:
# Z_L = ZZZZZ (the product of all Z)
# Actually, for the [[5,1,3]] code, Z_L and X_L are:
# X_L = XXXXX, Z_L = ZZZZZ

# Measure: for each physical qubit, what is <Z_i> for |0>_L vs |1>_L?
# If Z_L is "spread", individual Z_i should show no signal.

print("\n[[5,1,3]] code:")
_, state_0L, state_1L = build_five_qubit_code()
rho_0L = np.outer(state_0L, state_0L.conj())
rho_1L = np.outer(state_1L, state_1L.conj())

z_expectations_0 = []
z_expectations_1 = []
for q in range(5):
    Zq = kron_list([Z if i == q else I2 for i in range(5)])
    exp_0 = float(np.real(np.trace(Zq @ rho_0L)))
    exp_1 = float(np.real(np.trace(Zq @ rho_1L)))
    z_expectations_0.append(round(exp_0, 4))
    z_expectations_1.append(round(exp_1, 4))
    print(f"  qubit {q}: <Z>_0L = {exp_0:.4f}, <Z>_1L = {exp_1:.4f}")

# Also check 2-body ZZ
print("  2-body ZZ:")
zz_expectations_0 = {}
for i, j in combinations(range(5), 2):
    ZZ = kron_list([Z if k in (i, j) else I2 for k in range(5)])
    exp_0 = float(np.real(np.trace(ZZ @ rho_0L)))
    zz_expectations_0[f"{i}-{j}"] = round(exp_0, 4)
    if abs(exp_0) > 0.01:
        print(f"    <Z{i}Z{j}>_0L = {exp_0:.4f}")
print(f"  All ZZ expectations: {zz_expectations_0}")

results['5qubit_z_spread'] = {
    'Z_single_0L': z_expectations_0,
    'Z_single_1L': z_expectations_1,
    'ZZ_pairs_0L': zz_expectations_0
}

# 3-qubit code comparison
print("\n3-qubit code:")
state_0_3q = np.array([1, 0, 0, 0, 0, 0, 0, 0], dtype=complex)  # |000>
state_1_3q = np.array([0, 0, 0, 0, 0, 0, 0, 1], dtype=complex)  # |111>
rho_0_3q = np.outer(state_0_3q, state_0_3q.conj())
rho_1_3q = np.outer(state_1_3q, state_1_3q.conj())

for q in range(3):
    Zq = kron_list([Z if i == q else I2 for i in range(3)])
    exp_0 = float(np.real(np.trace(Zq @ rho_0_3q)))
    exp_1 = float(np.real(np.trace(Zq @ rho_1_3q)))
    print(f"  qubit {q}: <Z>_0L = {exp_0:.4f}, <Z>_1L = {exp_1:.4f}")


# ============================================================
# Experiment 3: Decoupling — how many qubits needed for recovery?
# ============================================================
print("\n=== Quantum Singleton bound and recovery threshold ===")
print("For [[n,k,d]] code: need n-d+1 qubits to recover k logical qubits")
print("[[5,1,3]]: need 5-3+1 = 3 qubits (majority of 5)")
print("3-qubit:   need 3-1+1 = 3 qubits (all 3 — no redundancy against general errors)")

# Verify: MI(ref : any 3 of 5) should be 2.0 (full recovery)
print("\n[[5,1,3]] verification:")
for k in [2, 3, 4]:
    mis = []
    for subset in combinations(range(1, 6), k):
        mi = mutual_information_AB(rho_5q, [0], list(subset), n_total_5)
        mis.append(round(float(mi), 4))
    print(f"  MI(ref : any {k} of 5) = {mis}")

results['singleton_verification'] = {
    '5qubit_3of5': [mutual_information_AB(rho_5q, [0], list(s), n_total_5) for s in combinations(range(1, 6), 3)],
    '5qubit_2of5': [mutual_information_AB(rho_5q, [0], list(s), n_total_5) for s in combinations(range(1, 6), 2)],
}
# Convert to serializable
results['singleton_verification'] = {k: [round(float(v), 4) for v in vals] for k, vals in results['singleton_verification'].items()}

# Summary
print("\n=== SUMMARY ===")
print("Information recovery threshold (MI=2.0 with reference):")
print("  [[5,1,3]]: k=3 qubits sufficient (ANY 3 of 5)")
print("  GHZ/3-qubit: k=n qubits needed (ALL qubits)")
print("  Product: k=1 qubit sufficient (trivially)")
print("\n[[5,1,3]] achieves Hayden-Preskill-like democratization WITHOUT scrambling:")
print("  Any 3 physical qubits contain the full logical qubit.")
print("  This is the static version of what scrambling does dynamically.")

# Save
with open('results/sprint_014c_qec_scrambling.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_014c_qec_scrambling.json")
