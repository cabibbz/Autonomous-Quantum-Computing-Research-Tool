"""Sprint 014a: Entanglement structure of QEC code states.

Implements the 3-qubit bit-flip code and [[5,1,3]] perfect code,
then characterizes their entanglement using MI, I3, and entanglement spectrum.
Compares to our known archetypes (GHZ, W, Cluster).
"""

import numpy as np
import json
from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator
from itertools import combinations

sim = AerSimulator(method='statevector')

def get_statevector(qc):
    qc_copy = qc.copy()
    qc_copy.save_statevector()
    result = sim.run(qc_copy).result()
    return np.array(result.get_statevector(qc_copy))

def partial_trace(rho, keep, n):
    """Trace out all qubits except those in 'keep'.
    Uses tensor reshape approach. Qubit 0 = leftmost = most significant bit.
    """
    keep = sorted(keep)
    remove = sorted(set(range(n)) - set(keep))

    # Reshape rho into tensor with 2n indices (n for row, n for col)
    rho_t = rho.reshape([2]*n + [2]*n)

    # Trace out qubits in 'remove' by contracting row and col indices
    # We need to trace axis pairs. Each removal shifts indices.
    for i, q in enumerate(remove):
        ax_row = q - i  # row axis (shifts as we remove)
        ax_col = ax_row + (n - i)  # corresponding col axis
        rho_t = np.trace(rho_t, axis1=ax_row, axis2=ax_col)

    k = len(keep)
    return rho_t.reshape(2**k, 2**k)

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return float(-np.sum(evals * np.log2(evals)))

def mutual_information(rho_full, i, j, n):
    rho_i = partial_trace(rho_full, [i], n)
    rho_j = partial_trace(rho_full, [j], n)
    rho_ij = partial_trace(rho_full, [i, j], n)
    return von_neumann_entropy(rho_i) + von_neumann_entropy(rho_j) - von_neumann_entropy(rho_ij)

def tripartite_info(rho_full, i, j, k, n):
    S_i = von_neumann_entropy(partial_trace(rho_full, [i], n))
    S_j = von_neumann_entropy(partial_trace(rho_full, [j], n))
    S_k = von_neumann_entropy(partial_trace(rho_full, [k], n))
    S_ij = von_neumann_entropy(partial_trace(rho_full, [i, j], n))
    S_ik = von_neumann_entropy(partial_trace(rho_full, [i, k], n))
    S_jk = von_neumann_entropy(partial_trace(rho_full, [j, k], n))
    S_ijk = von_neumann_entropy(partial_trace(rho_full, [i, j, k], n))
    return S_i + S_j + S_k - S_ij - S_ik - S_jk + S_ijk

def negativity(rho_full, subsys_a, n):
    """Compute negativity across bipartition A|B via partial transpose."""
    subsys_a = sorted(subsys_a)
    subsys_b = sorted(set(range(n)) - set(subsys_a))
    n_a = len(subsys_a)
    n_b = len(subsys_b)

    # Reorder qubits: A first, then B
    perm = subsys_a + subsys_b

    # Build permutation matrix for statevector
    dim = 2**n
    d_a = 2**n_a
    d_b = 2**n_b

    # Reshape rho into tensor form, permute, then partial transpose
    rho_t = rho_full.reshape([2]*n + [2]*n)

    # Permute: move qubits to order [A... B...]
    # Row indices: 0..n-1, Col indices: n..2n-1
    row_perm = perm
    col_perm = [p + n for p in perm]
    full_perm = row_perm + col_perm
    rho_t = np.transpose(rho_t, full_perm)

    # Now reshape as (d_a, d_b, d_a, d_b) and partial transpose on A
    rho_t = rho_t.reshape(d_a, d_b, d_a, d_b)
    rho_pt = rho_t.transpose(2, 1, 0, 3)  # swap A row/col indices
    rho_pt = rho_pt.reshape(dim, dim)

    evals = np.linalg.eigvalsh(rho_pt)
    return float((np.sum(np.abs(evals)) - 1.0) / 2.0)


# ============================================================
# Build QEC code states
# ============================================================

# --- 3-qubit bit-flip code ---
def three_qubit_code(logical='plus'):
    qc = QuantumCircuit(3)
    if logical == 'plus':
        qc.h(0)
    elif logical == 'one':
        qc.x(0)
    # else: |0> on qubit 0 -> |000>_L
    qc.cx(0, 1)
    qc.cx(0, 2)
    return get_statevector(qc)

# --- [[5,1,3]] perfect code ---
def five_qubit_code(logical='zero'):
    """Encode via stabilizer projection."""
    n = 5
    dim = 2**n

    I2 = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    pauli_map = {'I': I2, 'X': X, 'Z': Z}

    def pauli_string(s):
        result = pauli_map[s[0]]
        for c in s[1:]:
            result = np.kron(result, pauli_map[c])
        return result

    # Standard [[5,1,3]] stabilizers (cyclic)
    stabilizers = ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']

    # Project |00000> onto +1 eigenspace of all stabilizers
    state = np.zeros(dim, dtype=complex)
    state[0] = 1.0

    for stab in stabilizers:
        S = pauli_string(stab)
        projector = (np.eye(dim) + S) / 2.0
        state = projector @ state
        norm = np.linalg.norm(state)
        if norm < 1e-12:
            raise ValueError(f"Stabilizer {stab} annihilated the state!")
        state = state / norm

    if logical == 'one':
        X_all = pauli_string('XXXXX')
        state = X_all @ state
    elif logical == 'plus':
        state_zero = state.copy()
        X_all = pauli_string('XXXXX')
        state_one = X_all @ state_zero
        state = (state_zero + state_one) / np.sqrt(2)
        state = state / np.linalg.norm(state)

    return state


# ============================================================
# Archetype states for comparison (n=5)
# ============================================================
def ghz_state(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(1, n):
        qc.cx(0, i)
    return get_statevector(qc)

def w_state(n):
    dim = 2**n
    sv = np.zeros(dim, dtype=complex)
    for i in range(n):
        idx = 1 << (n - 1 - i)
        sv[idx] = 1.0 / np.sqrt(n)
    return sv

def cluster_1d(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        qc.cz(i, i + 1)
    return get_statevector(qc)


# ============================================================
# Sanity check: partial_trace
# ============================================================
print("=== Sanity checks ===")
# Bell state: partial trace should give I/2
qc = QuantumCircuit(2); qc.h(0); qc.cx(0, 1)
sv = get_statevector(qc)
rho = np.outer(sv, sv.conj())
rho_0 = partial_trace(rho, [0], 2)
print(f"Bell state, trace over q1: {np.diag(rho_0).real}  (expect [0.5, 0.5])")
mi_01 = mutual_information(rho, 0, 1, 2)
print(f"Bell MI(0,1) = {mi_01:.4f}  (expect 2.0)")

# ============================================================
# Analyze all states
# ============================================================
def analyze_state(name, sv, n):
    rho = np.outer(sv, sv.conj())

    # Single-qubit entropies
    sq_ent = [round(von_neumann_entropy(partial_trace(rho, [i], n)), 4) for i in range(n)]

    # Pairwise MI
    mi_pairs = {}
    for i, j in combinations(range(n), 2):
        mi_pairs[f"{i}-{j}"] = round(mutual_information(rho, i, j, n), 4)

    # I3 for all triples
    i3_triples = {}
    for i, j, k in combinations(range(n), 3):
        i3_triples[f"{i}-{j}-{k}"] = round(tripartite_info(rho, i, j, k, n), 4)

    # Negativity spectrum
    neg_spectrum = {}
    for size in range(1, n // 2 + 1):
        negs = []
        for subsys in combinations(range(n), size):
            neg = round(negativity(rho, list(subsys), n), 4)
            negs.append(neg)
        neg_spectrum[f'|A|={size}'] = negs

    total_mi = round(sum(mi_pairs.values()), 4)
    avg_i3 = round(float(np.mean(list(i3_triples.values()))), 4)

    print(f"\n  {name} (n={n}):")
    print(f"    Single-qubit S: {sq_ent}")
    print(f"    Total MI: {total_mi}, avg I3: {avg_i3}")
    print(f"    MI: {mi_pairs}")
    print(f"    I3: {i3_triples}")
    print(f"    Neg spectrum: {neg_spectrum}")

    return {
        'single_qubit_entropies': sq_ent,
        'MI': mi_pairs, 'total_MI': total_mi,
        'I3': i3_triples, 'avg_I3': avg_i3,
        'negativity_spectrum': neg_spectrum
    }

results = {}

# 3-qubit code
print("\n=== 3-qubit bit-flip code ===")
for logical in ['zero', 'plus', 'one']:
    sv = three_qubit_code(logical)
    results[f'3qubit_{logical}'] = analyze_state(f'3qubit |{logical}>_L', sv, 3)

# 5-qubit code
print("\n=== [[5,1,3]] perfect code ===")
for logical in ['zero', 'one', 'plus']:
    sv = five_qubit_code(logical)
    results[f'5qubit_{logical}'] = analyze_state(f'5qubit |{logical}>_L', sv, 5)

# Archetypes at n=5
print("\n=== Archetypes (n=5) ===")
for name, sv in [('GHZ', ghz_state(5)), ('W', w_state(5)), ('Cluster_1D', cluster_1d(5))]:
    results[f'archetype_{name}'] = analyze_state(name, sv, 5)

# Save
with open('results/sprint_014a_qec_entanglement.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_014a_qec_entanglement.json")
