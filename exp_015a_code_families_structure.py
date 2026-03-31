"""Sprint 015a: Entanglement structure of three distance-3 code families.

Builds Steane [[7,1,3]] CSS code and Shor [[9,1,3]] concatenated code,
measures MI, I3, negativity spectrum, and compares to [[5,1,3]] perfect code.
"""

import numpy as np
import json
import time
from itertools import combinations

# ============================================================
# Helper functions (same as Sprint 014)
# ============================================================

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
    subsys_a = sorted(subsys_a)
    subsys_b = sorted(set(range(n)) - set(subsys_a))
    n_a = len(subsys_a)
    n_b = len(subsys_b)
    perm = subsys_a + subsys_b
    dim = 2**n
    d_a = 2**n_a
    d_b = 2**n_b
    rho_t = rho_full.reshape([2]*n + [2]*n)
    row_perm = perm
    col_perm = [p + n for p in perm]
    full_perm = row_perm + col_perm
    rho_t = np.transpose(rho_t, full_perm)
    rho_t = rho_t.reshape(d_a, d_b, d_a, d_b)
    rho_pt = rho_t.transpose(2, 1, 0, 3)
    rho_pt = rho_pt.reshape(dim, dim)
    evals = np.linalg.eigvalsh(rho_pt)
    return float((np.sum(np.abs(evals)) - 1.0) / 2.0)

# ============================================================
# Code construction via stabilizer projection
# ============================================================

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
pauli_map = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}

def pauli_string(s):
    result = pauli_map[s[0]]
    for c in s[1:]:
        result = np.kron(result, pauli_map[c])
    return result

def encode_via_stabilizers(stabilizers, n, logical='zero'):
    """Encode by projecting |0...0> onto +1 eigenspace of all stabilizers."""
    dim = 2**n
    state = np.zeros(dim, dtype=complex)
    state[0] = 1.0  # |0...0>

    for stab in stabilizers:
        S = pauli_string(stab)
        projector = (np.eye(dim) + S) / 2.0
        state = projector @ state
        norm = np.linalg.norm(state)
        if norm < 1e-12:
            raise ValueError(f"Stabilizer {stab} annihilated the state!")
        state = state / norm

    if logical == 'one':
        # Apply logical X (product of all X for these codes)
        X_all = pauli_string('X' * n)
        state = X_all @ state
    elif logical == 'plus':
        state_zero = state.copy()
        X_all = pauli_string('X' * n)
        state_one = X_all @ state_zero
        state = (state_zero + state_one) / np.sqrt(2)
        state = state / np.linalg.norm(state)

    return state

# ============================================================
# Define the three codes
# ============================================================

# [[5,1,3]] perfect code
stab_5 = ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ']

# Steane [[7,1,3]] CSS code
# X stabilizers from Hamming parity check
# Z stabilizers from same parity check
stab_7 = [
    'IIIXXXX',   # X on qubits 3,4,5,6
    'IXXIIXX',   # X on qubits 1,2,5,6
    'XIXIXIX',   # X on qubits 0,2,4,6
    'IIIZZZZ',   # Z on qubits 3,4,5,6
    'IZZIIZZ',   # Z on qubits 1,2,5,6
    'ZIZIZIZ',   # Z on qubits 0,2,4,6
]

# Shor [[9,1,3]] concatenated code
# 6 bit-flip stabilizers (ZZ within blocks) + 2 phase-flip stabilizers (XXXXXX between blocks)
stab_9 = [
    'ZZIIIIIII',  # q0,q1
    'IZZIIIIII',  # q1,q2
    'IIIZZIIII',  # q3,q4
    'IIIIZZIII',  # q4,q5
    'IIIIIIZZI',  # q6,q7
    'IIIIIIIZZ',  # q7,q8
    'XXXXXXIII',  # block 0-1
    'IIIXXXXXX',  # block 1-2
]

# ============================================================
# Timing test: make sure 9-qubit analysis fits in time budget
# ============================================================
print("=== Timing test (9-qubit partial_trace) ===")
t0 = time.time()
sv_test = encode_via_stabilizers(stab_9, 9, 'zero')
rho_test = np.outer(sv_test, sv_test.conj())
_ = partial_trace(rho_test, [0, 1], 9)
t1 = time.time()
print(f"9-qubit encode + partial_trace: {t1-t0:.3f}s")

# Estimate full analysis time for 9 qubits
# MI: 36 pairs x ~0.003s each ≈ 0.1s
# I3: 84 triples x 7 partial_traces x ~0.003s ≈ 1.8s
# Negativity: sizes 1-4, ~sum C(9,k) for k=1..4 = 9+36+84+126 = 255 eigendecomps
t0 = time.time()
_ = negativity(rho_test, [0], 9)
t1 = time.time()
print(f"9-qubit negativity: {t1-t0:.3f}s")
print(f"Estimated total for 9-qubit analysis: ~{255*(t1-t0) + 2.0:.1f}s")

# ============================================================
# Analyze each code
# ============================================================

def analyze_code(name, sv, n):
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

    # Negativity spectrum (up to half system size)
    neg_spectrum = {}
    for size in range(1, n // 2 + 1):
        negs = []
        for subsys in combinations(range(n), size):
            neg = round(negativity(rho, list(subsys), n), 4)
            negs.append(neg)
        neg_spectrum[f'|A|={size}'] = negs

    total_mi = round(sum(mi_pairs.values()), 4)
    avg_i3 = round(float(np.mean(list(i3_triples.values()))), 4)
    mi_values = list(mi_pairs.values())
    i3_values = list(i3_triples.values())

    print(f"\n  {name} (n={n}):")
    print(f"    Single-qubit S: {sq_ent}")
    print(f"    Total MI: {total_mi}, avg I3: {avg_i3}")
    print(f"    MI range: [{min(mi_values):.4f}, {max(mi_values):.4f}]")
    print(f"    I3 range: [{min(i3_values):.4f}, {max(i3_values):.4f}]")
    for sz, negs in neg_spectrum.items():
        print(f"    Neg {sz}: range [{min(negs):.4f}, {max(negs):.4f}], mean {np.mean(negs):.4f}")

    return {
        'n': n,
        'single_qubit_entropies': sq_ent,
        'MI': mi_pairs, 'total_MI': total_mi,
        'I3': i3_triples, 'avg_I3': avg_i3,
        'negativity_spectrum': neg_spectrum,
        'mi_range': [min(mi_values), max(mi_values)],
        'i3_range': [min(i3_values), max(i3_values)],
    }

results = {}

# --- [[5,1,3]] perfect code ---
print("\n=== [[5,1,3]] Perfect Code ===")
for logical in ['zero', 'plus']:
    t0 = time.time()
    sv = encode_via_stabilizers(stab_5, 5, logical)
    results[f'five_qubit_{logical}'] = analyze_code(f'[[5,1,3]] |{logical}>_L', sv, 5)
    print(f"    Time: {time.time()-t0:.2f}s")

# --- Steane [[7,1,3]] CSS code ---
print("\n=== Steane [[7,1,3]] CSS Code ===")
for logical in ['zero', 'plus']:
    t0 = time.time()
    sv = encode_via_stabilizers(stab_7, 7, logical)
    results[f'steane_{logical}'] = analyze_code(f'Steane [[7,1,3]] |{logical}>_L', sv, 7)
    print(f"    Time: {time.time()-t0:.2f}s")

# --- Shor [[9,1,3]] concatenated code ---
print("\n=== Shor [[9,1,3]] Concatenated Code ===")
for logical in ['zero', 'plus']:
    t0 = time.time()
    sv = encode_via_stabilizers(stab_9, 9, logical)
    results[f'shor_{logical}'] = analyze_code(f'Shor [[9,1,3]] |{logical}>_L', sv, 9)
    print(f"    Time: {time.time()-t0:.2f}s")

# ============================================================
# Summary comparison table
# ============================================================
print("\n\n=== COMPARISON TABLE (|0>_L for each code) ===")
print(f"{'Property':<25} {'[[5,1,3]]':<20} {'Steane [[7,1,3]]':<20} {'Shor [[9,1,3]]':<20}")
print("-" * 85)

for code_key, label in [('five_qubit_zero', '[[5,1,3]]'), ('steane_zero', 'Steane'), ('shor_zero', 'Shor')]:
    r = results[code_key]
    sq = r['single_qubit_entropies']
    print(f"  S(single) range: {label:<10} [{min(sq):.3f}, {max(sq):.3f}]")

print()
for code_key, label in [('five_qubit_zero', '[[5,1,3]]'), ('steane_zero', 'Steane'), ('shor_zero', 'Shor')]:
    r = results[code_key]
    print(f"  MI range: {label:<10} [{r['mi_range'][0]:.4f}, {r['mi_range'][1]:.4f}], total={r['total_MI']:.4f}")

print()
for code_key, label in [('five_qubit_zero', '[[5,1,3]]'), ('steane_zero', 'Steane'), ('shor_zero', 'Shor')]:
    r = results[code_key]
    print(f"  I3 range: {label:<10} [{r['i3_range'][0]:.4f}, {r['i3_range'][1]:.4f}], avg={r['avg_I3']:.4f}")

# Save results
with open('results/sprint_015a_code_families_structure.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_015a_code_families_structure.json")
