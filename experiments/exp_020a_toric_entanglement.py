"""
Sprint 020a: Toric code entanglement structure
Prepare [[8,2,2]] toric code on 2x2 torus, compute MI, I3, negativity.
Compare with algebraic codes ([[5,1,3]], Steane, Shor) from Sprint 015.
"""

import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator
from itertools import combinations
import json, time

start = time.time()

def get_statevector(qc):
    sim = AerSimulator(method='statevector')
    qc_copy = qc.copy()
    qc_copy.save_statevector()
    result = sim.run(qc_copy, shots=1).result()
    return np.array(result.get_statevector(qc_copy))

def density_matrix(sv):
    sv = sv.reshape(-1, 1)
    return sv @ sv.conj().T

def partial_trace(rho, keep, n):
    """Trace out all qubits except those in 'keep'."""
    d = 2
    rho = rho.reshape([d] * (2 * n))
    trace_out = sorted(set(range(n)) - set(keep))
    for i, q in enumerate(sorted(trace_out)):
        # Adjust index for already-traced qubits
        q_adj = q - i
        n_cur = n - i
        # Trace: contract axis q_adj with axis q_adj + n_cur
        rho = np.trace(rho, axis1=q_adj, axis2=q_adj + n_cur)
    return rho.reshape(d**len(keep), d**len(keep))

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return -np.sum(evals * np.log2(evals))

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
    """Compute negativity across bipartition A|rest."""
    d = 2
    n_a = len(subsys_a)
    n_b = n - n_a
    # Reorder: A qubits first, then B qubits
    subsys_b = sorted(set(range(n)) - set(subsys_a))
    order = list(subsys_a) + subsys_b
    rho = rho_full.reshape([d] * (2 * n))
    # Transpose to get A|B ordering
    perm = order + [o + n for o in order]
    rho = rho.transpose(perm).reshape(d**n_a * d**n_b, d**n_a * d**n_b)
    # Partial transpose over A
    rho_pt = rho.reshape(d**n_a, d**n_b, d**n_a, d**n_b)
    rho_pt = rho_pt.transpose(2, 1, 0, 3).reshape(d**n_a * d**n_b, d**n_a * d**n_b)
    evals = np.linalg.eigvalsh(rho_pt)
    return (np.sum(np.abs(evals)) - 1) / 2

# ============================================================
# Toric code on 2x2 torus: 8 qubits (edges)
# ============================================================
# Layout: 2x2 grid with periodic boundary conditions
# Vertices: (0,0), (0,1), (1,0), (1,1)
# Edges (qubits):
#   Horizontal: q0=(0,0)-(0,1), q1=(0,1)-(0,0)[wrap], q2=(1,0)-(1,1), q3=(1,1)-(1,0)[wrap]
#   Vertical:   q4=(0,0)-(1,0), q5=(0,1)-(1,1), q6=(1,0)-(0,0)[wrap], q7=(1,1)-(0,1)[wrap]
#
# Wait - for a 2x2 torus, horizontal edges per row = 2, vertical edges per column = 2
# Total = 2*2 + 2*2 = 8 edges = 8 qubits
#
# Vertex stabilizers (X on all edges touching vertex):
#   A_v(0,0): X on q0, q1(wrap), q4, q6(wrap) -> but q1 = h(0,1)-(0,0) wraps back
#
# Let me use a cleaner labeling:
# Horizontal edges: h[row][col] = edge from (row,col) to (row, col+1 mod 2)
#   h00 = q0, h01 = q1, h10 = q2, h11 = q3
# Vertical edges: v[row][col] = edge from (row,col) to (row+1 mod 2, col)
#   v00 = q4, v01 = q5, v10 = q6, v11 = q7
#
# Vertex (r,c) stabilizer: X on h[r][c], h[r][(c-1)%2], v[r][c], v[(r-1)%2][c]
#   A(0,0): X on h00, h01, v00, v10 -> q0, q1, q4, q6
#   A(0,1): X on h01, h00, v01, v11 -> q1, q0, q5, q7
#   A(1,0): X on h10, h11, v10, v00 -> q2, q3, q6, q4
#   A(1,1): X on h11, h10, v11, v01 -> q3, q2, q7, q5
#
# Plaquette (r,c) stabilizer: Z on top, bottom, left, right edges of plaquette
#   Plaquette (r,c): h[r][c], h[(r+1)%2][c], v[r][c], v[r][(c+1)%2]
#   B(0,0): Z on h00, h10, v00, v01 -> q0, q2, q4, q5
#   B(0,1): Z on h01, h11, v01, v00 -> q1, q3, q5, q4
#   B(1,0): Z on h10, h00, v10, v11 -> q2, q0, q6, q7
#   B(1,1): Z on h11, h01, v11, v10 -> q3, q1, q7, q6
#
# Note: Product of all vertex stabs = I, product of all plaquette stabs = I
# So 3 independent vertex + 3 independent plaquette = 6 stabilizers
# 8 qubits - 6 stabilizers = 2 logical qubits

# Stabilizer generators (choosing 3 vertex + 3 plaquette):
vertex_stabs = [
    [0, 1, 4, 6],  # A(0,0)
    [1, 0, 5, 7],  # A(0,1) -- note: same as {0,1,5,7}
    [2, 3, 6, 4],  # A(1,0) -- same as {2,3,4,6}
]
plaquette_stabs = [
    [0, 2, 4, 5],  # B(0,0)
    [1, 3, 5, 4],  # B(0,1) -- same as {1,3,4,5}
    [2, 0, 6, 7],  # B(1,0) -- same as {0,2,6,7}
]

# To prepare the toric code ground state (|00>_L), we:
# 1. Start with |+>^8 (all qubits in + state)
# 2. Apply Z-stabilizer projections (plaquettes)
# OR use stabilizer-based preparation:
# Start from |0>^8, apply H to all, then project onto +1 eigenspace of all stabilizers

# More efficient: direct circuit construction
# The toric code ground state is the equal superposition of all states in the
# common +1 eigenspace of all stabilizers.
#
# Method: Start with |00000000>, apply stabilizer encoding circuit.
# For stabilizer codes, we can use the standard encoding procedure:
# 1. Prepare ancilla in |+> for each stabilizer generator
# 2. Apply controlled-stabilizer operations
# 3. This projects into the code space
#
# Simpler approach for small codes: build the projector and apply it.

def build_toric_ground_state():
    """Build toric code ground state by projecting onto stabilizer eigenspace."""
    n = 8
    dim = 2**n

    # Build Pauli matrices
    I2 = np.eye(2)
    X = np.array([[0, 1], [1, 0]])
    Z = np.array([[1, 0], [0, -1]])

    def pauli_string(ops_dict, n):
        """Build n-qubit Pauli operator. ops_dict maps qubit -> matrix."""
        result = np.eye(1)
        for q in range(n):
            result = np.kron(result, ops_dict.get(q, I2))
        return result

    # Build projector P = prod_i (I + S_i)/2
    proj = np.eye(dim)

    # Vertex stabilizers (X-type)
    for qubits in vertex_stabs:
        S = pauli_string({q: X for q in qubits}, n)
        proj = proj @ (np.eye(dim) + S) / 2

    # Plaquette stabilizers (Z-type)
    for qubits in plaquette_stabs:
        S = pauli_string({q: Z for q in qubits}, n)
        proj = proj @ (np.eye(dim) + S) / 2

    # The projection gives us the code space (4-dimensional for 2 logical qubits)
    # Get the ground state: project |0...0> and normalize
    state_0 = np.zeros(dim)
    state_0[0] = 1.0

    ground = proj @ state_0
    norm = np.linalg.norm(ground)
    if norm < 1e-10:
        # Try different initial state
        state_plus = np.ones(dim) / np.sqrt(dim)
        ground = proj @ state_plus
        norm = np.linalg.norm(ground)

    ground = ground / norm
    return ground

# Build ground state
print("Building toric code ground state...")
sv_toric = build_toric_ground_state()
rho_toric = density_matrix(sv_toric)
n_toric = 8

print(f"State norm: {np.linalg.norm(sv_toric):.6f}")
print(f"State purity: {np.real(np.trace(rho_toric @ rho_toric)):.6f}")

# ============================================================
# Also prepare [[5,1,3]] for comparison
# ============================================================
def build_five_qubit_code():
    """Build [[5,1,3]] code ground state (logical |0>)."""
    n = 5
    dim = 2**n
    I2 = np.eye(2)
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]])

    def pauli_string(ops_list, n):
        result = np.eye(1)
        for q in range(n):
            result = np.kron(result, ops_list[q] if q in ops_list else I2)
        return result

    # [[5,1,3]] stabilizers: XZZXI, IXZZX, XIXZZ, ZXIXZ
    stabs = [
        {0: X, 1: Z, 2: Z, 3: X},
        {1: X, 2: Z, 3: Z, 4: X},
        {0: X, 2: X, 3: Z, 4: Z},
        {0: Z, 1: X, 3: X, 4: Z},
    ]

    proj = np.eye(dim)
    for s in stabs:
        S = pauli_string(s, n)
        proj = proj @ (np.eye(dim) + S) / 2

    state_0 = np.zeros(dim)
    state_0[0] = 1.0
    ground = proj @ state_0
    norm = np.linalg.norm(ground)
    if norm < 1e-10:
        state_0 = np.ones(dim) / np.sqrt(dim)
        ground = proj @ state_0
        norm = np.linalg.norm(ground)
    return ground / norm

print("\nBuilding [[5,1,3]] ground state for comparison...")
sv_513 = build_five_qubit_code()
rho_513 = density_matrix(sv_513)
n_513 = 5

# ============================================================
# Compute MI matrices
# ============================================================
print("\n=== Mutual Information ===")

# Toric code MI
print("\nToric code [[8,2,2]] pairwise MI:")
mi_toric = {}
for i, j in combinations(range(n_toric), 2):
    mi = mutual_information(rho_toric, i, j, n_toric)
    mi_toric[f"{i}-{j}"] = round(mi, 4)
    if abs(mi) > 0.01:
        print(f"  MI({i},{j}) = {mi:.4f}")

total_mi_toric = sum(mi_toric.values())
nonzero_mi_toric = sum(1 for v in mi_toric.values() if abs(v) > 0.01)
print(f"  Total pairwise MI: {total_mi_toric:.4f}")
print(f"  Nonzero pairs: {nonzero_mi_toric}/{len(mi_toric)}")

# [[5,1,3]] MI
print("\n[[5,1,3]] pairwise MI:")
mi_513 = {}
for i, j in combinations(range(n_513), 2):
    mi = mutual_information(rho_513, i, j, n_513)
    mi_513[f"{i}-{j}"] = round(mi, 4)
    if abs(mi) > 0.01:
        print(f"  MI({i},{j}) = {mi:.4f}")

total_mi_513 = sum(mi_513.values())
print(f"  Total pairwise MI: {total_mi_513:.4f}")

# ============================================================
# Compute I3 for all triples
# ============================================================
print("\n=== Tripartite Information ===")

print("\nToric code I3:")
i3_toric = {}
neg_count_toric = 0
pos_count_toric = 0
for i, j, k in combinations(range(n_toric), 3):
    i3 = tripartite_info(rho_toric, i, j, k, n_toric)
    i3_toric[f"{i}-{j}-{k}"] = round(i3, 4)
    if i3 < -0.01:
        neg_count_toric += 1
    elif i3 > 0.01:
        pos_count_toric += 1

i3_values_toric = list(i3_toric.values())
print(f"  Negative I3 triples: {neg_count_toric}/{len(i3_toric)}")
print(f"  Positive I3 triples: {pos_count_toric}/{len(i3_toric)}")
print(f"  I3 range: [{min(i3_values_toric):.4f}, {max(i3_values_toric):.4f}]")
print(f"  Mean I3: {np.mean(i3_values_toric):.4f}")

# Show a few examples
print("  Sample negative I3 triples:")
for key, val in sorted(i3_toric.items(), key=lambda x: x[1])[:5]:
    print(f"    I3({key}) = {val}")

print("\n[[5,1,3]] I3:")
i3_513 = {}
neg_count_513 = 0
for i, j, k in combinations(range(n_513), 3):
    i3 = tripartite_info(rho_513, i, j, k, n_513)
    i3_513[f"{i}-{j}-{k}"] = round(i3, 4)
    if i3 < -0.01:
        neg_count_513 += 1

print(f"  Negative I3 triples: {neg_count_513}/{len(i3_513)}")
print(f"  I3 range: [{min(i3_513.values()):.4f}, {max(i3_513.values()):.4f}]")

# ============================================================
# Compute negativity spectrum
# ============================================================
print("\n=== Negativity Spectrum ===")

print("\nToric code negativity (all bipartition sizes):")
neg_toric = {}
for size in range(1, n_toric // 2 + 1):
    negs = []
    for subset in combinations(range(n_toric), size):
        neg = negativity(rho_toric, list(subset), n_toric)
        negs.append(round(neg, 4))
    neg_toric[f"|A|={size}"] = {
        'mean': round(np.mean(negs), 4),
        'min': round(np.min(negs), 4),
        'max': round(np.max(negs), 4),
        'std': round(np.std(negs), 4),
        'count': len(negs)
    }
    print(f"  |A|={size}: mean={np.mean(negs):.4f}, range=[{np.min(negs):.4f}, {np.max(negs):.4f}], std={np.std(negs):.4f} (n={len(negs)})")

print("\n[[5,1,3]] negativity:")
neg_513 = {}
for size in range(1, n_513 // 2 + 1):
    negs = []
    for subset in combinations(range(n_513), size):
        neg = negativity(rho_513, list(subset), n_513)
        negs.append(round(neg, 4))
    neg_513[f"|A|={size}"] = {
        'mean': round(np.mean(negs), 4),
        'min': round(np.min(negs), 4),
        'max': round(np.max(negs), 4),
        'std': round(np.std(negs), 4)
    }
    print(f"  |A|={size}: mean={np.mean(negs):.4f}, range=[{np.min(negs):.4f}, {np.max(negs):.4f}], std={np.std(negs):.4f}")

# ============================================================
# Single-qubit entropies
# ============================================================
print("\n=== Single-Qubit Entropies ===")
print("Toric code:")
sq_ent_toric = []
for q in range(n_toric):
    rho_q = partial_trace(rho_toric, [q], n_toric)
    s = von_neumann_entropy(rho_q)
    sq_ent_toric.append(round(s, 4))
    print(f"  S(q{q}) = {s:.4f}")

print("\n[[5,1,3]]:")
sq_ent_513 = []
for q in range(n_513):
    rho_q = partial_trace(rho_513, [q], n_513)
    s = von_neumann_entropy(rho_q)
    sq_ent_513.append(round(s, 4))
    print(f"  S(q{q}) = {s:.4f}")

# Half-cut entropy
print("\n=== Half-Cut Entropy ===")
half_toric = list(range(n_toric // 2))
rho_half = partial_trace(rho_toric, half_toric, n_toric)
s_half_toric = von_neumann_entropy(rho_half)
print(f"Toric [[8,2,2]] half-cut entropy: {s_half_toric:.4f}")

half_513 = list(range(n_513 // 2 + 1))  # 3 qubits
rho_half_513 = partial_trace(rho_513, half_513, n_513)
s_half_513 = von_neumann_entropy(rho_half_513)
print(f"[[5,1,3]] half-cut entropy: {s_half_513:.4f}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save results
results = {
    'sprint': '020a',
    'description': 'Toric code entanglement structure',
    'toric_code': {
        'n': n_toric, 'k': 2, 'd': 2,
        'pairwise_MI': mi_toric,
        'total_pairwise_MI': round(total_mi_toric, 4),
        'nonzero_MI_pairs': nonzero_mi_toric,
        'I3': i3_toric,
        'negative_I3_count': neg_count_toric,
        'positive_I3_count': pos_count_toric,
        'I3_range': [min(i3_values_toric), max(i3_values_toric)],
        'I3_mean': round(np.mean(i3_values_toric), 4),
        'negativity_spectrum': neg_toric,
        'single_qubit_entropy': sq_ent_toric,
        'half_cut_entropy': round(s_half_toric, 4)
    },
    'five_qubit_code': {
        'n': n_513, 'k': 1, 'd': 3,
        'pairwise_MI': mi_513,
        'total_pairwise_MI': round(total_mi_513, 4),
        'I3': i3_513,
        'negative_I3_count': neg_count_513,
        'negativity_spectrum': neg_513,
        'single_qubit_entropy': sq_ent_513,
        'half_cut_entropy': round(s_half_513, 4)
    },
    'elapsed_seconds': round(elapsed, 1)
}

with open('results/sprint_020a_toric_entanglement.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_020a_toric_entanglement.json")
