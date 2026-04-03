"""
Sprint 020b: Topological Degeneracy and Local Indistinguishability
The hallmark of topological order: ground states are locally identical
but globally distinct. Prepare all 4 toric code ground states, verify
local indistinguishability, and show global distinguishability via
non-contractible loop operators.
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

def multi_pauli(ops_dict, n):
    """ops_dict: qubit -> Pauli matrix"""
    result = np.eye(1)
    for q in range(n):
        result = np.kron(result, ops_dict.get(q, I2))
    return result

# Stabilizer generators
vertex_stabs = [[0, 1, 4, 6], [0, 1, 5, 7], [2, 3, 4, 6]]
plaquette_stabs = [[0, 2, 4, 5], [1, 3, 4, 5], [0, 2, 6, 7]]

# Build projector onto code space
proj = np.eye(dim)
for qubits in vertex_stabs:
    S = pauli_string(qubits, X, n)
    proj = proj @ (np.eye(dim) + S) / 2
for qubits in plaquette_stabs:
    S = pauli_string(qubits, Z, n)
    proj = proj @ (np.eye(dim) + S) / 2

# The code space is 4-dimensional (2 logical qubits)
# Logical operators for the toric code on 2x2 torus:
# X_L1: X on all horizontal edges in row 0 = X on q0, q1
# X_L2: X on all vertical edges in column 0 = X on q4, q6
# Z_L1: Z on all vertical edges in column 0 = Z on q4, q6
# Z_L2: Z on all horizontal edges in row 0 = Z on q0, q1
#
# Wait, need to be more careful. Non-contractible loops:
# Horizontal loop (row 0): edges q0, q1 -> Z_L1 = Z_q0 Z_q1 (plaquette-type)
# Horizontal loop (row 1): edges q2, q3 -> equivalent
# Vertical loop (col 0): edges q4, q6 -> Z_L2 = Z_q4 Z_q6
# Vertical loop (col 1): edges q5, q7 -> equivalent
#
# X-type logical operators (dual):
# X_L1 acts on a non-contractible cut crossing the horizontal loop
# = X on vertical edges crossing the cut = X_q4 X_q6 (for vertical cut)
# X_L2 = X_q0 X_q1 (for horizontal cut)
#
# Actually the logical operators are:
# Z_L1 = Z on edges along horizontal non-contractible loop = Z_q0 Z_q1
# Z_L2 = Z on edges along vertical non-contractible loop = Z_q4 Z_q6
# X_L1 = X on edges crossing horizontal loop = X_q4 X_q6 (perpendicular)
# X_L2 = X on edges crossing vertical loop = X_q0 X_q1 (perpendicular)
#
# Hmm, this would make Z_L1 and X_L2 the same qubits, which breaks commutation.
# Let me reconsider.
#
# For the 2x2 toric code:
# Non-contractible loops go around the torus in two directions.
# Direction 1 (horizontal): q0, q1 (row 0 horizontal edges)
# Direction 2 (vertical): q4, q6 (col 0 vertical edges)
#
# Z-string around horizontal loop: Z_q0 Z_q1 -> this is Z_L1
# X-string around vertical loop perpendicular to horizontal: X_q4 X_q5
# Wait, I need X on edges that CROSS the horizontal loop.
# The horizontal loop goes through vertices (0,0)-(0,1)-(0,0) [wrapping].
# An X-string perpendicular to this crosses from row 0 to row 1.
# The edges crossing are v00=q4 and v01=q5 (or v10=q6, v11=q7).
# So X_L1 = X_q4 X_q5 (or equivalently X_q6 X_q7)
#
# Z_L2 = Z_q4 Z_q6 (vertical non-contractible loop)
# X_L2 = X_q0 X_q2 (crosses vertical loop) or equivalently X_q1 X_q3
#
# Check commutation:
# [Z_L1, X_L1] = [Z_q0 Z_q1, X_q4 X_q5] = 0 (different qubits) -> commutes? No!
# Z_L1 and X_L1 should ANTICOMMUTE.
#
# Let me reconsider. Z_L1 = Z_q0 Z_q1, X_L1 = X_q4 X_q5.
# These act on different qubits, so they commute. That's wrong.
#
# The correct pairing: logical X and Z for the SAME logical qubit should anticommute.
# Z_L1 measures the "horizontal parity" (Z-string around horizontal loop)
# X_L1 should create a "horizontal excitation" (X-string perpendicular to horizontal loop)
# For anticommutation, they must share an odd number of qubits.
#
# Let me use the standard construction:
# Logical qubit 1: associated with horizontal non-contractible cycle
#   Z_L1 = Z-string along horizontal cycle = Z_q0 Z_q1
#   X_L1 = X-string along a PATH that crosses the horizontal cycle once
#         = X on vertical edges in one column = X_q4 X_q6
#   Check: Z_q0 Z_q1 and X_q4 X_q6 share no qubits -> commutes -> WRONG
#
# I think the issue is that for the toric code, the logical X for logical qubit 1
# is a string of X operators along the DUAL lattice that intersects the Z-string once.
# But on a 2x2 torus, every non-contractible loop intersects every other non-contractible
# loop an even number of times!
#
# Actually on a 2x2 torus (L=2), every non-contractible loop has length 2.
# A horizontal loop and a vertical loop intersect at exactly 1 point on the torus,
# but the strings are on edges, not vertices...
#
# Let me just construct the logical operators explicitly and verify.
# The standard construction:
# For a torus with L=2, the lattice has 4 vertices and 8 edges.
# Logical Z1: Z on edges of horizontal cycle (row 0) = Z_q0 Z_q1
# Logical X1: X on edges dual to vertical cycle = X on vertical edges in col 0 = X_q4 X_q6
#
# [Z_q0 Z_q1, X_q4 X_q6]: these operators act on completely different qubits,
# so they commute. But logical X and Z should anticommute!
#
# The issue: on a 2x2 torus, a horizontal Z-loop and a vertical X-loop
# share ZERO edges. They cross at vertices, not edges.
#
# I think the correct pairing for the toric code is:
# Logical qubit 1:
#   Z_L1 = product of Z along one non-contractible cycle
#   X_L1 = product of X along the SAME DIRECTION cycle (on the dual lattice)
# But dual lattice edges correspond to FACES of the original lattice...
#
# Let me just use the projector approach to find the 4 ground states.

# Find the 4 ground states by diagonalizing within the code space
print("Finding the 4 toric code ground states...")

# The projector has rank 4 (code space dimension)
eigenvalues, eigenvectors = np.linalg.eigh(proj)
# Code space basis = eigenvectors with eigenvalue 1
code_basis_idx = np.where(eigenvalues > 0.5)[0]
print(f"Code space dimension: {len(code_basis_idx)}")
code_basis = eigenvectors[:, code_basis_idx]

# These are 4 orthonormal states spanning the code space
# To get the 4 "canonical" ground states, we need the logical operators.
# Let me find them by trying candidate operators and checking.

# Candidate logical Z operators: products of Z along non-contractible loops
# Horizontal loops: {q0, q1}, {q2, q3}
# Vertical loops: {q4, q6}, {q5, q7}
# These should be related by stabilizers: Z_q0 Z_q1 * B(0,0) * B(1,0) = Z_q2 Z_q3
# So Z_q0 Z_q1 and Z_q2 Z_q3 are equivalent logical operators.

Z_L1 = pauli_string([0, 1], Z, n)  # horizontal loop
Z_L2 = pauli_string([4, 6], Z, n)  # vertical loop

# Project into code space
Z_L1_code = code_basis.conj().T @ Z_L1 @ code_basis
Z_L2_code = code_basis.conj().T @ Z_L2 @ code_basis

print(f"\nZ_L1 in code space:\n{np.round(Z_L1_code.real, 3)}")
print(f"\nZ_L2 in code space:\n{np.round(Z_L2_code.real, 3)}")

# Check if they commute (they should — two independent Z logicals)
commutator = Z_L1_code @ Z_L2_code - Z_L2_code @ Z_L1_code
print(f"\n[Z_L1, Z_L2] commutator norm: {np.linalg.norm(commutator):.6f}")

# Find simultaneous eigenstates of Z_L1 and Z_L2
# Joint diagonalization
# First diagonalize Z_L1
evals1, evecs1 = np.linalg.eigh(Z_L1_code)
print(f"\nZ_L1 eigenvalues: {np.round(evals1, 4)}")

# For each eigenspace of Z_L1, diagonalize Z_L2
ground_states = []
labels = []
for eval_target in [1.0, -1.0]:
    idx = np.where(np.abs(evals1 - eval_target) < 0.1)[0]
    sub_evecs = evecs1[:, idx]
    Z_L2_sub = sub_evecs.conj().T @ Z_L2_code @ sub_evecs
    evals2, evecs2 = np.linalg.eigh(Z_L2_sub)
    for j, ev2 in enumerate(evals2):
        vec_code = sub_evecs @ evecs2[:, j]
        vec_full = code_basis @ vec_code
        ground_states.append(vec_full)
        labels.append(f"|{'+' if eval_target > 0 else '-'}{'+' if ev2 > 0 else '-'}>")
        print(f"  Ground state {labels[-1]}: Z_L1={eval_target:+.0f}, Z_L2={ev2:+.1f}")

# ============================================================
# Local indistinguishability: reduced density matrices
# ============================================================
print("\n=== Local Indistinguishability ===")
print("Compare single-qubit and 2-qubit RDMs across all 4 ground states.\n")

# Single-qubit RDMs
print("Single-qubit trace distances between ground states:")
max_single_dist = 0
for q in range(n):
    dists = []
    for i in range(4):
        rho_i = density_matrix(ground_states[i])
        rdm_i = partial_trace(rho_i, [q], n)
        for j in range(i+1, 4):
            rho_j = density_matrix(ground_states[j])
            rdm_j = partial_trace(rho_j, [q], n)
            # Trace distance
            diff = rdm_i - rdm_j
            evals = np.linalg.eigvalsh(diff)
            td = np.sum(np.abs(evals)) / 2
            dists.append(td)
            max_single_dist = max(max_single_dist, td)
    print(f"  q{q}: max trace distance = {max(dists):.6f}")

print(f"\nMax single-qubit trace distance across all states: {max_single_dist:.6f}")
if max_single_dist < 1e-6:
    print("  => Ground states are PERFECTLY indistinguishable at single-qubit level!")

# Two-qubit RDMs
print("\nTwo-qubit trace distances between ground states:")
max_two_dist = 0
two_qubit_dists = {}
for q1, q2 in combinations(range(n), 2):
    dists = []
    for i in range(4):
        rho_i = density_matrix(ground_states[i])
        rdm_i = partial_trace(rho_i, [q1, q2], n)
        for j in range(i+1, 4):
            rho_j = density_matrix(ground_states[j])
            rdm_j = partial_trace(rho_j, [q1, q2], n)
            diff = rdm_i - rdm_j
            evals = np.linalg.eigvalsh(diff)
            td = np.sum(np.abs(evals)) / 2
            dists.append(td)
            max_two_dist = max(max_two_dist, td)
    two_qubit_dists[f"{q1}-{q2}"] = round(max(dists), 6)

# Only print nonzero
nonzero_2q = {k: v for k, v in two_qubit_dists.items() if v > 1e-6}
print(f"  Nonzero 2-qubit distinguishability: {len(nonzero_2q)}/{len(two_qubit_dists)} pairs")
for k, v in sorted(nonzero_2q.items(), key=lambda x: -x[1])[:10]:
    print(f"    ({k}): trace distance = {v:.6f}")
if len(nonzero_2q) == 0:
    print("  => Ground states are PERFECTLY indistinguishable at 2-qubit level!")

# Three-qubit RDMs (sample a few)
print("\nThree-qubit trace distances (sample):")
max_three_dist = 0
three_qubit_dists = {}
for q1, q2, q3 in list(combinations(range(n), 3))[:20]:
    dists = []
    for i in range(4):
        rho_i = density_matrix(ground_states[i])
        rdm_i = partial_trace(rho_i, [q1, q2, q3], n)
        for j in range(i+1, 4):
            rho_j = density_matrix(ground_states[j])
            rdm_j = partial_trace(rho_j, [q1, q2, q3], n)
            diff = rdm_i - rdm_j
            evals = np.linalg.eigvalsh(diff)
            td = np.sum(np.abs(evals)) / 2
            dists.append(td)
            max_three_dist = max(max_three_dist, td)
    three_qubit_dists[f"{q1}-{q2}-{q3}"] = round(max(dists), 6)

nonzero_3q = {k: v for k, v in three_qubit_dists.items() if v > 1e-6}
print(f"  Nonzero 3-qubit distinguishability: {len(nonzero_3q)}/{len(three_qubit_dists)} triples")
for k, v in sorted(nonzero_3q.items(), key=lambda x: -x[1])[:10]:
    print(f"    ({k}): trace distance = {v:.6f}")

# ============================================================
# Global distinguishability: non-local operators needed
# ============================================================
print("\n=== Global Distinguishability ===")
print("The ground states ARE globally distinct (orthogonal):")
for i in range(4):
    for j in range(i+1, 4):
        overlap = abs(np.dot(ground_states[i].conj(), ground_states[j]))
        print(f"  <{labels[i]}|{labels[j]}> = {overlap:.6f}")

# What subset size is needed to distinguish?
print("\nMinimum subset size to distinguish ground states:")
for size in range(1, n // 2 + 1):
    max_dist = 0
    distinguishing_subsets = 0
    total_subsets = 0
    for subset in combinations(range(n), size):
        total_subsets += 1
        for i in range(4):
            rho_i = density_matrix(ground_states[i])
            rdm_i = partial_trace(rho_i, list(subset), n)
            for j in range(i+1, 4):
                rho_j = density_matrix(ground_states[j])
                rdm_j = partial_trace(rho_j, list(subset), n)
                diff = rdm_i - rdm_j
                evals = np.linalg.eigvalsh(diff)
                td = np.sum(np.abs(evals)) / 2
                if td > 1e-6:
                    distinguishing_subsets += 1
                    max_dist = max(max_dist, td)
    print(f"  |A|={size}: {distinguishing_subsets} distinguishing subsets out of {total_subsets * 6} pairs, max distance = {max_dist:.4f}")
    if distinguishing_subsets > 0 and size <= 2:
        break  # Found minimum distinguishing size, continue to show full picture
    if size >= 4:
        break  # Don't go too far

# Full scan for all sizes
print("\nFull local indistinguishability analysis:")
for size in range(1, 5):
    any_distinguishable = False
    for subset in combinations(range(n), size):
        for i in range(4):
            rho_i = density_matrix(ground_states[i])
            rdm_i = partial_trace(rho_i, list(subset), n)
            for j in range(i+1, 4):
                rho_j = density_matrix(ground_states[j])
                rdm_j = partial_trace(rho_j, list(subset), n)
                diff = rdm_i - rdm_j
                evals = np.linalg.eigvalsh(diff)
                td = np.sum(np.abs(evals)) / 2
                if td > 1e-6:
                    any_distinguishable = True
                    break
            if any_distinguishable:
                break
        if any_distinguishable:
            break
    status = "DISTINGUISHABLE" if any_distinguishable else "INDISTINGUISHABLE"
    print(f"  |A|={size}: {status}")

# ============================================================
# Entropy of each ground state (should be identical)
# ============================================================
print("\n=== Entropy Comparison Across Ground States ===")
for gs_idx, (gs, label) in enumerate(zip(ground_states, labels)):
    rho = density_matrix(gs)
    half = list(range(n // 2))
    s_half = von_neumann_entropy(partial_trace(rho, half, n))
    sq_entropies = [von_neumann_entropy(partial_trace(rho, [q], n)) for q in range(n)]
    print(f"  {label}: half-cut S = {s_half:.4f}, single-qubit S = {[round(s, 3) for s in sq_entropies]}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save results
results = {
    'sprint': '020b',
    'description': 'Topological degeneracy and local indistinguishability',
    'ground_states': labels,
    'logical_operators': {
        'Z_L1': 'Z_q0 Z_q1 (horizontal loop)',
        'Z_L2': 'Z_q4 Z_q6 (vertical loop)',
    },
    'local_indistinguishability': {
        'max_single_qubit_trace_distance': round(max_single_dist, 6),
        'two_qubit_trace_distances': two_qubit_dists,
        'nonzero_2q_pairs': len(nonzero_2q),
        'three_qubit_trace_distances': three_qubit_dists,
        'nonzero_3q_triples': len(nonzero_3q),
    },
    'elapsed_seconds': round(elapsed, 1)
}

with open('results/sprint_020b_topo_degeneracy.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_020b_topo_degeneracy.json")
