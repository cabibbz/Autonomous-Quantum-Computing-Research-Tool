"""
Sprint 020b: Topological Entanglement Entropy (TEE)
Compute the Kitaev-Preskill TEE for the toric code and compare with
algebraic codes. TEE = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC.
For toric code: expect S_topo = -ln(2) ≈ -0.693 (in nats) = -1.0 (in bits).
For algebraic codes: expect S_topo = 0 (no topological order).
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

def build_stabilizer_state(n, stabs_x, stabs_z):
    """Build code ground state from X-type and Z-type stabilizers."""
    dim = 2**n
    I2 = np.eye(2)
    X = np.array([[0, 1], [1, 0]])
    Z = np.array([[1, 0], [0, -1]])

    def pauli_string(qubits, pauli, n):
        result = np.eye(1)
        for q in range(n):
            result = np.kron(result, pauli if q in qubits else I2)
        return result

    proj = np.eye(dim)
    for qubits in stabs_x:
        S = pauli_string(qubits, X, n)
        proj = proj @ (np.eye(dim) + S) / 2
    for qubits in stabs_z:
        S = pauli_string(qubits, Z, n)
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

# ============================================================
# Build code states
# ============================================================

# Toric code [[8,2,2]]
toric_vertex_stabs = [[0, 1, 4, 6], [0, 1, 5, 7], [2, 3, 4, 6]]
toric_plaquette_stabs = [[0, 2, 4, 5], [1, 3, 4, 5], [0, 2, 6, 7]]
sv_toric = build_stabilizer_state(8, toric_vertex_stabs, toric_plaquette_stabs)
rho_toric = density_matrix(sv_toric)
n_toric = 8

# [[5,1,3]] code
n_513 = 5
dim_513 = 2**n_513
I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])

def pauli_string_general(ops_dict, n):
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
    S = pauli_string_general(s, n_513)
    proj_513 = proj_513 @ (np.eye(dim_513) + S) / 2
state_0 = np.zeros(dim_513); state_0[0] = 1.0
sv_513 = proj_513 @ state_0
sv_513 = sv_513 / np.linalg.norm(sv_513)
rho_513 = density_matrix(sv_513)

# GHZ-8 for comparison (non-topological, known I3 structure)
sv_ghz8 = np.zeros(2**8)
sv_ghz8[0] = 1/np.sqrt(2)
sv_ghz8[255] = 1/np.sqrt(2)
rho_ghz8 = density_matrix(sv_ghz8)

# ============================================================
# Kitaev-Preskill TEE: S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC
# ============================================================
# The key: regions A, B, C must tile a region of the torus such that
# their union covers the boundary effects correctly.
# For the toric code on 2×2 torus with 8 qubits:
# We need to choose A, B, C as regions of the torus.
#
# Qubit layout on the torus:
# h00=q0, h01=q1 (horizontal edges, row 0)
# h10=q2, h11=q3 (horizontal edges, row 1)
# v00=q4, v01=q5 (vertical edges, row 0)
# v10=q6, v11=q7 (vertical edges, row 1)
#
# For TEE we need a tripartition of the qubits into A, B, C
# such that the boundaries create the right cancellation.
# Multiple partitions to average over for robustness.

def compute_tee(rho, A, B, C, n):
    """Compute Kitaev-Preskill TEE = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC."""
    AB = sorted(set(A) | set(B))
    BC = sorted(set(B) | set(C))
    AC = sorted(set(A) | set(C))
    ABC = sorted(set(A) | set(B) | set(C))

    S_A = von_neumann_entropy(partial_trace(rho, list(A), n))
    S_B = von_neumann_entropy(partial_trace(rho, list(B), n))
    S_C = von_neumann_entropy(partial_trace(rho, list(C), n))
    S_AB = von_neumann_entropy(partial_trace(rho, AB, n))
    S_BC = von_neumann_entropy(partial_trace(rho, BC, n))
    S_AC = von_neumann_entropy(partial_trace(rho, AC, n))
    S_ABC = von_neumann_entropy(partial_trace(rho, ABC, n))

    tee = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC
    return tee, {
        'S_A': round(S_A, 4), 'S_B': round(S_B, 4), 'S_C': round(S_C, 4),
        'S_AB': round(S_AB, 4), 'S_BC': round(S_BC, 4), 'S_AC': round(S_AC, 4),
        'S_ABC': round(S_ABC, 4), 'TEE': round(tee, 4)
    }

# ============================================================
# TEE for toric code with multiple partitions
# ============================================================
print("=== Topological Entanglement Entropy ===\n")

# Try multiple tripartitions of the 8 toric code qubits
# Partition by spatial regions on the torus
partitions_toric = [
    # Partition 1: by edge type and position
    ([0, 1], [2, 3], [4, 5, 6, 7]),  # horizontal row0 | horizontal row1 | all vertical
    ([0, 4], [1, 5], [2, 3, 6, 7]),  # left-col edges | right-col edges | rest
    ([0, 2, 4, 6], [1, 3, 5, 7], []),  # can't do empty...
    # Partition 2: geographic regions
    ([0, 4], [1, 5, 3, 7], [2, 6]),  # column 0 | column 1 | row1-horizontal+vertical
    ([0, 4, 5], [2, 6, 7], [1, 3]),  # top-left region | bottom-right | horizontal wraps
    ([0, 1, 4, 5], [2, 3, 6, 7], []),  # can't use empty
    # Partition 3: more balanced
    ([0, 2, 4], [1, 3, 5], [6, 7]),
    ([0, 5, 6], [1, 2, 7], [3, 4]),
    ([0, 3, 4, 7], [1, 2], [5, 6]),
]

# Filter out partitions with empty regions
partitions_toric = [p for p in partitions_toric if all(len(r) > 0 for r in p)]

print("Toric code [[8,2,2]] TEE with various tripartitions:")
tee_toric_results = []
for i, (A, B, C) in enumerate(partitions_toric):
    tee, details = compute_tee(rho_toric, A, B, C, n_toric)
    tee_toric_results.append({'partition': [A, B, C], 'tee': round(tee, 4), 'details': details})
    print(f"  Partition {i+1}: A={A}, B={B}, C={C}")
    print(f"    TEE = {tee:.4f}")

mean_tee_toric = np.mean([r['tee'] for r in tee_toric_results])
print(f"\n  Mean TEE: {mean_tee_toric:.4f}")
print(f"  Expected (topological): -1.0000 bits (-ln2 in nats)")

# ============================================================
# TEE for [[5,1,3]] — should be 0 (no topological order)
# ============================================================
partitions_513 = [
    ([0, 1], [2, 3], [4]),
    ([0, 2], [1, 4], [3]),
    ([0, 3], [1, 2], [4]),
    ([0, 4], [2, 3], [1]),
    ([1, 3], [0, 4], [2]),
]

print("\n[[5,1,3]] TEE with various tripartitions:")
tee_513_results = []
for i, (A, B, C) in enumerate(partitions_513):
    tee, details = compute_tee(rho_513, A, B, C, n_513)
    tee_513_results.append({'partition': [A, B, C], 'tee': round(tee, 4), 'details': details})
    print(f"  Partition {i+1}: A={A}, B={B}, C={C}")
    print(f"    TEE = {tee:.4f}")

mean_tee_513 = np.mean([r['tee'] for r in tee_513_results])
print(f"\n  Mean TEE: {mean_tee_513:.4f}")
print(f"  Expected (algebraic, no topology): 0.0000")

# ============================================================
# TEE for GHZ-8 — should be 0 (no topological order)
# ============================================================
partitions_ghz = partitions_toric  # same partitions for fair comparison

print("\nGHZ-8 TEE with various tripartitions:")
tee_ghz_results = []
for i, (A, B, C) in enumerate(partitions_ghz):
    tee, details = compute_tee(rho_ghz8, A, B, C, 8)
    tee_ghz_results.append({'partition': [A, B, C], 'tee': round(tee, 4), 'details': details})
    print(f"  Partition {i+1}: A={A}, B={B}, C={C}")
    print(f"    TEE = {tee:.4f}")

mean_tee_ghz = np.mean([r['tee'] for r in tee_ghz_results])
print(f"\n  Mean TEE: {mean_tee_ghz:.4f}")
print(f"  Expected: 0.0000")

# ============================================================
# Entropy scaling: S(A) vs |A| for all codes
# ============================================================
print("\n=== Entropy Scaling: S(|A|) ===")

def entropy_by_size(rho, n, max_size=None):
    """Compute mean entropy for each subsystem size."""
    if max_size is None:
        max_size = n // 2
    results = {}
    for size in range(1, max_size + 1):
        entropies = []
        for subset in combinations(range(n), size):
            rho_sub = partial_trace(rho, list(subset), n)
            entropies.append(von_neumann_entropy(rho_sub))
        results[size] = {
            'mean': round(np.mean(entropies), 4),
            'min': round(np.min(entropies), 4),
            'max': round(np.max(entropies), 4),
            'std': round(np.std(entropies), 4)
        }
    return results

print("\nToric code [[8,2,2]]:")
ent_scaling_toric = entropy_by_size(rho_toric, n_toric)
for size, vals in ent_scaling_toric.items():
    print(f"  |A|={size}: S={vals['mean']:.4f} ± {vals['std']:.4f} [{vals['min']:.4f}, {vals['max']:.4f}]")

print("\n[[5,1,3]]:")
ent_scaling_513 = entropy_by_size(rho_513, n_513)
for size, vals in ent_scaling_513.items():
    print(f"  |A|={size}: S={vals['mean']:.4f} ± {vals['std']:.4f} [{vals['min']:.4f}, {vals['max']:.4f}]")

print("\nGHZ-8:")
ent_scaling_ghz = entropy_by_size(rho_ghz8, 8)
for size, vals in ent_scaling_ghz.items():
    print(f"  |A|={size}: S={vals['mean']:.4f} ± {vals['std']:.4f} [{vals['min']:.4f}, {vals['max']:.4f}]")

# The TEE should show up as a constant offset in the entropy scaling:
# S(A) = α|∂A| - γ, where γ is the TEE
print("\n=== Entropy vs boundary analysis ===")
print("For the toric code, S(A) = α|∂A| - γ where γ = TEE")
print("If TEE = -1.0 bits, then S(A) has a -1.0 constant correction")
print(f"Toric |A|=1: S = {ent_scaling_toric[1]['mean']:.4f} (boundary=4 edges on torus)")
print(f"Toric |A|=2: S = {ent_scaling_toric[2]['mean']:.4f}")
print(f"Toric |A|=3: S = {ent_scaling_toric[3]['mean']:.4f}")
print(f"Toric |A|=4: S = {ent_scaling_toric[4]['mean']:.4f}")

elapsed = time.time() - start
print(f"\nTotal time: {elapsed:.1f}s")

# Save results
results = {
    'sprint': '020b',
    'description': 'Topological entanglement entropy',
    'toric_code_tee': {
        'partitions': tee_toric_results,
        'mean_tee': round(mean_tee_toric, 4),
        'expected': -1.0
    },
    'five_qubit_tee': {
        'partitions': tee_513_results,
        'mean_tee': round(mean_tee_513, 4),
        'expected': 0.0
    },
    'ghz8_tee': {
        'partitions': tee_ghz_results,
        'mean_tee': round(mean_tee_ghz, 4),
        'expected': 0.0
    },
    'entropy_scaling': {
        'toric': ent_scaling_toric,
        'five_qubit': ent_scaling_513,
        'ghz8': ent_scaling_ghz
    },
    'elapsed_seconds': round(elapsed, 1)
}

with open('results/sprint_020b_topological_entropy.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_020b_topological_entropy.json")
