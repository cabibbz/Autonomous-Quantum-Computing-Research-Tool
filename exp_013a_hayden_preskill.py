"""
Sprint 013a: Hayden-Preskill Protocol — MI Recovery vs Scrambling Depth

Setup (8 qubits total):
- Reference R (qubit 0): entangled with Alice's message A (qubit 1)
- Black hole B (qubits 1,2,3): qubit 1 is Alice's message + qubits 2,3 are BH
- Early radiation B' (qubits 4,5,6): initially entangled with B qubits 2,3,... wait

Revised labeling for clarity:
- Qubits 0-3: "inside" system (qubit 0 = Alice's message, qubits 1-3 = black hole)
- Qubits 4-6: "early radiation" B' (Bob holds, entangled with BH qubits 1-3)
- Qubit 7: "reference" R (entangled with Alice's message qubit 0)

Protocol:
1. Entangle R-A: Bell pair on qubits 7,0
2. Entangle B-B': Bell pairs on (1,4), (2,5), (3,6)
3. Apply scrambling unitary U on qubits 0-3 (A+B)
4. Measure MI(R : B' ∪ D) where D is some subset of output qubits 0-3

Track MI(R : B'∪D) as scrambling depth increases.
D = {qubit 0} (1 qubit of late radiation).
"""

import numpy as np
import json
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy, random_unitary

def von_neumann_entropy(sv, qubits_to_keep, n):
    """Compute S(subsystem) by tracing out complement."""
    qubits_to_trace = [q for q in range(n) if q not in qubits_to_keep]
    if not qubits_to_trace:
        return 0.0  # full system, pure state
    rdm = partial_trace(sv, qubits_to_trace)
    return float(entropy(rdm, base=2))

def mutual_information_subsystems(sv, subsys_a, subsys_b, n):
    """MI(A:B) = S(A) + S(B) - S(AB)"""
    sa = von_neumann_entropy(sv, subsys_a, n)
    sb = von_neumann_entropy(sv, subsys_b, n)
    sab = von_neumann_entropy(sv, subsys_a + subsys_b, n)
    return sa + sb - sab

def apply_random_layer_alltoall(qc, qubits):
    """Apply random 2-qubit gates on random pairs within specified qubits."""
    available = list(qubits)
    np.random.shuffle(available)
    for idx in range(0, len(available) - 1, 2):
        u = random_unitary(4)
        qc.unitary(u, [available[idx], available[idx + 1]])

def apply_random_layer_1d(qc, qubits, layer_idx):
    """Brick-wall layer on specified qubits."""
    offset = layer_idx % 2
    for i in range(offset, len(qubits) - 1, 2):
        u = random_unitary(4)
        qc.unitary(u, [qubits[i], qubits[i + 1]])

# Parameters
n = 8  # total qubits
n_bh = 4  # black hole + Alice's message (qubits 0-3)
max_layers = 12
n_instances = 8
np.random.seed(42)

# Qubit assignments
R = [7]        # reference
A = [0]        # Alice's message (part of BH after throwing in)
BH = [0, 1, 2, 3]  # scrambling system (A + black hole)
B_prime = [4, 5, 6]  # early radiation (Bob holds)
D1 = [0]       # 1 qubit of late radiation (Bob captures)
D2 = [0, 1]    # 2 qubits of late radiation

results = {
    'n_qubits': n,
    'n_bh': n_bh,
    'max_layers': max_layers,
    'n_instances': n_instances,
    'qubit_map': {
        'reference_R': R,
        'alice_message_A': A,
        'black_hole_system': BH,
        'early_radiation_B_prime': B_prime,
        'late_radiation_D1': D1,
        'late_radiation_D2': D2,
    },
    'layer_data': []
}

print("=== Hayden-Preskill Protocol: MI Recovery vs Scrambling Depth ===\n")
print(f"Setup: {n} qubits, BH system = qubits {BH}, B' = qubits {B_prime}, R = qubit {R[0]}")
print(f"Scrambling: all-to-all random unitaries on qubits {BH}\n")

for layer in range(max_layers + 1):
    mi_r_bprime = []       # MI(R : B') — no late radiation
    mi_r_bprime_d1 = []    # MI(R : B'∪D1) — 1 qubit late radiation
    mi_r_bprime_d2 = []    # MI(R : B'∪D2) — 2 qubits late radiation
    mi_r_d1 = []           # MI(R : D1) — late radiation only, no early
    entropy_r = []

    for inst in range(n_instances):
        np.random.seed(42 + inst * 1000)

        qc = QuantumCircuit(n)

        # Step 1: Create Bell pair R-A (qubits 7, 0)
        qc.h(7)
        qc.cx(7, 0)

        # Step 2: Create Bell pairs B-B' (qubits 1-4, 2-5, 3-6)
        for i, j in [(1, 4), (2, 5), (3, 6)]:
            qc.h(i)
            qc.cx(i, j)

        # Step 3: Apply scrambling unitary on BH system (qubits 0-3)
        for l in range(layer):
            apply_random_layer_alltoall(qc, BH)

        # Get statevector
        sv = Statevector.from_instruction(qc)

        # Compute various mutual informations
        mi_1 = mutual_information_subsystems(sv, R, B_prime, n)
        mi_2 = mutual_information_subsystems(sv, R, B_prime + D1, n)
        mi_3 = mutual_information_subsystems(sv, R, B_prime + D2, n)
        mi_4 = mutual_information_subsystems(sv, R, D1, n)
        s_r = von_neumann_entropy(sv, R, n)

        mi_r_bprime.append(mi_1)
        mi_r_bprime_d1.append(mi_2)
        mi_r_bprime_d2.append(mi_3)
        mi_r_d1.append(mi_4)
        entropy_r.append(s_r)

    layer_result = {
        'layer': layer,
        'MI_R_Bprime_mean': float(np.mean(mi_r_bprime)),
        'MI_R_Bprime_std': float(np.std(mi_r_bprime)),
        'MI_R_Bprime_D1_mean': float(np.mean(mi_r_bprime_d1)),
        'MI_R_Bprime_D1_std': float(np.std(mi_r_bprime_d1)),
        'MI_R_Bprime_D2_mean': float(np.mean(mi_r_bprime_d2)),
        'MI_R_Bprime_D2_std': float(np.std(mi_r_bprime_d2)),
        'MI_R_D1_only_mean': float(np.mean(mi_r_d1)),
        'MI_R_D1_only_std': float(np.std(mi_r_d1)),
        'S_R_mean': float(np.mean(entropy_r)),
    }
    results['layer_data'].append(layer_result)

    print(f"Layer {layer:2d}: MI(R:B')={layer_result['MI_R_Bprime_mean']:.3f}  "
          f"MI(R:B'D1)={layer_result['MI_R_Bprime_D1_mean']:.3f}  "
          f"MI(R:B'D2)={layer_result['MI_R_Bprime_D2_mean']:.3f}  "
          f"MI(R:D1)={layer_result['MI_R_D1_only_mean']:.3f}")

# Theoretical expectations
print("\n--- Theoretical limits ---")
print("Before scrambling: MI(R:B') = 0 (R entangled only with A, not B')")
print("After perfect scrambling: MI(R:B'∪D) → 2.0 (full recovery from B' + 1 qubit)")
print("MI(R:D1 only) should stay small — late radiation alone can't recover info")

# Save
with open('results/sprint_013a_hayden_preskill.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_013a_hayden_preskill.json")
