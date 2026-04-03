"""
Sprint 013b: Per-Qubit Recovery — The Hayden-Preskill Miracle

The key insight: BEFORE scrambling, recovering Alice's info requires
capturing her specific qubit. AFTER scrambling, ANY qubit of late
radiation is equally useful (combined with early radiation B').

For each output qubit q in {0,1,2,3}, compute MI(R : B' ∪ {q}).
Also compute cumulative MI as Bob captures 0,1,2,3,4 qubits of late radiation.

Compare: no scrambling vs 5 layers vs 10 layers of scrambling.
"""

import numpy as np
import json
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy, random_unitary

def von_neumann_entropy(sv, qubits_to_keep, n):
    qubits_to_trace = [q for q in range(n) if q not in qubits_to_keep]
    if not qubits_to_trace:
        return 0.0
    rdm = partial_trace(sv, qubits_to_trace)
    return float(entropy(rdm, base=2))

def mutual_information_subsystems(sv, subsys_a, subsys_b, n):
    sa = von_neumann_entropy(sv, subsys_a, n)
    sb = von_neumann_entropy(sv, subsys_b, n)
    sab = von_neumann_entropy(sv, subsys_a + subsys_b, n)
    return sa + sb - sab

def apply_random_layer_alltoall(qc, qubits):
    available = list(qubits)
    np.random.shuffle(available)
    for idx in range(0, len(available) - 1, 2):
        u = random_unitary(4)
        qc.unitary(u, [available[idx], available[idx + 1]])

n = 8
BH = [0, 1, 2, 3]
B_prime = [4, 5, 6]
R = [7]
n_instances = 10

scrambling_depths = [0, 1, 3, 6, 10]

results = {
    'n_qubits': n,
    'scrambling_depths': scrambling_depths,
    'n_instances': n_instances,
    'per_qubit': {},   # MI(R:B'∪{q}) for each output qubit q
    'cumulative': {},  # MI(R:B'∪D) as |D| grows
}

print("=== Per-Qubit Recovery: Which Qubit Matters? ===\n")

for depth in scrambling_depths:
    print(f"--- Scrambling depth = {depth} layers ---")

    per_qubit_mi = {q: [] for q in BH}
    cumulative_mi = {k: [] for k in range(5)}  # 0 to 4 qubits captured

    for inst in range(n_instances):
        np.random.seed(42 + inst * 1000)

        qc = QuantumCircuit(n)
        # Bell pairs
        qc.h(7); qc.cx(7, 0)
        for i, j in [(1, 4), (2, 5), (3, 6)]:
            qc.h(i); qc.cx(i, j)

        # Scrambling
        for l in range(depth):
            apply_random_layer_alltoall(qc, BH)

        sv = Statevector.from_instruction(qc)

        # Per-qubit: MI(R : B' ∪ {q}) for each output qubit
        for q in BH:
            mi = mutual_information_subsystems(sv, R, B_prime + [q], n)
            per_qubit_mi[q].append(mi)

        # Cumulative: MI(R : B' ∪ D) as |D| grows (capture qubits in order 0,1,2,3)
        for k in range(5):
            D = BH[:k]
            bob_system = B_prime + D
            mi = mutual_information_subsystems(sv, R, bob_system, n)
            cumulative_mi[k].append(mi)

    # Per-qubit results
    print(f"  Per-qubit MI(R:B'∪{{q}}):")
    per_qubit_summary = {}
    for q in BH:
        mean_mi = float(np.mean(per_qubit_mi[q]))
        std_mi = float(np.std(per_qubit_mi[q]))
        label = "(Alice's qubit)" if q == 0 else ""
        print(f"    q={q}: MI = {mean_mi:.3f} ± {std_mi:.3f}  {label}")
        per_qubit_summary[str(q)] = {'mean': mean_mi, 'std': std_mi}

    # Spread = max - min (measure of position-dependence)
    means = [per_qubit_summary[str(q)]['mean'] for q in BH]
    spread = max(means) - min(means)
    print(f"  Spread (max-min): {spread:.3f}  {'← position matters!' if spread > 0.5 else '← equalized!'}")

    # Cumulative results
    print(f"  Cumulative MI(R:B'∪D[0..k]):")
    cumulative_summary = {}
    for k in range(5):
        mean_mi = float(np.mean(cumulative_mi[k]))
        std_mi = float(np.std(cumulative_mi[k]))
        print(f"    |D|={k}: MI = {mean_mi:.3f} ± {std_mi:.3f}")
        cumulative_summary[str(k)] = {'mean': mean_mi, 'std': std_mi}

    results['per_qubit'][str(depth)] = {
        'per_qubit': per_qubit_summary,
        'spread': float(spread),
        'cumulative': cumulative_summary,
    }
    print()

# Summary
print("=== SUMMARY: Information Delocalization ===")
print(f"{'Depth':>5}  {'Spread':>7}  {'q=0 (Alice)':>12}  {'q=1':>7}  {'q=2':>7}  {'q=3':>7}")
for depth in scrambling_depths:
    d = results['per_qubit'][str(depth)]
    pq = d['per_qubit']
    print(f"{depth:5d}  {d['spread']:7.3f}  {pq['0']['mean']:12.3f}  "
          f"{pq['1']['mean']:7.3f}  {pq['2']['mean']:7.3f}  {pq['3']['mean']:7.3f}")

with open('results/sprint_013b_per_qubit_recovery.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_013b_per_qubit_recovery.json")
