"""
Sprint 013c: Geometry Comparison — 1D vs All-to-All Recovery

Compare how fast scrambling geometry equalizes per-qubit recovery.
Track:
1. Spread (max-min of per-qubit MI) — measures position-dependence
2. MI(R:B'∪D2) — recovery with 2 qubits of late radiation
3. Per-qubit MI for Alice's qubit vs farthest qubit

Also test: what if Alice's qubit enters at different positions in the 1D chain?
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

def apply_random_layer_1d(qc, qubits, layer_idx):
    offset = layer_idx % 2
    for i in range(offset, len(qubits) - 1, 2):
        u = random_unitary(4)
        qc.unitary(u, [qubits[i], qubits[i + 1]])

n = 8
BH = [0, 1, 2, 3]
B_prime = [4, 5, 6]
R = [7]
max_layers = 12
n_instances = 10

results = {
    'n_qubits': n,
    'max_layers': max_layers,
    'n_instances': n_instances,
    'geometries': {}
}

print("=== Geometry Comparison: 1D vs All-to-All Recovery ===\n")

for geom in ['1D_brickwall', 'all_to_all']:
    print(f"--- {geom} ---")
    layer_data = []

    for layer in range(max_layers + 1):
        spreads = []
        mi_d2_vals = []
        per_qubit_means = {q: [] for q in BH}

        for inst in range(n_instances):
            np.random.seed(42 + inst * 1000)

            qc = QuantumCircuit(n)
            qc.h(7); qc.cx(7, 0)
            for i, j in [(1, 4), (2, 5), (3, 6)]:
                qc.h(i); qc.cx(i, j)

            for l in range(layer):
                if geom == '1D_brickwall':
                    apply_random_layer_1d(qc, BH, l)
                else:
                    apply_random_layer_alltoall(qc, BH)

            sv = Statevector.from_instruction(qc)

            # Per-qubit MI
            mis = []
            for q in BH:
                mi = mutual_information_subsystems(sv, R, B_prime + [q], n)
                per_qubit_means[q].append(mi)
                mis.append(mi)

            spreads.append(max(mis) - min(mis))

            # MI with 2 qubits of late radiation
            mi_d2 = mutual_information_subsystems(sv, R, B_prime + [0, 1], n)
            mi_d2_vals.append(mi_d2)

        layer_result = {
            'layer': layer,
            'spread_mean': float(np.mean(spreads)),
            'spread_std': float(np.std(spreads)),
            'MI_R_Bprime_D2_mean': float(np.mean(mi_d2_vals)),
            'MI_R_Bprime_D2_std': float(np.std(mi_d2_vals)),
            'per_qubit': {str(q): {'mean': float(np.mean(per_qubit_means[q])),
                                    'std': float(np.std(per_qubit_means[q]))}
                         for q in BH}
        }
        layer_data.append(layer_result)

        if layer % 2 == 0 or layer <= 3:
            print(f"  Layer {layer:2d}: spread={layer_result['spread_mean']:.3f}  "
                  f"MI(D2)={layer_result['MI_R_Bprime_D2_mean']:.3f}  "
                  f"q0={layer_result['per_qubit']['0']['mean']:.3f}  "
                  f"q3={layer_result['per_qubit']['3']['mean']:.3f}")

    results['geometries'][geom] = layer_data
    print()

# Find equalization depth (spread < 0.1)
print("=== Equalization Depth (spread < 0.1) ===")
for geom in ['1D_brickwall', 'all_to_all']:
    for entry in results['geometries'][geom]:
        if entry['spread_mean'] < 0.1:
            print(f"  {geom}: layer {entry['layer']}")
            break
    else:
        print(f"  {geom}: not reached by layer {max_layers}")

# Part 2: Alice enters at edge vs center of 1D chain
print("\n=== Alice's Entry Position (1D only, depth=10) ===")
# Remap: Alice at position 0 (edge) vs position 1 (center-ish)
# For Alice at position p: swap qubit 0 and p before scrambling
entry_positions = {
    'edge (q=0)': [0, 1, 2, 3],   # Alice at edge of chain
    'center (q=1)': [1, 0, 2, 3],  # Alice in center of chain (swap 0,1)
}

entry_results = {}
for name, qubit_order in entry_positions.items():
    spreads = []
    per_qubit = {q: [] for q in BH}

    for inst in range(n_instances):
        np.random.seed(42 + inst * 1000)

        qc = QuantumCircuit(n)
        # Bell pair: R with Alice's qubit (which is now at qubit_order[0])
        qc.h(7); qc.cx(7, qubit_order[0])
        # Bell pairs: BH with B'
        bh_qubits = [q for q in qubit_order if q != qubit_order[0]]
        for bh_q, bp_q in zip(bh_qubits, B_prime):
            qc.h(bh_q); qc.cx(bh_q, bp_q)

        # 1D scrambling using the qubit_order as chain
        for l in range(10):
            apply_random_layer_1d(qc, qubit_order, l)

        sv = Statevector.from_instruction(qc)

        mis = []
        for q in BH:
            mi = mutual_information_subsystems(sv, R, B_prime + [q], n)
            per_qubit[q].append(mi)
            mis.append(mi)
        spreads.append(max(mis) - min(mis))

    entry_results[name] = {
        'spread_mean': float(np.mean(spreads)),
        'per_qubit': {str(q): float(np.mean(per_qubit[q])) for q in BH}
    }
    print(f"  {name}: spread={entry_results[name]['spread_mean']:.3f}  "
          f"per-qubit={[f'{v:.3f}' for v in entry_results[name]['per_qubit'].values()]}")

results['entry_position'] = entry_results

with open('results/sprint_013c_geometry_recovery.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_013c_geometry_recovery.json")
