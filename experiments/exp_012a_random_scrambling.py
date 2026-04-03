"""
Sprint 012a: Random Circuit Scrambling Dynamics

Track half-cut entropy, total pairwise MI, and tripartite information
as layers of random 2-qubit gates are applied.

Compare two geometries:
- 1D brick-wall (nearest-neighbor, alternating even/odd layers)
- All-to-all (random pairs from all qubits)

n=6 qubits, up to 20 layers, 10 random instances averaged.
"""

import numpy as np
import json
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy, random_unitary

def get_rdm(sv, qubits_to_keep, n):
    """Get reduced density matrix by tracing out complement."""
    qubits_to_trace = [q for q in range(n) if q not in qubits_to_keep]
    if not qubits_to_trace:
        return sv.to_operator()
    return partial_trace(sv, qubits_to_trace)

def von_neumann_entropy(sv, qubits, n):
    """Compute von Neumann entropy of subsystem."""
    rdm = get_rdm(sv, qubits, n)
    return entropy(rdm, base=2)

def mutual_information(sv, i, j, n):
    """MI(i:j) = S(i) + S(j) - S(i,j)"""
    si = von_neumann_entropy(sv, [i], n)
    sj = von_neumann_entropy(sv, [j], n)
    sij = von_neumann_entropy(sv, [i, j], n)
    return si + sj - sij

def tripartite_info(sv, i, j, k, n):
    """I3(i:j:k) = I(i:j) + I(i:k) - I(i:jk)"""
    mi_ij = mutual_information(sv, i, j, n)
    mi_ik = mutual_information(sv, i, k, n)
    # I(i:jk) = S(i) + S(jk) - S(ijk)
    si = von_neumann_entropy(sv, [i], n)
    sjk = von_neumann_entropy(sv, [j, k], n)
    sijk = von_neumann_entropy(sv, [i, j, k], n)
    mi_ijk = si + sjk - sijk
    return mi_ij + mi_ik - mi_ijk

def apply_random_layer_1d(qc, n, layer_idx):
    """Apply one brick-wall layer of random 2-qubit gates."""
    offset = layer_idx % 2  # alternate even/odd
    for i in range(offset, n - 1, 2):
        u = random_unitary(4)
        qc.unitary(u, [i, i + 1])

def apply_random_layer_alltoall(qc, n):
    """Apply n//2 random 2-qubit gates on random non-overlapping pairs."""
    available = list(range(n))
    np.random.shuffle(available)
    for idx in range(0, n - 1, 2):
        u = random_unitary(4)
        qc.unitary(u, [available[idx], available[idx + 1]])

def measure_scrambling(sv, n):
    """Compute scrambling metrics."""
    # Half-cut entropy
    half = list(range(n // 2))
    half_entropy = von_neumann_entropy(sv, half, n)

    # Page value (maximum expected entropy for random state)
    # For n qubits with half-cut: S_page ≈ n/2 * log2 - 1/(2*ln2 * 2^n)
    # Simplified: for n=6, half-cut max is 3.0 bits
    max_entropy = n / 2  # max possible

    # Total pairwise MI
    total_mi = 0.0
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            total_mi += mutual_information(sv, i, j, n)
            count += 1

    # Average I3 over all triples
    total_i3 = 0.0
    n_triples = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                total_i3 += tripartite_info(sv, i, j, k, n)
                n_triples += 1

    return {
        'half_entropy': float(half_entropy),
        'half_entropy_fraction': float(half_entropy / max_entropy),
        'total_pairwise_mi': float(total_mi),
        'avg_pairwise_mi': float(total_mi / count),
        'avg_i3': float(total_i3 / n_triples),
    }

# Parameters
n = 6
max_layers = 15
n_instances = 10
np.random.seed(42)

results = {
    'n_qubits': n,
    'max_layers': max_layers,
    'n_instances': n_instances,
    'geometries': {}
}

for geom_name in ['1D_brickwall', 'all_to_all']:
    print(f"\n=== {geom_name} ===")

    # Store per-layer averages
    layer_data = []

    for layer in range(max_layers + 1):
        entropies = []
        mis = []
        i3s = []
        ent_fracs = []

        for inst in range(n_instances):
            # Build circuit up to this layer
            qc = QuantumCircuit(n)
            np.random.seed(42 + inst * 1000)  # reproducible per instance

            for l in range(layer):
                if geom_name == '1D_brickwall':
                    apply_random_layer_1d(qc, n, l)
                else:
                    apply_random_layer_alltoall(qc, n)

            sv = Statevector.from_instruction(qc)
            metrics = measure_scrambling(sv, n)
            entropies.append(metrics['half_entropy'])
            ent_fracs.append(metrics['half_entropy_fraction'])
            mis.append(metrics['total_pairwise_mi'])
            i3s.append(metrics['avg_i3'])

        layer_result = {
            'layer': layer,
            'half_entropy_mean': float(np.mean(entropies)),
            'half_entropy_std': float(np.std(entropies)),
            'half_entropy_frac_mean': float(np.mean(ent_fracs)),
            'total_mi_mean': float(np.mean(mis)),
            'total_mi_std': float(np.std(mis)),
            'avg_i3_mean': float(np.mean(i3s)),
            'avg_i3_std': float(np.std(i3s)),
        }
        layer_data.append(layer_result)

        print(f"  Layer {layer:2d}: S={layer_result['half_entropy_mean']:.3f}±{layer_result['half_entropy_std']:.3f}  "
              f"MI={layer_result['total_mi_mean']:.3f}  I3={layer_result['avg_i3_mean']:.4f}")

    results['geometries'][geom_name] = layer_data

# Page entropy for comparison (random state average)
# For n=6, d=64, half-cut d_A=d_B=8
# S_Page = log2(d_A) - d_A/(2*ln2*d_B*d_A) ≈ 3 - 1/(2*ln2*64) ≈ 2.989
d_A = 2**(n//2)
d_B = 2**(n - n//2)
s_page = np.log2(d_A) - d_A / (2 * np.log(2) * d_A * d_B)
results['page_entropy'] = float(s_page)
print(f"\nPage entropy (random state limit): {s_page:.4f}")

# Save
with open('results/sprint_012a_random_scrambling.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_012a_random_scrambling.json")
