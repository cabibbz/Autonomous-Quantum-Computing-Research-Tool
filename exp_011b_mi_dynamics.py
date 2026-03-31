"""Sprint 011b: Mutual information dynamics — track pairwise MI at each gate layer.

Instead of the coarse half-cut, we track ALL pairwise mutual information values
as gates are applied. This reveals how the correlation structure builds up —
does it appear gradually or suddenly?

n=6 qubits throughout.
"""

import numpy as np
import json
from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.quantum_info import DensityMatrix, partial_trace, entropy

n = 6
sim = AerSimulator(method='statevector')

def get_dm(qc):
    qc_copy = qc.copy()
    qc_copy.save_density_matrix()
    result = sim.run(qc_copy).result()
    return DensityMatrix(result.data()['density_matrix'])

def pairwise_mi(dm, n):
    """Compute mutual information I(i:j) for all pairs."""
    # Single-qubit entropies
    s1 = {}
    for i in range(n):
        trace_out = [j for j in range(n) if j != i]
        rho_i = partial_trace(dm, trace_out)
        s1[i] = float(entropy(rho_i, base=2))

    # Two-qubit entropies and MI
    mi_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            trace_out = [k for k in range(n) if k != i and k != j]
            rho_ij = partial_trace(dm, trace_out)
            s_ij = float(entropy(rho_ij, base=2))
            mi = s1[i] + s1[j] - s_ij
            mi_matrix[i][j] = mi
            mi_matrix[j][i] = mi

    return mi_matrix, s1

def total_mi(mi_matrix, n):
    """Sum of all pairwise MI."""
    return sum(mi_matrix[i][j] for i in range(n) for j in range(i+1, n))

def mi_sparsity(mi_matrix, n, threshold=0.01):
    """Fraction of pairs with MI > threshold."""
    total_pairs = n * (n - 1) // 2
    nonzero = sum(1 for i in range(n) for j in range(i+1, n) if mi_matrix[i][j] > threshold)
    return nonzero / total_pairs

# Reuse layer functions from 011a
def ghz_layers(n):
    circuits, labels = [], []
    qc = QuantumCircuit(n)
    circuits.append(qc.copy()); labels.append("|0>^n")
    qc.h(0)
    circuits.append(qc.copy()); labels.append("H(0)")
    for i in range(1, n):
        qc.cx(0, i)
        circuits.append(qc.copy()); labels.append(f"CX(0,{i})")
    return circuits, labels

def w_layers(n):
    circuits, labels = [], []
    qc = QuantumCircuit(n)
    circuits.append(qc.copy()); labels.append("|0>^n")
    qc.x(0)
    circuits.append(qc.copy()); labels.append("X(0)")
    for k in range(n - 1):
        theta = 2 * np.arccos(np.sqrt(1.0 / (n - k)))
        qc.cry(theta, k, k + 1)
        qc.cx(k + 1, k)
        circuits.append(qc.copy()); labels.append(f"dist→q{k+1}")
    return circuits, labels

def cluster1d_layers(n):
    circuits, labels = [], []
    qc = QuantumCircuit(n)
    circuits.append(qc.copy()); labels.append("|0>^n")
    for i in range(n):
        qc.h(i)
    circuits.append(qc.copy()); labels.append("H⊗n")
    for i in range(n - 1):
        qc.cz(i, i + 1)
        circuits.append(qc.copy()); labels.append(f"CZ({i},{i+1})")
    return circuits, labels

def cluster2d_layers(n):
    circuits, labels = [], []
    qc = QuantumCircuit(n)
    circuits.append(qc.copy()); labels.append("|0>^n")
    for i in range(n):
        qc.h(i)
    circuits.append(qc.copy()); labels.append("H⊗n")
    edges = [(0,1), (1,2), (3,4), (4,5), (0,3), (1,4), (2,5)]
    for (a, b) in edges:
        qc.cz(a, b)
        circuits.append(qc.copy()); labels.append(f"CZ({a},{b})")
    return circuits, labels


results = {}

for name, layer_fn in [("GHZ", ghz_layers), ("W", w_layers),
                         ("Cluster_1D", cluster1d_layers), ("Cluster_2D", cluster2d_layers)]:
    print(f"\n=== {name} ===")
    circuits, labels = layer_fn(n)

    layer_data = []

    for idx, (qc, label) in enumerate(zip(circuits, labels)):
        dm = get_dm(qc)
        mi_mat, s1_dict = pairwise_mi(dm, n)

        t_mi = total_mi(mi_mat, n)
        sparsity = mi_sparsity(mi_mat, n)
        max_mi = float(np.max(mi_mat))

        # Which pairs have significant MI?
        sig_pairs = [(i, j, round(mi_mat[i][j], 4))
                     for i in range(n) for j in range(i+1, n) if mi_mat[i][j] > 0.01]

        layer_data.append({
            "layer": idx,
            "label": label,
            "total_mi": round(t_mi, 4),
            "sparsity": round(sparsity, 4),
            "max_mi": round(max_mi, 4),
            "significant_pairs": sig_pairs,
            "single_qubit_entropies": {str(k): round(v, 4) for k, v in s1_dict.items()},
            "mi_matrix": [[round(mi_mat[i][j], 4) for j in range(n)] for i in range(n)]
        })

        print(f"  L{idx} [{label}]: total_MI={t_mi:.3f}, sparsity={sparsity:.2f}, max={max_mi:.3f}, pairs={len(sig_pairs)}")

    results[name] = layer_data

# Save
with open("results/sprint_011b_mi_dynamics.json", "w") as f:
    json.dump(results, f, indent=2)

# Analysis
print("\n\n=== COMPARATIVE ANALYSIS ===")
for name in results:
    data = results[name]
    total_mis = [d["total_mi"] for d in data]
    sparsities = [d["sparsity"] for d in data]

    # Find first layer with any MI
    first_mi = next((i for i, t in enumerate(total_mis) if t > 0.01), None)

    # Check for non-monotonic total MI
    diffs = [total_mis[i+1] - total_mis[i] for i in range(len(total_mis)-1)]
    decreases = [i for i, d in enumerate(diffs) if d < -0.01]

    # Check if MI spreads gradually (increasing sparsity) or all at once
    final_sparsity = sparsities[-1]
    half_sparsity_layer = next((i for i, s in enumerate(sparsities) if s >= final_sparsity / 2), None)

    print(f"\n{name}:")
    print(f"  First MI at layer: {first_mi}")
    print(f"  Total MI trajectory: {[round(t, 2) for t in total_mis]}")
    print(f"  Final sparsity: {final_sparsity:.2f} ({round(final_sparsity * n*(n-1)//2)} of {n*(n-1)//2} pairs)")
    print(f"  Half-sparsity reached at layer: {half_sparsity_layer}")
    if decreases:
        print(f"  NON-MONOTONIC! Decreases at layers: {decreases}")

    # Final MI structure
    final = data[-1]
    print(f"  Final structure: {len(final['significant_pairs'])} significant pairs")
    for (i, j, mi) in final['significant_pairs'][:5]:
        print(f"    ({i},{j}): MI={mi}")

print("\nDone! Results saved to results/sprint_011b_mi_dynamics.json")
