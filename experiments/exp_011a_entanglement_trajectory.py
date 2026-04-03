"""Sprint 011a: Entanglement trajectory — track entropy & negativity at each gate layer.

For each archetype (GHZ, W, Cluster 1D, Cluster 2D), we build the circuit gate by gate,
capturing the density matrix after each layer and computing:
- Half-cut entanglement entropy
- Half-cut logarithmic negativity
- Purity (to verify we're in pure states)

n=6 qubits throughout.
"""

import numpy as np
import json
from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.quantum_info import DensityMatrix, partial_trace, entropy, negativity

n = 6
sim = AerSimulator(method='statevector')

def get_dm(qc):
    """Get density matrix from circuit."""
    qc_copy = qc.copy()
    qc_copy.save_density_matrix()
    result = sim.run(qc_copy).result()
    return DensityMatrix(result.data()['density_matrix'])

def half_cut_entropy(dm, n):
    """Entropy of tracing out second half."""
    keep = list(range(n // 2))
    rho_A = partial_trace(dm, [i for i in range(n) if i not in keep])
    return float(entropy(rho_A, base=2))

def half_cut_negativity(dm, n):
    """Logarithmic negativity across half cut."""
    keep_A = list(range(n // 2))
    keep_B = list(range(n // 2, n))
    return float(negativity(dm, keep_A))

def purity(dm):
    return float(np.real(np.trace(dm.data @ dm.data)))

# ===== GHZ construction: H(0), then CNOT(0,1), CNOT(0,2), ..., CNOT(0,5) =====
def ghz_layers(n):
    """Return list of circuits at each layer."""
    circuits = []
    labels = []

    qc = QuantumCircuit(n)
    circuits.append(qc.copy())
    labels.append("initial |0...0>")

    qc.h(0)
    circuits.append(qc.copy())
    labels.append("H(0)")

    for i in range(1, n):
        qc.cx(0, i)
        circuits.append(qc.copy())
        labels.append(f"CNOT(0,{i})")

    return circuits, labels

# ===== W state construction =====
def w_layers(n):
    """Build W state incrementally using rotation+CNOT cascade."""
    circuits = []
    labels = []

    qc = QuantumCircuit(n)
    circuits.append(qc.copy())
    labels.append("initial |0...0>")

    # Start with |1> on qubit 0
    qc.x(0)
    circuits.append(qc.copy())
    labels.append("X(0)")

    # Distribute the excitation: at step k, rotate qubit k to share with remaining
    for k in range(n - 1):
        theta = 2 * np.arccos(np.sqrt(1.0 / (n - k)))
        qc.cry(theta, k, k + 1)
        qc.cx(k + 1, k)
        circuits.append(qc.copy())
        labels.append(f"distribute to qubit {k+1}")

    return circuits, labels

# ===== 1D Cluster construction: H on all, then CZ chain =====
def cluster1d_layers(n):
    circuits = []
    labels = []

    qc = QuantumCircuit(n)
    circuits.append(qc.copy())
    labels.append("initial |0...0>")

    # All H gates (one layer)
    for i in range(n):
        qc.h(i)
    circuits.append(qc.copy())
    labels.append("H on all (|+>^n)")

    # CZ gates one by one
    for i in range(n - 1):
        qc.cz(i, i + 1)
        circuits.append(qc.copy())
        labels.append(f"CZ({i},{i+1})")

    return circuits, labels

# ===== 2D Cluster (2x3 grid) construction =====
def cluster2d_layers(n):
    """2x3 grid: qubits 0-2 top row, 3-5 bottom row."""
    circuits = []
    labels = []

    qc = QuantumCircuit(n)
    circuits.append(qc.copy())
    labels.append("initial |0...0>")

    for i in range(n):
        qc.h(i)
    circuits.append(qc.copy())
    labels.append("H on all (|+>^n)")

    # Horizontal edges: 0-1, 1-2, 3-4, 4-5
    # Vertical edges: 0-3, 1-4, 2-5
    edges = [(0,1), (1,2), (3,4), (4,5), (0,3), (1,4), (2,5)]
    for (a, b) in edges:
        qc.cz(a, b)
        circuits.append(qc.copy())
        labels.append(f"CZ({a},{b})")

    return circuits, labels


# ===== Run all =====
results = {}

for name, layer_fn in [("GHZ", ghz_layers), ("W", w_layers),
                         ("Cluster_1D", cluster1d_layers), ("Cluster_2D", cluster2d_layers)]:
    print(f"\n=== {name} ===")
    circuits, labels = layer_fn(n)

    entropies = []
    negativities = []
    purities = []

    for i, (qc, label) in enumerate(zip(circuits, labels)):
        dm = get_dm(qc)
        e = half_cut_entropy(dm, n)
        neg = half_cut_negativity(dm, n)
        p = purity(dm)
        entropies.append(e)
        negativities.append(neg)
        purities.append(p)
        print(f"  Layer {i} [{label}]: entropy={e:.4f}, negativity={neg:.4f}, purity={p:.4f}")

    results[name] = {
        "labels": labels,
        "entropy": entropies,
        "negativity": negativities,
        "purity": purities
    }

# Save
with open("results/sprint_011a_entanglement_trajectory.json", "w") as f:
    json.dump(results, f, indent=2)

print("\n=== ANALYSIS ===")
for name in results:
    ent = results[name]["entropy"]
    neg = results[name]["negativity"]
    n_layers = len(ent)
    # Find first layer where entropy > 0.5
    first_entangled = next((i for i, e in enumerate(ent) if e > 0.1), None)
    # Find layer of max entropy
    max_layer = np.argmax(ent)
    # Check for non-monotonic behavior
    diffs = [ent[i+1] - ent[i] for i in range(len(ent)-1)]
    decreases = [i for i, d in enumerate(diffs) if d < -0.01]

    print(f"\n{name}:")
    print(f"  Layers: {n_layers}, Final entropy: {ent[-1]:.4f}, Final negativity: {neg[-1]:.4f}")
    print(f"  First entangled (S>0.1): layer {first_entangled}")
    print(f"  Max entropy at layer: {max_layer} (S={ent[max_layer]:.4f})")
    if decreases:
        print(f"  NON-MONOTONIC! Entropy decreases at layers: {decreases}")
    else:
        print(f"  Monotonically increasing entropy")

print("\nDone! Results saved to results/sprint_011a_entanglement_trajectory.json")
