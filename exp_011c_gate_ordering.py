"""Sprint 011c: Gate ordering dependence — same final state, different construction paths.

For 1D cluster state (n=6), try different orderings of CZ gates:
  - Forward: CZ(0,1), CZ(1,2), ..., CZ(4,5)
  - Reverse: CZ(4,5), CZ(3,4), ..., CZ(0,1)
  - Outside-in: CZ(0,1), CZ(4,5), CZ(1,2), CZ(3,4), CZ(2,3)
  - Inside-out: CZ(2,3), CZ(1,2), CZ(3,4), CZ(0,1), CZ(4,5)

Track: total MI, pairwise MI sparsity, max single-pair MI, half-cut entropy at each layer.
All orderings produce the SAME final state (CZ gates commute), but take different paths.

Also test GHZ with different CNOT orderings (CNOTs from qubit 0 commute too).
"""

import numpy as np
import json
from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.quantum_info import DensityMatrix, partial_trace, entropy, negativity

n = 6
sim = AerSimulator(method='statevector')

def get_dm(qc):
    qc_copy = qc.copy()
    qc_copy.save_density_matrix()
    result = sim.run(qc_copy).result()
    return DensityMatrix(result.data()['density_matrix'])

def compute_metrics(dm, n):
    """Compute half-cut entropy, total pairwise MI, MI sparsity, max MI."""
    # Single-qubit entropies
    s1 = {}
    for i in range(n):
        trace_out = [j for j in range(n) if j != i]
        rho_i = partial_trace(dm, trace_out)
        s1[i] = float(entropy(rho_i, base=2))

    # Pairwise MI
    mi_pairs = {}
    total_mi = 0.0
    max_mi = 0.0
    n_sig = 0
    for i in range(n):
        for j in range(i+1, n):
            trace_out = [k for k in range(n) if k != i and k != j]
            rho_ij = partial_trace(dm, trace_out)
            s_ij = float(entropy(rho_ij, base=2))
            mi = s1[i] + s1[j] - s_ij
            mi_pairs[(i,j)] = mi
            total_mi += mi
            max_mi = max(max_mi, mi)
            if mi > 0.01:
                n_sig += 1

    sparsity = n_sig / (n * (n-1) // 2)

    # Half-cut entropy
    keep = list(range(n // 2))
    rho_A = partial_trace(dm, [i for i in range(n) if i not in keep])
    hc_entropy = float(entropy(rho_A, base=2))

    return {
        "total_mi": round(total_mi, 4),
        "max_mi": round(max_mi, 4),
        "sparsity": round(sparsity, 4),
        "hc_entropy": round(hc_entropy, 4),
        "n_significant_pairs": n_sig
    }


def build_with_ordering(n, edge_ordering, state_name):
    """Build cluster state with given CZ ordering, track metrics at each step."""
    circuits = []
    labels = []

    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    circuits.append(qc.copy())
    labels.append("H⊗n")

    for (a, b) in edge_ordering:
        qc.cz(a, b)
        circuits.append(qc.copy())
        labels.append(f"CZ({a},{b})")

    return circuits, labels


def build_ghz_with_ordering(n, target_ordering):
    """Build GHZ with given CNOT ordering."""
    circuits = []
    labels = []

    qc = QuantumCircuit(n)
    qc.h(0)
    circuits.append(qc.copy())
    labels.append("H(0)")

    for t in target_ordering:
        qc.cx(0, t)
        circuits.append(qc.copy())
        labels.append(f"CX(0,{t})")

    return circuits, labels


# ===== Cluster 1D orderings =====
cluster_orderings = {
    "forward":    [(0,1), (1,2), (2,3), (3,4), (4,5)],
    "reverse":    [(4,5), (3,4), (2,3), (1,2), (0,1)],
    "outside_in": [(0,1), (4,5), (1,2), (3,4), (2,3)],
    "inside_out": [(2,3), (1,2), (3,4), (0,1), (4,5)],
}

# ===== GHZ orderings =====
ghz_orderings = {
    "forward":    [1, 2, 3, 4, 5],
    "reverse":    [5, 4, 3, 2, 1],
    "outside_in": [1, 5, 2, 4, 3],
    "alternating": [3, 1, 5, 2, 4],
}

results = {"cluster_1d": {}, "ghz": {}}

# Run cluster orderings
for ord_name, edges in cluster_orderings.items():
    print(f"\n=== Cluster 1D: {ord_name} ===")
    circuits, labels = build_with_ordering(n, edges, "cluster_1d")
    trajectory = []
    for idx, (qc, label) in enumerate(zip(circuits, labels)):
        dm = get_dm(qc)
        m = compute_metrics(dm, n)
        m["label"] = label
        m["layer"] = idx
        trajectory.append(m)
        print(f"  L{idx} [{label}]: MI={m['total_mi']:.2f}, spar={m['sparsity']:.2f}, hc_S={m['hc_entropy']:.2f}")
    results["cluster_1d"][ord_name] = trajectory

# Run GHZ orderings
for ord_name, targets in ghz_orderings.items():
    print(f"\n=== GHZ: {ord_name} ===")
    circuits, labels = build_ghz_with_ordering(n, targets)
    trajectory = []
    for idx, (qc, label) in enumerate(zip(circuits, labels)):
        dm = get_dm(qc)
        m = compute_metrics(dm, n)
        m["label"] = label
        m["layer"] = idx
        trajectory.append(m)
        print(f"  L{idx} [{label}]: MI={m['total_mi']:.2f}, spar={m['sparsity']:.2f}, hc_S={m['hc_entropy']:.2f}")
    results["ghz"][ord_name] = trajectory

# Save
with open("results/sprint_011c_gate_ordering.json", "w") as f:
    json.dump(results, f, indent=2)

# Analysis
print("\n\n=== KEY COMPARISONS ===")

print("\nCluster 1D — Total MI trajectory by ordering:")
for ord_name in cluster_orderings:
    traj = [d["total_mi"] for d in results["cluster_1d"][ord_name]]
    print(f"  {ord_name:12s}: {traj}")

print(f"\n  All end at same final MI: {results['cluster_1d']['forward'][-1]['total_mi']}")

# Peak MI during construction
print("\nCluster 1D — Peak total MI during construction:")
for ord_name in cluster_orderings:
    traj = [d["total_mi"] for d in results["cluster_1d"][ord_name]]
    peak = max(traj)
    peak_layer = traj.index(peak)
    print(f"  {ord_name:12s}: peak MI = {peak:.2f} at layer {peak_layer}")

# Half-cut entropy trajectories
print("\nCluster 1D — Half-cut entropy trajectories:")
for ord_name in cluster_orderings:
    traj = [d["hc_entropy"] for d in results["cluster_1d"][ord_name]]
    print(f"  {ord_name:12s}: {traj}")

print("\nGHZ — Total MI trajectory by ordering:")
for ord_name in ghz_orderings:
    traj = [d["total_mi"] for d in results["ghz"][ord_name]]
    print(f"  {ord_name:12s}: {traj}")

# Diversity metric: how different are the paths?
print("\nPath diversity (std of total MI at each layer):")
for state_type, orderings in [("cluster_1d", cluster_orderings), ("ghz", ghz_orderings)]:
    max_layers = max(len(results[state_type][o]) for o in orderings)
    for layer in range(max_layers):
        vals = [results[state_type][o][layer]["total_mi"]
                for o in orderings if layer < len(results[state_type][o])]
        if len(vals) > 1:
            std = np.std(vals)
            if std > 0.01:
                print(f"  {state_type} layer {layer}: mean={np.mean(vals):.2f}, std={std:.2f}")

print("\nDone! Results saved to results/sprint_011c_gate_ordering.json")
