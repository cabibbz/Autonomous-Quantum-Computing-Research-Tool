"""
Experiment 4a: Quantum Mutual Information Matrix
Compute pairwise mutual information I(i:j) = S(i) + S(j) - S(ij)
for all qubit pairs in GHZ, W, and Cluster states (6 qubits).

This gives a "correlation map" — which qubits are most correlated?
"""

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy
import numpy as np
import json, os, time

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

t0 = time.time()

N = 6

def make_ghz(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(1, n):
        qc.cx(0, i)
    return Statevector.from_instruction(qc)

def make_w(n):
    qc = QuantumCircuit(n)
    qc.x(0)
    for i in range(n - 1):
        theta = 2 * np.arccos(np.sqrt(1 / (n - i)))
        qc.cry(theta, i, i + 1)
        qc.cx(i + 1, i)
    return Statevector.from_instruction(qc)

def make_cluster(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        qc.cz(i, i + 1)
    return Statevector.from_instruction(qc)

def subsystem_entropy(sv, keep_qubits, n_total):
    """Entropy of a subsystem (specified by qubits to keep)."""
    trace_out = [i for i in range(n_total) if i not in keep_qubits]
    if not trace_out:
        return 0.0  # Pure state
    rho = partial_trace(sv, trace_out)
    return float(entropy(rho, base=2))

def mutual_info(sv, i, j, n_total):
    """I(i:j) = S(i) + S(j) - S(ij)"""
    Si = subsystem_entropy(sv, [i], n_total)
    Sj = subsystem_entropy(sv, [j], n_total)
    Sij = subsystem_entropy(sv, [i, j], n_total)
    return Si + Sj - Sij

results = {}

for state_name, sv in [("GHZ", make_ghz(N)), ("W", make_w(N)), ("Cluster", make_cluster(N))]:
    print(f"\n=== {state_name} (n={N}) — Mutual Information Matrix ===")

    # Pairwise MI matrix
    mi_matrix = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            mi = mutual_info(sv, i, j, N)
            mi_matrix[i, j] = mi
            mi_matrix[j, i] = mi

    # Also compute single-qubit entropies
    single_entropies = [subsystem_entropy(sv, [i], N) for i in range(N)]

    # Total correlation: sum of all pairwise MI
    total_mi = np.sum(mi_matrix) / 2  # divide by 2 since symmetric

    # Print matrix
    print("  MI matrix (bits):")
    header = "     " + "  ".join(f"  q{j}" for j in range(N))
    print(header)
    for i in range(N):
        row = f"  q{i}:"
        for j in range(N):
            if i == j:
                row += f" {single_entropies[i]:4.2f}"  # diagonal = single-qubit entropy
            else:
                row += f" {mi_matrix[i, j]:4.2f}"
        print(row)
    print(f"  Total pairwise MI: {total_mi:.4f}")

    results[state_name] = {
        "mi_matrix": mi_matrix.tolist(),
        "single_entropies": single_entropies,
        "total_pairwise_mi": float(total_mi),
    }

dt = time.time() - t0

# Comparison
print(f"\n=== Summary ===")
for name in ["GHZ", "W", "Cluster"]:
    mi = np.array(results[name]["mi_matrix"])
    # Off-diagonal elements only
    offdiag = mi[np.triu_indices(N, k=1)]
    print(f"{name:8s}: total MI={results[name]['total_pairwise_mi']:.3f}, "
          f"mean pair MI={np.mean(offdiag):.4f}, "
          f"std={np.std(offdiag):.4f}, "
          f"max={np.max(offdiag):.4f}, min={np.min(offdiag):.4f}")

results["runtime_seconds"] = dt

with open(os.path.join(RESULTS_DIR, "sprint_004a_mutual_info.json"), "w") as f:
    json.dump(results, f, indent=2)

print(f"\nRuntime: {dt:.2f}s")
print(f"Saved to results/sprint_004a_mutual_info.json")
