"""
Experiment 3a: Progressive Qubit Loss — Entanglement Dynamics
Start with 8-qubit entangled states. Trace out qubits one by one.
Track how half-cut entropy of the remaining system evolves.
Compare GHZ, W, and Cluster states.
"""

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy, DensityMatrix
import numpy as np
import json, os, time

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

t0 = time.time()

N = 8

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

def half_cut_entropy(dm_or_sv, n_qubits):
    """Compute entropy of tracing out the second half of qubits."""
    if n_qubits < 2:
        return 0.0
    half = n_qubits // 2
    rho_A = partial_trace(dm_or_sv, list(range(half, n_qubits)))
    return float(entropy(rho_A, base=2))

results = {}

for state_name, sv in [("GHZ", make_ghz(N)), ("W", make_w(N)), ("Cluster", make_cluster(N))]:
    print(f"\n=== {state_name} state, n={N} ===")
    trace = []

    # Progressive loss: trace out qubits 0, 1, 2, ... from the left (edge first)
    current_state = DensityMatrix(sv)
    remaining = N

    # Initial entropy
    ent = half_cut_entropy(current_state, remaining)
    trace.append({"qubits_lost": 0, "remaining": remaining, "half_cut_entropy": ent})
    print(f"  Lost 0: remaining={remaining}, half-cut entropy={ent:.4f}")

    for k in range(1, N - 1):  # Leave at least 2 qubits
        # Trace out qubit 0 (always the "leftmost" of what remains)
        current_state = partial_trace(current_state, [0])
        remaining = N - k
        ent = half_cut_entropy(current_state, remaining)
        trace.append({"qubits_lost": k, "remaining": remaining, "half_cut_entropy": ent})
        print(f"  Lost {k}: remaining={remaining}, half-cut entropy={ent:.4f}")

    results[state_name] = {
        "initial_qubits": N,
        "progressive_loss": trace,
    }

dt = time.time() - t0

# Analysis
print(f"\n=== Comparison ===")
for state_name in ["GHZ", "W", "Cluster"]:
    entropies = [d["half_cut_entropy"] for d in results[state_name]["progressive_loss"]]
    peak = max(entropies)
    peak_at = entropies.index(peak)
    print(f"{state_name:8s}: peak entropy={peak:.4f} at {peak_at} qubits lost, "
          f"final={entropies[-1]:.4f}")

results["runtime_seconds"] = dt

with open(os.path.join(RESULTS_DIR, "sprint_003a_progressive_loss.json"), "w") as f:
    json.dump(results, f, indent=2)

print(f"\nRuntime: {dt:.2f}s")
print(f"Saved to results/sprint_003a_progressive_loss.json")
