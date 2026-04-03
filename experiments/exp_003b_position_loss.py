"""
Experiment 3b: Position-Dependent Qubit Loss
For an 8-qubit cluster state, does it matter WHICH qubit you trace out?
Edge qubits (0, 7) vs interior qubits (3, 4) vs off-center (1, 2).
Also test: what happens when you lose TWO specific qubits?
"""

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy, DensityMatrix
import numpy as np
import json, os, time
from itertools import combinations

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

t0 = time.time()

N = 8

def make_cluster(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        qc.cz(i, i + 1)
    return Statevector.from_instruction(qc)

def entropy_after_loss(sv, qubits_to_lose, n_total):
    """Trace out specified qubits, return half-cut entropy of remainder."""
    dm = partial_trace(sv, qubits_to_lose)
    remaining = n_total - len(qubits_to_lose)
    if remaining < 2:
        return 0.0
    half = remaining // 2
    rho_A = partial_trace(dm, list(range(half, remaining)))
    return float(entropy(rho_A, base=2))

sv = make_cluster(N)

# Part 1: Lose a single qubit at each position
print("=== Single qubit loss (8-qubit linear cluster) ===")
single_loss = {}
for q in range(N):
    ent = entropy_after_loss(sv, [q], N)
    position = "edge" if q in [0, 7] else ("near-edge" if q in [1, 6] else "interior")
    single_loss[q] = {"entropy": ent, "position": position}
    print(f"  Lose qubit {q} ({position:10s}): half-cut entropy = {ent:.4f}")

# Part 2: Lose two qubits — all pairs
print(f"\n=== Two-qubit loss (all {N*(N-1)//2} pairs) ===")
two_loss = {}
for q1, q2 in combinations(range(N), 2):
    ent = entropy_after_loss(sv, [q1, q2], N)
    key = f"{q1},{q2}"
    two_loss[key] = ent

# Show as a matrix-like format
print("     ", "  ".join(f"  {j}" for j in range(1, N)))
for i in range(N - 1):
    row = f"  {i}:"
    for j in range(1, N):
        if j > i:
            row += f" {two_loss[f'{i},{j}']:.2f}"
        else:
            row += "     "
    print(row)

# Part 3: Lose adjacent vs non-adjacent pairs
print(f"\n=== Adjacent vs non-adjacent pair loss ===")
adjacent_ents = []
nonadjacent_ents = []
for key, ent in two_loss.items():
    q1, q2 = map(int, key.split(","))
    if abs(q1 - q2) == 1:
        adjacent_ents.append(ent)
    else:
        nonadjacent_ents.append(ent)

print(f"  Adjacent pairs: mean entropy = {np.mean(adjacent_ents):.4f} (n={len(adjacent_ents)})")
print(f"  Non-adjacent:   mean entropy = {np.mean(nonadjacent_ents):.4f} (n={len(nonadjacent_ents)})")

dt = time.time() - t0

results = {
    "experiment": "position_dependent_qubit_loss",
    "state": "linear_cluster_8",
    "single_loss": single_loss,
    "two_loss": two_loss,
    "adjacent_pair_mean_entropy": float(np.mean(adjacent_ents)),
    "nonadjacent_pair_mean_entropy": float(np.mean(nonadjacent_ents)),
    "runtime_seconds": dt,
}

with open(os.path.join(RESULTS_DIR, "sprint_003b_position_loss.json"), "w") as f:
    json.dump(results, f, indent=2)

print(f"\nRuntime: {dt:.2f}s")
print(f"Saved to results/sprint_003b_position_loss.json")
