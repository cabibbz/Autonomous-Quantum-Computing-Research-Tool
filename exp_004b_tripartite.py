"""
Experiment 4b: Tripartite Mutual Information
I3(A:B:C) = I(A:B) + I(A:C) - I(A:BC)
         = S(A) + S(B) + S(C) - S(AB) - S(AC) - S(BC) + S(ABC)

Key facts:
- I3 = 0 for states with only pairwise correlations
- I3 < 0 indicates genuinely tripartite correlations (monogamy of entanglement)
- I3 > 0 is possible classically but rare quantumly

Also compute I3 for ALL triples in each 6-qubit state.
"""

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy
import numpy as np
import json, os, time
from itertools import combinations

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

# Cache entropy computations
def get_entropy(sv, keep, n_total, cache):
    key = tuple(sorted(keep))
    if key not in cache:
        trace_out = [i for i in range(n_total) if i not in keep]
        if not trace_out:
            cache[key] = 0.0  # pure state
        else:
            rho = partial_trace(sv, trace_out)
            cache[key] = float(entropy(rho, base=2))
    return cache[key]

def tripartite_info(sv, a, b, c, n_total, cache):
    """I3(A:B:C) = S(A)+S(B)+S(C) - S(AB)-S(AC)-S(BC) + S(ABC)"""
    Sa = get_entropy(sv, [a], n_total, cache)
    Sb = get_entropy(sv, [b], n_total, cache)
    Sc = get_entropy(sv, [c], n_total, cache)
    Sab = get_entropy(sv, [a, b], n_total, cache)
    Sac = get_entropy(sv, [a, c], n_total, cache)
    Sbc = get_entropy(sv, [b, c], n_total, cache)
    Sabc = get_entropy(sv, [a, b, c], n_total, cache)
    return Sa + Sb + Sc - Sab - Sac - Sbc + Sabc

results = {}

for state_name, sv in [("GHZ", make_ghz(N)), ("W", make_w(N)), ("Cluster", make_cluster(N))]:
    print(f"\n=== {state_name} (n={N}) — Tripartite Information ===")
    cache = {}

    all_i3 = []
    for a, b, c in combinations(range(N), 3):
        i3 = tripartite_info(sv, a, b, c, N, cache)
        all_i3.append({"triple": [a, b, c], "I3": float(i3)})

    i3_values = [d["I3"] for d in all_i3]
    mean_i3 = np.mean(i3_values)
    min_i3 = min(i3_values)
    max_i3 = max(i3_values)
    n_negative = sum(1 for v in i3_values if v < -1e-10)
    n_zero = sum(1 for v in i3_values if abs(v) < 1e-10)
    n_positive = sum(1 for v in i3_values if v > 1e-10)

    results[state_name] = {
        "all_triples": all_i3,
        "mean_I3": float(mean_i3),
        "min_I3": float(min_i3),
        "max_I3": float(max_i3),
        "n_negative": n_negative,
        "n_zero": n_zero,
        "n_positive": n_positive,
        "total_triples": len(all_i3),
    }

    print(f"  {len(all_i3)} triples computed")
    print(f"  I3 range: [{min_i3:.4f}, {max_i3:.4f}]")
    print(f"  Mean I3: {mean_i3:.4f}")
    print(f"  Negative (genuine tripartite): {n_negative}/{len(all_i3)}")
    print(f"  Zero: {n_zero}/{len(all_i3)}")
    print(f"  Positive: {n_positive}/{len(all_i3)}")

    # Show most negative triple
    most_neg = min(all_i3, key=lambda d: d["I3"])
    most_pos = max(all_i3, key=lambda d: d["I3"])
    print(f"  Most negative: I3({most_neg['triple']}) = {most_neg['I3']:.4f}")
    if n_positive > 0:
        print(f"  Most positive: I3({most_pos['triple']}) = {most_pos['I3']:.4f}")

dt = time.time() - t0

results["runtime_seconds"] = dt

with open(os.path.join(RESULTS_DIR, "sprint_004b_tripartite.json"), "w") as f:
    json.dump(results, f, indent=2)

print(f"\nRuntime: {dt:.2f}s")
print(f"Saved to results/sprint_004b_tripartite.json")
