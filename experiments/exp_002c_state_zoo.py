"""
Experiment 2c: Entanglement Zoo — Comparing State Families
Compare entanglement entropy of different state types at 4-8 qubits.
All use partial_trace on <=8 qubits — safe for CPU.

State families:
- GHZ: (|00..0> + |11..1>)/sqrt(2) — "fragile" all-or-nothing entanglement
- W: (|100..0> + |010..0> + ... + |00..01>)/sqrt(n) — "robust" distributed entanglement
- Cluster: linear cluster state — resource for measurement-based quantum computation
- Product: |+>^n — no entanglement (baseline)
"""

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy
import numpy as np
import json, os, time

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

t0 = time.time()

def make_ghz(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(1, n):
        qc.cx(0, i)
    return qc

def make_w(n):
    """W state: equal superposition of all single-excitation states."""
    qc = QuantumCircuit(n)
    # Build W state using a cascade of controlled rotations
    qc.x(0)  # Start with |100...0>
    for i in range(n - 1):
        # Rotate qubit i to split amplitude with remaining qubits
        theta = 2 * np.arccos(np.sqrt(1 / (n - i)))
        qc.cry(theta, i, i + 1)
        qc.cx(i + 1, i)
    return qc

def make_cluster(n):
    """Linear cluster state: |+>^n then CZ between neighbors."""
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        qc.cz(i, i + 1)
    return qc

def make_product(n):
    """Product state |+>^n — no entanglement."""
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    return qc

results = {"states": {}}

for n in [4, 6, 8]:
    print(f"\n=== n = {n} qubits ===")
    states = {
        "GHZ": make_ghz(n),
        "W": make_w(n),
        "Cluster": make_cluster(n),
        "Product": make_product(n),
    }

    for name, qc in states.items():
        sv = Statevector.from_instruction(qc)

        # Compute half-cut entropy
        half = n // 2
        rho = partial_trace(sv, list(range(half, n)))
        ent_half = float(entropy(rho, base=2))

        # Compute single-qubit entropy (trace out all but qubit 0)
        rho_1 = partial_trace(sv, list(range(1, n)))
        ent_single = float(entropy(rho_1, base=2))

        # Robustness: entropy after "losing" one qubit (trace out qubit 0)
        rho_lost = partial_trace(sv, [0])
        # Then check entropy of half of remaining system
        remaining = n - 1
        half_r = remaining // 2
        rho_lost_half = partial_trace(rho_lost, list(range(half_r, remaining)))
        ent_after_loss = float(entropy(rho_lost_half, base=2))

        key = f"{name}_n{n}"
        results["states"][key] = {
            "name": name,
            "n_qubits": n,
            "half_cut_entropy": ent_half,
            "single_qubit_entropy": ent_single,
            "entropy_after_qubit_loss": ent_after_loss,
        }

        print(f"  {name:8s}: half-cut={ent_half:.4f}, single-qubit={ent_single:.4f}, "
              f"after-loss={ent_after_loss:.4f}")

dt = time.time() - t0

# Comparative analysis
print(f"\n=== Key Observations ===")
for n in [4, 6, 8]:
    ghz = results["states"][f"GHZ_n{n}"]
    w = results["states"][f"W_n{n}"]
    cl = results["states"][f"Cluster_n{n}"]
    print(f"n={n}: GHZ half={ghz['half_cut_entropy']:.3f}, W half={w['half_cut_entropy']:.3f}, "
          f"Cluster half={cl['half_cut_entropy']:.3f}")
    print(f"      After loss: GHZ={ghz['entropy_after_qubit_loss']:.3f}, "
          f"W={w['entropy_after_qubit_loss']:.3f}, Cluster={cl['entropy_after_qubit_loss']:.3f}")

results["runtime_seconds"] = dt

with open(os.path.join(RESULTS_DIR, "sprint_002c_state_zoo.json"), "w") as f:
    json.dump(results, f, indent=2)

print(f"\nRuntime: {dt:.2f}s")
print(f"Saved to results/sprint_002c_state_zoo.json")
