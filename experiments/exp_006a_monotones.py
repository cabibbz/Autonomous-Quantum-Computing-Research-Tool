"""Sprint 006a: Pairwise concurrence & negativity for GHZ, W, Cluster states.

Computes entanglement monotones for all qubit pairs across n=4,6,8.
Completes the archetype characterization table.
"""

import numpy as np
import json
import time
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, concurrence, negativity

def make_ghz(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(n - 1):
        qc.cx(i, i + 1)
    return Statevector.from_instruction(qc)

def make_w(n):
    """W state: equal superposition of all single-excitation states."""
    dim = 2**n
    vec = np.zeros(dim, dtype=complex)
    for i in range(n):
        idx = 1 << (n - 1 - i)
        vec[idx] = 1.0 / np.sqrt(n)
    return Statevector(vec)

def make_cluster(n):
    """1D linear cluster state."""
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        qc.cz(i, i + 1)
    return Statevector.from_instruction(qc)

def pairwise_metrics(sv, n):
    """Compute concurrence and negativity for all qubit pairs."""
    results = []
    for i in range(n):
        for j in range(i + 1, n):
            # Trace out all qubits except i, j
            keep = [i, j]
            trace_out = [k for k in range(n) if k not in keep]
            rho_ij = partial_trace(sv, trace_out)

            c = concurrence(rho_ij)
            neg = negativity(rho_ij, [0])  # partial transpose over first subsystem

            results.append({
                'qubits': [i, j],
                'distance': j - i,
                'concurrence': float(c),
                'negativity': float(neg),
            })
    return results

def run():
    all_results = {}

    for n in [4, 6, 8]:
        print(f"\n=== n={n} qubits ===")
        t0 = time.time()

        states = {
            'GHZ': make_ghz(n),
            'W': make_w(n),
            'Cluster': make_cluster(n),
        }

        for name, sv in states.items():
            pairs = pairwise_metrics(sv, n)
            key = f"{name}_n{n}"
            all_results[key] = {
                'state': name,
                'n_qubits': n,
                'pairs': pairs,
                'mean_concurrence': float(np.mean([p['concurrence'] for p in pairs])),
                'max_concurrence': float(np.max([p['concurrence'] for p in pairs])),
                'mean_negativity': float(np.mean([p['negativity'] for p in pairs])),
                'max_negativity': float(np.max([p['negativity'] for p in pairs])),
            }

            # Summary
            nonzero_c = sum(1 for p in pairs if p['concurrence'] > 1e-10)
            nonzero_n = sum(1 for p in pairs if p['negativity'] > 1e-10)
            print(f"  {name}: mean_C={all_results[key]['mean_concurrence']:.4f}, "
                  f"max_C={all_results[key]['max_concurrence']:.4f}, "
                  f"nonzero_pairs={nonzero_c}/{len(pairs)}, "
                  f"mean_neg={all_results[key]['mean_negativity']:.4f}")

        elapsed = time.time() - t0
        print(f"  (took {elapsed:.1f}s)")

    # Save results
    with open('results/sprint_006a_monotones.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    print("\nResults saved to results/sprint_006a_monotones.json")

if __name__ == '__main__':
    run()
