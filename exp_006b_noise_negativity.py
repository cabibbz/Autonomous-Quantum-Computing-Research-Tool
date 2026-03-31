"""Sprint 006b: Entanglement monotones under depolarizing noise.

Since GHZ and Cluster have zero PAIRWISE entanglement, we measure:
1. Half-cut negativity (global bipartition) under noise for all three archetypes
2. W pairwise concurrence under noise (only state with nonzero pairwise entanglement)

Goal: Find the critical noise threshold where entanglement vanishes for each archetype.
"""

import numpy as np
import json
import time
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import (
    Statevector, DensityMatrix, partial_trace, concurrence, negativity
)

def make_ghz(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(n - 1):
        qc.cx(i, i + 1)
    return Statevector.from_instruction(qc)

def make_w(n):
    dim = 2**n
    vec = np.zeros(dim, dtype=complex)
    for i in range(n):
        idx = 1 << (n - 1 - i)
        vec[idx] = 1.0 / np.sqrt(n)
    return Statevector(vec)

def make_cluster(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        qc.cz(i, i + 1)
    return Statevector.from_instruction(qc)

def depolarize(rho, p):
    """Apply depolarizing noise: rho -> (1-p)*rho + p*I/d"""
    d = np.array(rho).shape[0]
    noisy = (1 - p) * np.array(rho) + p * np.eye(d) / d
    return DensityMatrix(noisy)

def half_cut_negativity(rho_or_sv, n):
    """Negativity across the half-cut bipartition."""
    half = n // 2
    qargs = list(range(half))  # partial transpose over first half
    return negativity(rho_or_sv, qargs)

def run():
    n = 6  # Fixed size for noise study
    noise_levels = np.linspace(0, 1.0, 51)  # 0 to 100% in 2% steps

    states = {
        'GHZ': make_ghz(n),
        'W': make_w(n),
        'Cluster': make_cluster(n),
    }

    results = {
        'n_qubits': n,
        'noise_levels': noise_levels.tolist(),
        'half_cut_negativity': {},
        'w_pairwise_concurrence': {},
    }

    # 1. Half-cut negativity under noise
    print("=== Half-cut negativity under noise (n=6) ===")
    for name, sv in states.items():
        rho = DensityMatrix(sv)
        neg_values = []
        for p in noise_levels:
            noisy = depolarize(rho, p)
            neg = half_cut_negativity(noisy, n)
            neg_values.append(float(neg))

        results['half_cut_negativity'][name] = neg_values

        # Find critical threshold (where negativity first hits ~0)
        threshold = None
        for i, neg in enumerate(neg_values):
            if neg < 1e-6:
                threshold = float(noise_levels[i])
                break

        clean_neg = neg_values[0]
        print(f"  {name}: clean={clean_neg:.4f}, threshold={threshold}")

    # 2. W pairwise concurrence under noise
    print("\n=== W pairwise concurrence under noise (n=6) ===")
    w_rho = DensityMatrix(states['W'])
    mean_c_values = []
    for p in noise_levels:
        noisy = depolarize(w_rho, p)
        # Compute mean concurrence over all pairs
        concurrences = []
        for i in range(n):
            for j in range(i + 1, n):
                trace_out = [k for k in range(n) if k not in [i, j]]
                rho_ij = partial_trace(noisy, trace_out)
                c = concurrence(rho_ij)
                concurrences.append(float(c))
        mean_c_values.append(float(np.mean(concurrences)))

    results['w_pairwise_concurrence'] = mean_c_values

    # Find W concurrence threshold
    w_threshold = None
    for i, c in enumerate(mean_c_values):
        if c < 1e-6:
            w_threshold = float(noise_levels[i])
            break
    print(f"  W pairwise concurrence: clean={mean_c_values[0]:.4f}, threshold={w_threshold}")

    # Save
    with open('results/sprint_006b_noise_negativity.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nResults saved to results/sprint_006b_noise_negativity.json")

if __name__ == '__main__':
    t0 = time.time()
    run()
    print(f"\nTotal time: {time.time()-t0:.1f}s")
