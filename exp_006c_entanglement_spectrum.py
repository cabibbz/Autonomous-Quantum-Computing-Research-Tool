"""Sprint 006c: Entanglement spectrum — negativity across ALL bipartitions.

For n=6, enumerate all non-trivial bipartitions (1-vs-5, 2-vs-4, 3-vs-3)
and compute negativity for each. This creates an "entanglement spectrum"
that should distinguish the three archetypes even though the half-cut
negativity is identical.

Also: concurrence monogamy. For the W state, check the Coffman-Kundu-Wootters
inequality: C²(A|BC) ≥ C²(A|B) + C²(A|C). The residual tangle τ measures
genuine tripartite entanglement.
"""

import numpy as np
import json
import time
from itertools import combinations
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, DensityMatrix, partial_trace, negativity, concurrence

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

def all_bipartitions(n):
    """Generate all bipartitions (A, B) where |A| <= |B|, up to size n//2."""
    partitions = []
    for k in range(1, n // 2 + 1):
        for combo in combinations(range(n), k):
            A = list(combo)
            B = [i for i in range(n) if i not in A]
            # For k = n//2, only keep one of each complementary pair
            if k == n // 2 and A[0] != 0:
                continue
            partitions.append((A, B))
    return partitions

def run():
    n = 6
    states = {
        'GHZ': make_ghz(n),
        'W': make_w(n),
        'Cluster': make_cluster(n),
    }

    results = {'n_qubits': n, 'bipartition_negativity': {}, 'summary': {}}

    print(f"=== Entanglement spectrum (n={n}) ===")
    partitions = all_bipartitions(n)
    print(f"Total bipartitions: {len(partitions)}")

    for name, sv in states.items():
        neg_by_size = {}  # group by |A|
        all_negs = []

        for A, B in partitions:
            k = len(A)
            neg = float(negativity(sv, A))
            all_negs.append({'A': A, 'B': B, 'size_A': k, 'negativity': neg})

            if k not in neg_by_size:
                neg_by_size[k] = []
            neg_by_size[k].append(neg)

        results['bipartition_negativity'][name] = all_negs

        # Summary statistics by partition size
        print(f"\n  {name}:")
        summary = {}
        for k in sorted(neg_by_size.keys()):
            vals = neg_by_size[k]
            mean_neg = np.mean(vals)
            std_neg = np.std(vals)
            min_neg = np.min(vals)
            max_neg = np.max(vals)
            summary[str(k)] = {
                'mean': float(mean_neg),
                'std': float(std_neg),
                'min': float(min_neg),
                'max': float(max_neg),
                'count': len(vals),
            }
            print(f"    |A|={k}: mean={mean_neg:.4f} ± {std_neg:.4f} "
                  f"[{min_neg:.4f}, {max_neg:.4f}] ({len(vals)} partitions)")
        results['summary'][name] = summary

    # Concurrence monogamy for W state (n=4 for simplicity of CKW)
    print("\n=== Concurrence monogamy (W state, n=4) ===")
    w4 = make_w(4)
    rho_w4 = DensityMatrix(w4)

    # For qubit 0 vs rest: C(0|123)
    # Need to compute concurrence of qubit 0 vs qubits 1,2,3
    # This requires the 2x8 -> we need to compute the negativity instead
    # Actually concurrence is only defined for 2-qubit states
    # Use tangle: τ_A = C²(A|B) for pure states where B is the complement
    # For a pure state, C(A|B) = 2*sqrt(det(ρ_A)) for qubit A

    # Single-qubit reduced states
    tangles = []
    pairwise_c_sq_sum = []
    for i in range(4):
        rho_i = partial_trace(w4, [j for j in range(4) if j != i])
        rho_arr = np.array(rho_i)
        # Tangle = 4*det(ρ_i) for a single qubit from pure state
        tangle_i = 4 * np.real(np.linalg.det(rho_arr))

        # Sum of squared pairwise concurrences
        c_sq_sum = 0
        for j in range(4):
            if j == i:
                continue
            trace_out = [k for k in range(4) if k not in [i, j]]
            rho_ij = partial_trace(w4, trace_out)
            c = concurrence(rho_ij)
            c_sq_sum += c**2

        residual = tangle_i - c_sq_sum
        tangles.append(float(tangle_i))
        pairwise_c_sq_sum.append(float(c_sq_sum))
        print(f"  Qubit {i}: tangle={tangle_i:.4f}, ΣC²={c_sq_sum:.4f}, "
              f"residual τ={residual:.4f}")

    results['monogamy'] = {
        'state': 'W',
        'n': 4,
        'tangles': tangles,
        'pairwise_c_sq_sums': pairwise_c_sq_sum,
        'residual_tangles': [t - s for t, s in zip(tangles, pairwise_c_sq_sum)],
    }

    # Save
    with open('results/sprint_006c_entanglement_spectrum.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nResults saved to results/sprint_006c_entanglement_spectrum.json")

if __name__ == '__main__':
    t0 = time.time()
    run()
    print(f"\nTotal time: {time.time()-t0:.1f}s")
