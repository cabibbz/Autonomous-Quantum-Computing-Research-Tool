"""Sprint 008b: Phi under depolarizing noise.

For mixed states (noisy), I(A:Ā) = S(A) + S(Ā) - S(AĀ), and S(AĀ) ≠ 0.
Phi = min I(A:Ā) over all bipartitions.

Compare Phi degradation to half-cut negativity degradation (Sprint 006b).
Key question: Does integration die faster or slower than entanglement?
"""

import numpy as np
import json
import time
from itertools import combinations
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import (Statevector, DensityMatrix, partial_trace,
                                  entropy)


def make_ghz(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(1, n):
        qc.cx(0, i)
    return Statevector.from_instruction(qc)


def make_w(n):
    sv = np.zeros(2**n, dtype=complex)
    for i in range(n):
        sv[1 << i] = 1.0 / np.sqrt(n)
    return Statevector(sv)


def make_cluster_1d(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        qc.cz(i, i + 1)
    return Statevector.from_instruction(qc)


def make_cluster_2d(rows, cols):
    n = rows * cols
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for r in range(rows):
        for c in range(cols - 1):
            qc.cz(r * cols + c, r * cols + c + 1)
    for r in range(rows - 1):
        for c in range(cols):
            qc.cz(r * cols + c, (r + 1) * cols + c)
    return Statevector.from_instruction(qc)


def apply_depolarizing(dm, p):
    """Apply n-qubit depolarizing noise: rho -> (1-p)*rho + p*I/2^n."""
    d = dm.data.shape[0]
    identity = np.eye(d) / d
    noisy = (1 - p) * np.array(dm) + p * identity
    return DensityMatrix(noisy)


def all_bipartitions(n):
    """All non-trivial bipartitions, avoiding double-counting at |A|=n//2."""
    qubits = list(range(n))
    partitions = []
    for size in range(1, n // 2 + 1):
        for subset in combinations(qubits, size):
            A = list(subset)
            A_bar = [q for q in qubits if q not in A]
            if size == n // 2 and A[0] != 0:
                continue
            partitions.append((A, A_bar))
    return partitions


def compute_phi_mixed(dm, n):
    """Phi for a mixed state. I(A:Ā) = S(A) + S(Ā) - S(AĀ)."""
    partitions = all_bipartitions(n)
    s_total = float(entropy(dm, base=2))
    mi_min = float('inf')
    for A, A_bar in partitions:
        rho_A = partial_trace(dm, A_bar)
        rho_Abar = partial_trace(dm, A)
        s_A = float(entropy(rho_A, base=2))
        s_Abar = float(entropy(rho_Abar, base=2))
        mi = s_A + s_Abar - s_total
        if mi < mi_min:
            mi_min = mi
    return mi_min


def partial_transpose_manual(dm_array, n, subsys_a):
    """Manual partial transpose over subsystem A."""
    d = 2**n
    d_a = 2**len(subsys_a)
    d_b = d // d_a
    # Reshape into tensor, transpose A indices, reshape back
    # For half-cut, subsys_a = first n//2 qubits
    rho = dm_array.reshape([2]*2*n)
    # Swap the A bra/ket indices
    perm = list(range(2*n))
    for q in subsys_a:
        perm[q], perm[q + n] = perm[q + n], perm[q]
    rho_pt = np.transpose(rho, perm).reshape(d, d)
    return rho_pt


def compute_half_cut_negativity(dm, n):
    """Negativity across the half-cut."""
    half = n // 2
    rho_pt = partial_transpose_manual(np.array(dm), n, list(range(half)))
    eigs = np.linalg.eigvalsh(rho_pt)
    return float(sum(abs(e) for e in eigs if e < -1e-12))


def run():
    n = 6
    noise_levels = np.concatenate([
        np.arange(0, 0.1, 0.01),
        np.arange(0.1, 0.5, 0.05),
        np.arange(0.5, 1.01, 0.1),
    ])

    states = {
        'GHZ': make_ghz(n),
        'W': make_w(n),
        'Cluster_1D': make_cluster_1d(n),
        'Cluster_2D': make_cluster_2d(2, 3),
    }

    results = {'n_qubits': n, 'noise_levels': [round(float(p), 4) for p in noise_levels]}

    for name, sv in states.items():
        print(f"\n=== {name} ===")
        dm_pure = DensityMatrix(sv)
        phi_values = []
        neg_values = []

        for p in noise_levels:
            dm_noisy = apply_depolarizing(dm_pure, p)
            phi = compute_phi_mixed(dm_noisy, n)
            neg = compute_half_cut_negativity(dm_noisy, n)
            phi_values.append(round(float(phi), 6))
            neg_values.append(round(float(neg), 6))

        # Find death thresholds (where value first drops below 0.01)
        phi_death = None
        neg_death = None
        for i, p in enumerate(noise_levels):
            if phi_death is None and phi_values[i] < 0.01:
                phi_death = round(float(p), 4)
            if neg_death is None and neg_values[i] < 0.01:
                neg_death = round(float(p), 4)

        print(f"  Phi(p=0)={phi_values[0]:.3f}, Neg(p=0)={neg_values[0]:.3f}")
        print(f"  Phi death: p≈{phi_death}, Neg death: p≈{neg_death}")

        # Sample a few noise levels
        for i in [0, 5, 10, 15, 20]:
            if i < len(noise_levels):
                print(f"  p={noise_levels[i]:.2f}: Phi={phi_values[i]:.4f}, "
                      f"Neg={neg_values[i]:.4f}")

        results[name] = {
            'phi': phi_values,
            'negativity': neg_values,
            'phi_death': phi_death,
            'neg_death': neg_death,
        }

    # Summary
    print("\n=== DEATH THRESHOLDS ===")
    print(f"{'State':12s} {'Phi_death':>10s} {'Neg_death':>10s} {'Phi_dies_first':>15s}")
    for name in states:
        r = results[name]
        phi_d = r['phi_death'] or '>1.0'
        neg_d = r['neg_death'] or '>1.0'
        first = 'Phi' if (r['phi_death'] or 2) < (r['neg_death'] or 2) else 'Neg'
        if r['phi_death'] == r['neg_death']:
            first = 'same'
        print(f"{name:12s} {str(phi_d):>10s} {str(neg_d):>10s} {first:>15s}")

    with open('results/sprint_008b_phi_noise.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nResults saved.")


if __name__ == '__main__':
    t0 = time.time()
    run()
    print(f"\nTotal time: {time.time()-t0:.1f}s")
