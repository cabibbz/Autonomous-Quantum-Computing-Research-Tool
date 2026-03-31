"""Sprint 009a: Amplitude damping (T1 energy relaxation).

Each qubit independently undergoes |1⟩→|0⟩ decay with probability γ.
Kraus operators: K0 = |0⟩⟨0| + sqrt(1-γ)|1⟩⟨1|, K1 = sqrt(γ)|0⟩⟨1|

This models spontaneous emission / T1 decay in superconducting qubits.
Unlike depolarizing noise, it's asymmetric: biased toward |0⟩.

Compare Phi and half-cut negativity decay for GHZ, W, Cluster 1D, Cluster 2D.
"""

import numpy as np
import json
import time
from itertools import combinations
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import DensityMatrix, partial_trace, entropy


def make_ghz(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(1, n):
        qc.cx(0, i)
    return DensityMatrix.from_instruction(qc)


def make_w(n):
    sv = np.zeros(2**n, dtype=complex)
    for i in range(n):
        sv[1 << i] = 1.0 / np.sqrt(n)
    from qiskit.quantum_info import Statevector
    return DensityMatrix(Statevector(sv))


def make_cluster_1d(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        qc.cz(i, i + 1)
    return DensityMatrix.from_instruction(qc)


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
    return DensityMatrix.from_instruction(qc)


def apply_amplitude_damping(dm, gamma):
    """Apply single-qubit amplitude damping to each qubit independently.

    Kraus operators per qubit:
      K0 = [[1, 0], [0, sqrt(1-γ)]]
      K1 = [[0, sqrt(γ)], [0, 0]]
    """
    n = int(np.log2(dm.data.shape[0]))
    rho = np.array(dm.data, dtype=complex)

    for qubit in range(n):
        # Build n-qubit Kraus operators for this qubit
        K0_single = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
        K1_single = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)

        # Tensor product: I ⊗ ... ⊗ K ⊗ ... ⊗ I
        K0 = np.eye(1, dtype=complex)
        K1 = np.eye(1, dtype=complex)
        for i in range(n):
            if i == qubit:
                K0 = np.kron(K0, K0_single)
                K1 = np.kron(K1, K1_single)
            else:
                K0 = np.kron(K0, np.eye(2))
                K1 = np.kron(K1, np.eye(2))

        rho = K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T

    return DensityMatrix(rho)


def all_bipartitions(n):
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
    """Phi = min mutual information over all bipartitions."""
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
    d = 2**n
    rho = dm_array.reshape([2]*2*n)
    perm = list(range(2*n))
    for q in subsys_a:
        perm[q], perm[q + n] = perm[q + n], perm[q]
    rho_pt = np.transpose(rho, perm).reshape(d, d)
    return rho_pt


def compute_half_cut_negativity(dm, n):
    half = n // 2
    rho_pt = partial_transpose_manual(np.array(dm), n, list(range(half)))
    eigs = np.linalg.eigvalsh(rho_pt)
    return float(sum(abs(e) for e in eigs if e < -1e-12))


def run():
    n = 6
    gammas = np.concatenate([
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

    results = {
        'n_qubits': n,
        'noise_type': 'amplitude_damping',
        'gammas': [round(float(g), 4) for g in gammas],
    }

    for name, dm in states.items():
        print(f"\n=== {name} ===")
        phi_values = []
        neg_values = []

        for g in gammas:
            dm_noisy = apply_amplitude_damping(dm, g)
            phi = compute_phi_mixed(dm_noisy, n)
            neg = compute_half_cut_negativity(dm_noisy, n)
            phi_values.append(round(float(phi), 6))
            neg_values.append(round(float(neg), 6))

        # Find death thresholds
        phi_death = None
        neg_death = None
        for i, g in enumerate(gammas):
            if phi_death is None and phi_values[i] < 0.01:
                phi_death = round(float(g), 4)
            if neg_death is None and neg_values[i] < 0.01:
                neg_death = round(float(g), 4)

        print(f"  Phi(γ=0)={phi_values[0]:.3f}, Neg(γ=0)={neg_values[0]:.3f}")
        print(f"  Phi death: γ≈{phi_death}, Neg death: γ≈{neg_death}")
        for i in [0, 5, 10, 15, 20]:
            if i < len(gammas):
                print(f"  γ={gammas[i]:.2f}: Phi={phi_values[i]:.4f}, Neg={neg_values[i]:.4f}")

        results[name] = {
            'phi': phi_values,
            'negativity': neg_values,
            'phi_death': phi_death,
            'neg_death': neg_death,
        }

    # Summary
    print("\n=== DEATH THRESHOLDS (Amplitude Damping) ===")
    print(f"{'State':12s} {'Phi_death':>10s} {'Neg_death':>10s}")
    for name in states:
        r = results[name]
        phi_d = r['phi_death'] or '>1.0'
        neg_d = r['neg_death'] or '>1.0'
        print(f"{name:12s} {str(phi_d):>10s} {str(neg_d):>10s}")

    with open('results/sprint_009a_amplitude_damping.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nResults saved.")


if __name__ == '__main__':
    t0 = time.time()
    run()
    print(f"\nTotal time: {time.time()-t0:.1f}s")
