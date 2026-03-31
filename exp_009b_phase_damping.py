"""Sprint 009b: Phase damping (T2 dephasing).

Each qubit independently loses phase coherence without energy loss.
Kraus operators: K0 = sqrt(1-λ) I, K1 = sqrt(λ) |0⟩⟨0|, K2 = sqrt(λ) |1⟩⟨1|
Equivalently: off-diagonal elements of each qubit decay by (1-λ).

This is the DOMINANT noise in superconducting qubits (T2 < T1 always).
It kills superpositions but preserves populations — opposite character to
amplitude damping which changes populations but preserves some coherence.

Key question: Do states that survive amplitude damping also survive dephasing?
Or is noise resilience noise-type-specific?
"""

import numpy as np
import json
import time
from itertools import combinations
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import DensityMatrix, partial_trace, entropy, Statevector


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


def apply_phase_damping(dm, lam):
    """Apply single-qubit phase damping to each qubit independently.

    Kraus operators per qubit:
      K0 = sqrt(1-λ) * I
      K1 = sqrt(λ) * |0⟩⟨0| = sqrt(λ) * [[1,0],[0,0]]
      K2 = sqrt(λ) * |1⟩⟨1| = sqrt(λ) * [[0,0],[0,1]]
    """
    n = int(np.log2(dm.data.shape[0]))
    rho = np.array(dm.data, dtype=complex)

    for qubit in range(n):
        K0_s = np.sqrt(1 - lam) * np.eye(2, dtype=complex)
        K1_s = np.sqrt(lam) * np.array([[1, 0], [0, 0]], dtype=complex)
        K2_s = np.sqrt(lam) * np.array([[0, 0], [0, 1]], dtype=complex)

        K0 = np.eye(1, dtype=complex)
        K1 = np.eye(1, dtype=complex)
        K2 = np.eye(1, dtype=complex)
        for i in range(n):
            if i == qubit:
                K0 = np.kron(K0, K0_s)
                K1 = np.kron(K1, K1_s)
                K2 = np.kron(K2, K2_s)
            else:
                K0 = np.kron(K0, np.eye(2))
                K1 = np.kron(K1, np.eye(2))
                K2 = np.kron(K2, np.eye(2))

        rho = K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T + K2 @ rho @ K2.conj().T

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
    lambdas = np.concatenate([
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
        'noise_type': 'phase_damping',
        'lambdas': [round(float(l), 4) for l in lambdas],
    }

    for name, dm in states.items():
        print(f"\n=== {name} ===")
        phi_values = []
        neg_values = []

        for lam in lambdas:
            dm_noisy = apply_phase_damping(dm, lam)
            phi = compute_phi_mixed(dm_noisy, n)
            neg = compute_half_cut_negativity(dm_noisy, n)
            phi_values.append(round(float(phi), 6))
            neg_values.append(round(float(neg), 6))

        phi_death = None
        neg_death = None
        for i, lam in enumerate(lambdas):
            if phi_death is None and phi_values[i] < 0.01:
                phi_death = round(float(lam), 4)
            if neg_death is None and neg_values[i] < 0.01:
                neg_death = round(float(lam), 4)

        print(f"  Phi(λ=0)={phi_values[0]:.3f}, Neg(λ=0)={neg_values[0]:.3f}")
        print(f"  Phi death: λ≈{phi_death}, Neg death: λ≈{neg_death}")
        for i in [0, 5, 10, 15, 20]:
            if i < len(lambdas):
                print(f"  λ={lambdas[i]:.2f}: Phi={phi_values[i]:.4f}, Neg={neg_values[i]:.4f}")

        results[name] = {
            'phi': phi_values,
            'negativity': neg_values,
            'phi_death': phi_death,
            'neg_death': neg_death,
        }

    # Summary
    print("\n=== DEATH THRESHOLDS (Phase Damping) ===")
    print(f"{'State':12s} {'Phi_death':>10s} {'Neg_death':>10s}")
    for name in states:
        r = results[name]
        phi_d = r['phi_death'] or '>1.0'
        neg_d = r['neg_death'] or '>1.0'
        print(f"{name:12s} {str(phi_d):>10s} {str(neg_d):>10s}")

    with open('results/sprint_009b_phase_damping.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nResults saved.")


if __name__ == '__main__':
    t0 = time.time()
    run()
    print(f"\nTotal time: {time.time()-t0:.1f}s")
