"""Sprint 009c: Noise channel fingerprint — unified comparison.

Compare depolarizing, amplitude damping, and phase damping at matched
noise strengths. Build a noise-state discrimination matrix.

Key questions:
1. Do different states have different "noise fingerprints"?
2. Which state-noise combination is most/least resilient?
3. Does any state have a noise channel it's immune to?

We use two measures: Phi (integration) and half-cut negativity (entanglement).
"""

import numpy as np
import json
import time
from itertools import combinations
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import DensityMatrix, partial_trace, entropy, Statevector


# State construction
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


# Noise channels
def apply_depolarizing(dm, p):
    d = dm.data.shape[0]
    identity = np.eye(d) / d
    return DensityMatrix((1 - p) * np.array(dm) + p * identity)


def apply_amplitude_damping(dm, gamma):
    n = int(np.log2(dm.data.shape[0]))
    rho = np.array(dm.data, dtype=complex)
    for qubit in range(n):
        K0_s = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
        K1_s = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
        K0, K1 = np.eye(1, dtype=complex), np.eye(1, dtype=complex)
        for i in range(n):
            if i == qubit:
                K0, K1 = np.kron(K0, K0_s), np.kron(K1, K1_s)
            else:
                K0, K1 = np.kron(K0, np.eye(2)), np.kron(K1, np.eye(2))
        rho = K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T
    return DensityMatrix(rho)


def apply_phase_damping(dm, lam):
    n = int(np.log2(dm.data.shape[0]))
    rho = np.array(dm.data, dtype=complex)
    for qubit in range(n):
        K0_s = np.sqrt(1 - lam) * np.eye(2, dtype=complex)
        K1_s = np.sqrt(lam) * np.array([[1, 0], [0, 0]], dtype=complex)
        K2_s = np.sqrt(lam) * np.array([[0, 0], [0, 1]], dtype=complex)
        K0, K1, K2 = np.eye(1, dtype=complex), np.eye(1, dtype=complex), np.eye(1, dtype=complex)
        for i in range(n):
            if i == qubit:
                K0, K1, K2 = np.kron(K0, K0_s), np.kron(K1, K1_s), np.kron(K2, K2_s)
            else:
                K0, K1, K2 = np.kron(K0, np.eye(2)), np.kron(K1, np.eye(2)), np.kron(K2, np.eye(2))
        rho = K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T + K2 @ rho @ K2.conj().T
    return DensityMatrix(rho)


# Measures
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
    return np.transpose(rho, perm).reshape(d, d)


def compute_half_cut_negativity(dm, n):
    half = n // 2
    rho_pt = partial_transpose_manual(np.array(dm), n, list(range(half)))
    eigs = np.linalg.eigvalsh(rho_pt)
    return float(sum(abs(e) for e in eigs if e < -1e-12))


def run():
    n = 6
    # Use a smaller set of noise levels for 3x4 = 12 state-channel combos
    noise_levels = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]

    states = {
        'GHZ': make_ghz(n),
        'W': make_w(n),
        'Cluster_1D': make_cluster_1d(n),
        'Cluster_2D': make_cluster_2d(2, 3),
    }

    channels = {
        'depolarizing': apply_depolarizing,
        'amplitude_damping': apply_amplitude_damping,
        'phase_damping': apply_phase_damping,
    }

    results = {
        'n_qubits': n,
        'noise_levels': noise_levels,
    }

    # Compute Phi and Neg for every (state, channel, noise_level) triple
    for state_name, dm in states.items():
        results[state_name] = {}
        for ch_name, ch_fn in channels.items():
            phi_vals = []
            neg_vals = []
            for p in noise_levels:
                dm_noisy = ch_fn(dm, p)
                phi = compute_phi_mixed(dm_noisy, n)
                neg = compute_half_cut_negativity(dm_noisy, n)
                phi_vals.append(round(float(phi), 6))
                neg_vals.append(round(float(neg), 6))
            results[state_name][ch_name] = {
                'phi': phi_vals,
                'negativity': neg_vals,
            }
            print(f"{state_name:12s} + {ch_name:20s}: "
                  f"Phi[0]={phi_vals[0]:.2f}→Phi[0.5]={phi_vals[7]:.3f}  "
                  f"Neg[0]={neg_vals[0]:.2f}→Neg[0.5]={neg_vals[7]:.3f}")

    # Build discrimination matrix: Phi retention at p=0.3
    p_idx = noise_levels.index(0.3)
    print(f"\n=== PHI RETENTION AT p=0.3 (fraction of pure-state Phi) ===")
    print(f"{'State':12s} {'Depolar':>10s} {'Amp.Damp':>10s} {'Phase.D':>10s} {'Most_robust':>15s}")
    fingerprints = {}
    for state_name in states:
        phi0 = results[state_name]['depolarizing']['phi'][0]
        if phi0 == 0:
            phi0 = 1e-10
        row = {}
        for ch_name in channels:
            phi_at_p = results[state_name][ch_name]['phi'][p_idx]
            row[ch_name] = round(phi_at_p / phi0, 3)
        best = max(row, key=row.get)
        fingerprints[state_name] = row
        print(f"{state_name:12s} {row['depolarizing']:>10.3f} {row['amplitude_damping']:>10.3f} "
              f"{row['phase_damping']:>10.3f} {best:>15s}")

    results['fingerprints_phi_retention_p03'] = fingerprints

    # Negativity retention at p=0.3
    print(f"\n=== NEGATIVITY RETENTION AT p=0.3 ===")
    print(f"{'State':12s} {'Depolar':>10s} {'Amp.Damp':>10s} {'Phase.D':>10s} {'Most_robust':>15s}")
    neg_fingerprints = {}
    for state_name in states:
        neg0 = results[state_name]['depolarizing']['negativity'][0]
        if neg0 == 0:
            neg0 = 1e-10
        row = {}
        for ch_name in channels:
            neg_at_p = results[state_name][ch_name]['negativity'][p_idx]
            row[ch_name] = round(neg_at_p / neg0, 3)
        best = max(row, key=row.get)
        neg_fingerprints[state_name] = row
        print(f"{state_name:12s} {row['depolarizing']:>10.3f} {row['amplitude_damping']:>10.3f} "
              f"{row['phase_damping']:>10.3f} {best:>15s}")

    results['fingerprints_neg_retention_p03'] = neg_fingerprints

    # Asymmetry score: how different is each state's response across channels?
    print(f"\n=== NOISE ASYMMETRY (std across channels / mean) ===")
    for state_name in states:
        vals = list(fingerprints[state_name].values())
        asym = np.std(vals) / (np.mean(vals) + 1e-10)
        print(f"  {state_name:12s}: {asym:.3f}")

    with open('results/sprint_009c_noise_fingerprint.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nResults saved.")


if __name__ == '__main__':
    t0 = time.time()
    run()
    print(f"\nTotal time: {time.time()-t0:.1f}s")
