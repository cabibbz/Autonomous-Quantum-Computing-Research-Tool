"""Sprint 008c: Phi under qubit loss — mixed-state integration.

Trace out 1 or 2 qubits from n=6 states, compute Phi on the remaining
mixed state. For mixed states: I(A:Ā) = S(A) + S(Ā) - S(AĀ), S(AĀ) ≠ 0.

Key questions:
- Does 2D cluster maintain higher Phi than 1D under loss?
- Is Phi loss position-dependent (like 1D entropy) or uniform (like 2D entropy)?
- How does Phi compare to entropy as a robustness measure?
"""

import numpy as np
import json
import time
from itertools import combinations
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, DensityMatrix, partial_trace, entropy


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
    """Phi for a mixed state. I(A:Ā) = S(A) + S(Ā) - S(AĀ).
    Returns (phi, mip_A) where mip_A is the minimum information partition."""
    partitions = all_bipartitions(n)
    s_total = float(entropy(dm, base=2))
    mi_min = float('inf')
    mip = None
    for A, A_bar in partitions:
        rho_A = partial_trace(dm, A_bar)
        rho_Abar = partial_trace(dm, A)
        s_A = float(entropy(rho_A, base=2))
        s_Abar = float(entropy(rho_Abar, base=2))
        mi = s_A + s_Abar - s_total
        if mi < mi_min:
            mi_min = mi
            mip = A
    return mi_min, mip


def trace_out_qubits(sv, lost_qubits, n):
    """Trace out specified qubits, return (DensityMatrix, n_remaining)."""
    remaining = [q for q in range(n) if q not in lost_qubits]
    # partial_trace takes the qubits to KEEP as complement
    dm = DensityMatrix(sv)
    dm_reduced = partial_trace(dm, lost_qubits)
    return dm_reduced, len(remaining)


def run():
    n = 6
    states = {
        'GHZ': make_ghz(n),
        'W': make_w(n),
        'Cluster_1D': make_cluster_1d(n),
        'Cluster_2D': make_cluster_2d(2, 3),
    }

    results = {'n_qubits': n}

    # --- Single qubit loss ---
    print("=" * 60)
    print("SINGLE QUBIT LOSS (n=6 -> n=5)")
    print("=" * 60)

    for name, sv in states.items():
        print(f"\n--- {name} ---")
        # Compute baseline Phi (pure state, no loss)
        dm_pure = DensityMatrix(sv)
        phi_baseline, mip_baseline = compute_phi_mixed(dm_pure, n)
        print(f"  Baseline Phi = {phi_baseline:.4f} (MIP={mip_baseline})")

        single_loss = {}
        for lost_q in range(n):
            dm_reduced, n_rem = trace_out_qubits(sv, [lost_q], n)
            phi, mip = compute_phi_mixed(dm_reduced, n_rem)
            s_total = float(entropy(dm_reduced, base=2))
            single_loss[lost_q] = {
                'phi': round(float(phi), 6),
                'mip': mip,
                's_total': round(s_total, 6),
            }
            print(f"  Lose q{lost_q}: Phi={phi:.4f}, S(total)={s_total:.4f}, MIP={mip}")

        phi_values = [v['phi'] for v in single_loss.values()]
        print(f"  Phi range: [{min(phi_values):.4f}, {max(phi_values):.4f}]")
        print(f"  Phi mean:  {np.mean(phi_values):.4f}")
        print(f"  Position-dependent: {max(phi_values) - min(phi_values) > 0.01}")

        results[f'{name}_single_loss'] = {
            'baseline_phi': round(float(phi_baseline), 6),
            'baseline_mip': mip_baseline,
            'per_qubit': {str(k): v for k, v in single_loss.items()},
            'phi_mean': round(float(np.mean(phi_values)), 6),
            'phi_min': round(float(min(phi_values)), 6),
            'phi_max': round(float(max(phi_values)), 6),
            'position_dependent': bool(max(phi_values) - min(phi_values) > 0.01),
        }

    # --- Two qubit loss ---
    print("\n" + "=" * 60)
    print("TWO QUBIT LOSS (n=6 -> n=4)")
    print("=" * 60)

    # Sample representative pairs to keep runtime manageable
    # For 6 qubits: C(6,2)=15 pairs, each with 4-qubit Phi (7 bipartitions) — fast
    for name, sv in states.items():
        print(f"\n--- {name} ---")
        two_loss = {}
        for q1, q2 in combinations(range(n), 2):
            dm_reduced, n_rem = trace_out_qubits(sv, [q1, q2], n)
            phi, mip = compute_phi_mixed(dm_reduced, n_rem)
            s_total = float(entropy(dm_reduced, base=2))
            key = f"{q1},{q2}"
            two_loss[key] = {
                'phi': round(float(phi), 6),
                'mip': mip,
                's_total': round(s_total, 6),
            }
            print(f"  Lose q{q1},q{q2}: Phi={phi:.4f}, S={s_total:.4f}")

        phi_values = [v['phi'] for v in two_loss.values()]
        print(f"  Phi range: [{min(phi_values):.4f}, {max(phi_values):.4f}]")
        print(f"  Phi mean:  {np.mean(phi_values):.4f}")

        results[f'{name}_two_loss'] = {
            'per_pair': two_loss,
            'phi_mean': round(float(np.mean(phi_values)), 6),
            'phi_min': round(float(min(phi_values)), 6),
            'phi_max': round(float(max(phi_values)), 6),
            'position_dependent': bool(max(phi_values) - min(phi_values) > 0.01),
        }

    # --- Summary comparison ---
    print("\n" + "=" * 60)
    print("SUMMARY: Phi retention under qubit loss")
    print("=" * 60)
    print(f"{'State':12s} {'Phi(0)':>7s} {'Phi(1L)':>8s} {'Phi(2L)':>8s} "
          f"{'Retain1':>8s} {'Retain2':>8s} {'PosDepend':>10s}")
    for name in states:
        phi0 = results[f'{name}_single_loss']['baseline_phi']
        phi1 = results[f'{name}_single_loss']['phi_mean']
        phi2 = results[f'{name}_two_loss']['phi_mean']
        retain1 = phi1 / phi0 if phi0 > 0 else 0
        retain2 = phi2 / phi0 if phi0 > 0 else 0
        pos1 = results[f'{name}_single_loss']['position_dependent']
        print(f"{name:12s} {phi0:7.3f} {phi1:8.3f} {phi2:8.3f} "
              f"{retain1:7.1%} {retain2:7.1%} {'YES' if pos1 else 'no':>10s}")

    with open('results/sprint_008c_phi_qubit_loss.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nResults saved to results/sprint_008c_phi_qubit_loss.json")


if __name__ == '__main__':
    t0 = time.time()
    run()
    print(f"\nTotal time: {time.time()-t0:.1f}s")
