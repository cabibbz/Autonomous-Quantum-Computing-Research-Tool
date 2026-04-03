"""Sprint 008a: Quantum Integrated Information (Phi) and MI spectrum.

Compute Phi = min I(A:Ā) over all bipartitions for GHZ, W, Cluster-1D,
Cluster-2D at n=6. Also compute full MI spectrum.

For pure states: I(A:Ā) = 2*S(A) since S(A)=S(Ā) and S(total)=0.
Phi = 2 * min S(A) = the weakest informational link.

Also compute normalized Phi: Phi_norm = min [I(A:Ā) / min(|A|, |Ā|)]
to account for partition size asymmetry.
"""

import numpy as np
import json
import time
from itertools import combinations
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy


def make_ghz(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(1, n):
        qc.cx(0, i)
    return Statevector.from_instruction(qc)


def make_w(n):
    """W state: equal superposition of single-excitation states."""
    sv = np.zeros(2**n, dtype=complex)
    for i in range(n):
        idx = 1 << i
        sv[idx] = 1.0 / np.sqrt(n)
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
    """Generate all non-trivial bipartitions (A, Ā) where 1 <= |A| <= n//2."""
    qubits = list(range(n))
    partitions = []
    for size in range(1, n // 2 + 1):
        for subset in combinations(qubits, size):
            A = list(subset)
            A_bar = [q for q in qubits if q not in A]
            # Avoid double-counting when |A| = n//2
            if size == n // 2 and A[0] != 0:
                continue
            partitions.append((A, A_bar))
    return partitions


def compute_mi_spectrum(sv, n):
    """Compute MI for all bipartitions. For pure states, I(A:Ā) = 2*S(A)."""
    partitions = all_bipartitions(n)
    spectrum = []
    for A, A_bar in partitions:
        # Trace out A_bar to get rho_A, compute S(A)
        rho_A = partial_trace(sv, A_bar)
        s_A = float(entropy(rho_A, base=2))
        mi = 2 * s_A  # Pure state: I(A:Ā) = 2*S(A)
        size_A = len(A)
        spectrum.append({
            'partition_A': A,
            'size_A': size_A,
            'S_A': round(s_A, 6),
            'MI': round(mi, 6),
            'MI_normalized': round(mi / min(size_A, n - size_A), 6),
        })
    return spectrum


def compute_phi(spectrum):
    """Phi = min MI over all bipartitions."""
    mi_values = [x['MI'] for x in spectrum]
    mi_norm_values = [x['MI_normalized'] for x in spectrum]
    phi = min(mi_values)
    phi_norm = min(mi_norm_values)
    # Find the MIP (minimum information partition)
    mip_idx = mi_values.index(phi)
    mip = spectrum[mip_idx]['partition_A']
    return phi, phi_norm, mip


def run():
    n = 6
    states = {
        'GHZ': make_ghz(n),
        'W': make_w(n),
        'Cluster_1D': make_cluster_1d(n),
        'Cluster_2D': make_cluster_2d(2, 3),
    }

    results = {'n_qubits': n}

    for name, sv in states.items():
        print(f"\n=== {name} ===")
        spectrum = compute_mi_spectrum(sv, n)
        phi, phi_norm, mip = compute_phi(spectrum)

        # Summary statistics
        mi_values = [x['MI'] for x in spectrum]
        s_values = [x['S_A'] for x in spectrum]

        print(f"  Phi (raw):        {phi:.4f}")
        print(f"  Phi (normalized): {phi_norm:.4f}")
        print(f"  MIP:              {mip}")
        print(f"  MI range:         [{min(mi_values):.4f}, {max(mi_values):.4f}]")
        print(f"  MI mean:          {np.mean(mi_values):.4f}")
        print(f"  MI std:           {np.std(mi_values):.4f}")
        print(f"  S(A) range:       [{min(s_values):.4f}, {max(s_values):.4f}]")

        # Breakdown by partition size
        for size in range(1, n // 2 + 1):
            size_entries = [x for x in spectrum if x['size_A'] == size]
            if size_entries:
                mi_vals = [x['MI'] for x in size_entries]
                print(f"  |A|={size}: MI mean={np.mean(mi_vals):.4f}, "
                      f"range=[{min(mi_vals):.4f}, {max(mi_vals):.4f}], "
                      f"count={len(size_entries)}")

        results[name] = {
            'phi': round(phi, 6),
            'phi_normalized': round(phi_norm, 6),
            'mip': mip,
            'mi_mean': round(float(np.mean(mi_values)), 6),
            'mi_std': round(float(np.std(mi_values)), 6),
            'mi_min': round(min(mi_values), 6),
            'mi_max': round(max(mi_values), 6),
            'spectrum': spectrum,
        }

    # Comparison table
    print("\n=== COMPARISON ===")
    print(f"{'State':12s} {'Phi':>6s} {'Phi_n':>6s} {'MI_min':>7s} {'MI_max':>7s} "
          f"{'MI_mean':>8s} {'MI_std':>7s}")
    for name in states:
        r = results[name]
        print(f"{name:12s} {r['phi']:6.3f} {r['phi_normalized']:6.3f} "
              f"{r['mi_min']:7.3f} {r['mi_max']:7.3f} "
              f"{r['mi_mean']:8.3f} {r['mi_std']:7.3f}")

    with open('results/sprint_008a_phi_spectrum.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nResults saved.")


if __name__ == '__main__':
    t0 = time.time()
    run()
    print(f"\nTotal time: {time.time()-t0:.1f}s")
