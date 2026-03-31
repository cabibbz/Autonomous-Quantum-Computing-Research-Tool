"""Sprint 007a: W-state I3 (gap-fill) + 2D cluster state MI & I3.

Build 2D cluster states on a grid (2x3=6 qubits, 2x4=8 qubits).
A 2D cluster state: H on all qubits, then CZ on every grid edge.

Qubit layout (2x3):
  0 - 1 - 2
  |   |   |
  3 - 4 - 5

Qubit layout (2x4):
  0 - 1 - 2 - 3
  |   |   |   |
  4 - 5 - 6 - 7

Compute pairwise MI and tripartite information for all triples.
Compare 2D cluster to 1D chain cluster and fill W-state I3 gap.
"""

import numpy as np
import json
import time
from itertools import combinations
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy

def make_w(n):
    dim = 2**n
    vec = np.zeros(dim, dtype=complex)
    for i in range(n):
        idx = 1 << (n - 1 - i)
        vec[idx] = 1.0 / np.sqrt(n)
    return Statevector(vec)

def make_cluster_1d(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        qc.cz(i, i + 1)
    return Statevector.from_instruction(qc)

def make_cluster_2d(rows, cols):
    """2D cluster state on a rows x cols grid."""
    n = rows * cols
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    # Horizontal edges
    for r in range(rows):
        for c in range(cols - 1):
            qc.cz(r * cols + c, r * cols + c + 1)
    # Vertical edges
    for r in range(rows - 1):
        for c in range(cols):
            qc.cz(r * cols + c, (r + 1) * cols + c)
    return Statevector.from_instruction(qc)

def mutual_info(sv, i, j, n):
    """I(i:j) = S(i) + S(j) - S(i,j)"""
    all_but_i = [k for k in range(n) if k != i]
    all_but_j = [k for k in range(n) if k != j]
    all_but_ij = [k for k in range(n) if k not in (i, j)]

    rho_i = partial_trace(sv, all_but_i)
    rho_j = partial_trace(sv, all_but_j)
    rho_ij = partial_trace(sv, all_but_ij)

    return float(entropy(rho_i) + entropy(rho_j) - entropy(rho_ij))

def tripartite_info(sv, i, j, k, n):
    """I3(i:j:k) = I(i:j) + I(i:k) - I(i:jk)"""
    # I(i:jk) = S(i) + S(jk) - S(ijk)
    rest_i = [q for q in range(n) if q != i]
    rest_jk = [q for q in range(n) if q not in (j, k)]
    rest_ijk = [q for q in range(n) if q not in (i, j, k)]

    S_i = float(entropy(partial_trace(sv, rest_i)))
    S_jk = float(entropy(partial_trace(sv, rest_jk)))
    S_ijk = float(entropy(partial_trace(sv, rest_ijk)))
    I_i_jk = S_i + S_jk - S_ijk

    I_ij = mutual_info(sv, i, j, n)
    I_ik = mutual_info(sv, i, k, n)

    return I_ij + I_ik - I_i_jk

def grid_label(idx, cols):
    """Return (row, col) for a grid qubit index."""
    return (idx // cols, idx % cols)

def grid_distance(i, j, cols):
    """Manhattan distance on grid."""
    ri, ci = grid_label(i, cols)
    rj, cj = grid_label(j, cols)
    return abs(ri - rj) + abs(ci - cj)

def run():
    results = {}

    # --- Part 1: W-state I3 (gap-fill, n=6) ---
    print("=== W-state tripartite information (n=6) ===")
    n = 6
    w6 = make_w(n)
    w_i3 = []
    for triple in combinations(range(n), 3):
        i, j, k = triple
        val = tripartite_info(w6, i, j, k, n)
        w_i3.append({'triple': list(triple), 'I3': round(val, 6)})
        print(f"  I3({i},{j},{k}) = {val:.4f}")

    vals = [x['I3'] for x in w_i3]
    results['w_state_i3'] = {
        'n': n,
        'triples': w_i3,
        'mean': round(float(np.mean(vals)), 6),
        'min': round(float(np.min(vals)), 6),
        'max': round(float(np.max(vals)), 6),
    }
    print(f"  Mean I3 = {np.mean(vals):.4f}, range [{np.min(vals):.4f}, {np.max(vals):.4f}]")

    # --- Part 2: 2D cluster state (2x3) MI & I3 ---
    print("\n=== 2D Cluster State (2x3) ===")
    rows, cols = 2, 3
    n2d = rows * cols
    cluster_2d = make_cluster_2d(rows, cols)

    # Also compute 1D chain for comparison
    cluster_1d = make_cluster_1d(n2d)

    for label, sv in [('cluster_2d_2x3', cluster_2d), ('cluster_1d_6', cluster_1d)]:
        print(f"\n  --- {label} ---")
        mi_data = []
        for i, j in combinations(range(n2d), 2):
            mi = mutual_info(sv, i, j, n2d)
            dist = grid_distance(i, j, cols) if '2d' in label else abs(i - j)
            mi_data.append({
                'pair': [i, j],
                'MI': round(mi, 6),
                'distance': dist,
            })
            if mi > 0.01:
                print(f"    MI({i},{j}) = {mi:.4f} [dist={dist}]")

        i3_data = []
        for triple in combinations(range(n2d), 3):
            i, j, k = triple
            val = tripartite_info(sv, i, j, k, n2d)
            i3_data.append({'triple': list(triple), 'I3': round(val, 6)})

        i3_vals = [x['I3'] for x in i3_data]
        neg_i3 = [x for x in i3_data if x['I3'] < -0.01]
        pos_i3 = [x for x in i3_data if x['I3'] > 0.01]
        print(f"    I3: {len(neg_i3)} negative, {len(pos_i3)} positive, "
              f"range [{np.min(i3_vals):.4f}, {np.max(i3_vals):.4f}]")

        # Show most negative I3 triples
        sorted_i3 = sorted(i3_data, key=lambda x: x['I3'])
        for x in sorted_i3[:5]:
            print(f"    Most negative: I3{x['triple']} = {x['I3']:.4f}")

        results[label] = {
            'n': n2d,
            'MI': mi_data,
            'I3': i3_data,
            'I3_summary': {
                'mean': round(float(np.mean(i3_vals)), 6),
                'min': round(float(np.min(i3_vals)), 6),
                'max': round(float(np.max(i3_vals)), 6),
                'n_negative': len(neg_i3),
                'n_positive': len(pos_i3),
            }
        }

    # --- Part 3: 2D cluster (2x4) MI & I3 ---
    print("\n=== 2D Cluster State (2x4) ===")
    rows4, cols4 = 2, 4
    n2d4 = rows4 * cols4
    cluster_2d_8 = make_cluster_2d(rows4, cols4)
    cluster_1d_8 = make_cluster_1d(n2d4)

    for label, sv in [('cluster_2d_2x4', cluster_2d_8), ('cluster_1d_8', cluster_1d_8)]:
        print(f"\n  --- {label} ---")
        mi_data = []
        for i, j in combinations(range(n2d4), 2):
            mi = mutual_info(sv, i, j, n2d4)
            dist = grid_distance(i, j, cols4) if '2d' in label else abs(i - j)
            mi_data.append({
                'pair': [i, j],
                'MI': round(mi, 6),
                'distance': dist,
            })
            if mi > 0.01:
                print(f"    MI({i},{j}) = {mi:.4f} [dist={dist}]")

        # I3 for n=8 has C(8,3)=56 triples — feasible
        i3_data = []
        for triple in combinations(range(n2d4), 3):
            i, j, k = triple
            val = tripartite_info(sv, i, j, k, n2d4)
            i3_data.append({'triple': list(triple), 'I3': round(val, 6)})

        i3_vals = [x['I3'] for x in i3_data]
        neg_i3 = [x for x in i3_data if x['I3'] < -0.01]
        pos_i3 = [x for x in i3_data if x['I3'] > 0.01]
        print(f"    I3: {len(neg_i3)} negative, {len(pos_i3)} positive, "
              f"range [{np.min(i3_vals):.4f}, {np.max(i3_vals):.4f}]")

        sorted_i3 = sorted(i3_data, key=lambda x: x['I3'])
        for x in sorted_i3[:5]:
            print(f"    Most negative: I3{x['triple']} = {x['I3']:.4f}")

        results[label] = {
            'n': n2d4,
            'MI': mi_data,
            'I3': i3_data,
            'I3_summary': {
                'mean': round(float(np.mean(i3_vals)), 6),
                'min': round(float(np.min(i3_vals)), 6),
                'max': round(float(np.max(i3_vals)), 6),
                'n_negative': len(neg_i3),
                'n_positive': len(pos_i3),
            }
        }

    # Save
    with open('results/sprint_007a_2d_cluster_mi_i3.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nResults saved to results/sprint_007a_2d_cluster_mi_i3.json")

if __name__ == '__main__':
    t0 = time.time()
    run()
    print(f"\nTotal time: {time.time()-t0:.1f}s")
