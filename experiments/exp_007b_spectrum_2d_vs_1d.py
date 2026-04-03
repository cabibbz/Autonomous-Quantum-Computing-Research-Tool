"""Sprint 007b: Entanglement spectrum — 2D cluster vs 1D cluster.

Compare negativity across all bipartitions for:
- 1D chain cluster (6 qubits)
- 2D grid cluster (2x3 = 6 qubits)

Does 2D geometry create a richer, more varied entanglement spectrum?

Qubit layout (2x3):
  0 - 1 - 2
  |   |   |
  3 - 4 - 5
"""

import numpy as np
import json
import time
from itertools import combinations
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, negativity

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

def grid_label(idx, cols):
    return (idx // cols, idx % cols)

def all_bipartitions(n):
    """Generate all bipartitions (A, B) where |A| <= |B|."""
    partitions = []
    for k in range(1, n // 2 + 1):
        for combo in combinations(range(n), k):
            A = list(combo)
            B = [i for i in range(n) if i not in A]
            if k == n // 2 and A[0] != 0:
                continue
            partitions.append((A, B))
    return partitions

def classify_partition_2d(A, cols):
    """Classify a 2D partition by whether A is contiguous on the grid."""
    positions = [grid_label(i, cols) for i in A]
    # Check if all positions form a connected subgraph on the grid
    if len(A) == 1:
        return "single"
    # Simple adjacency check
    visited = {positions[0]}
    queue = [positions[0]]
    while queue:
        r, c = queue.pop(0)
        for dr, dc in [(-1,0),(1,0),(0,-1),(0,1)]:
            nb = (r+dr, c+dc)
            if nb in set(positions) and nb not in visited:
                visited.add(nb)
                queue.append(nb)
    connected = len(visited) == len(positions)
    return f"|A|={len(A)}_{'connected' if connected else 'disconnected'}"

def run():
    n = 6
    rows, cols = 2, 3

    cluster_1d = make_cluster_1d(n)
    cluster_2d = make_cluster_2d(rows, cols)

    partitions = all_bipartitions(n)
    print(f"Total bipartitions for n={n}: {len(partitions)}")

    results = {'n_qubits': n, 'grid': '2x3'}

    for label, sv in [('cluster_1d_6', cluster_1d), ('cluster_2d_2x3', cluster_2d)]:
        is_2d = '2d' in label
        neg_by_size = {}
        all_negs = []

        for A, B in partitions:
            k = len(A)
            neg = float(negativity(sv, A))
            entry = {
                'A': A, 'B': B, 'size_A': k,
                'negativity': round(neg, 6),
            }
            if is_2d:
                entry['topology'] = classify_partition_2d(A, cols)
            all_negs.append(entry)

            if k not in neg_by_size:
                neg_by_size[k] = []
            neg_by_size[k].append(neg)

        results[label] = {'bipartitions': all_negs}

        print(f"\n  === {label} ===")
        summary = {}
        for k in sorted(neg_by_size.keys()):
            vals = neg_by_size[k]
            s = {
                'mean': round(float(np.mean(vals)), 6),
                'std': round(float(np.std(vals)), 6),
                'min': round(float(np.min(vals)), 6),
                'max': round(float(np.max(vals)), 6),
                'count': len(vals),
            }
            summary[str(k)] = s
            print(f"    |A|={k}: mean={s['mean']:.4f} ± {s['std']:.4f} "
                  f"[{s['min']:.4f}, {s['max']:.4f}] ({s['count']} partitions)")
        results[label]['summary'] = summary

        # Show top 5 most entangled bipartitions
        sorted_negs = sorted(all_negs, key=lambda x: -x['negativity'])
        print(f"    Top-5 most entangled:")
        for x in sorted_negs[:5]:
            extra = f" [{x['topology']}]" if 'topology' in x else ""
            print(f"      A={x['A']}: neg={x['negativity']:.4f}{extra}")

        # For 2D: compare connected vs disconnected partitions
        if is_2d:
            connected = [x for x in all_negs if 'connected' in x.get('topology', '')]
            disconnected = [x for x in all_negs if 'disconnected' in x.get('topology', '')]
            if connected:
                mean_c = np.mean([x['negativity'] for x in connected])
                print(f"    Connected partitions: mean neg = {mean_c:.4f} (n={len(connected)})")
            if disconnected:
                mean_d = np.mean([x['negativity'] for x in disconnected])
                print(f"    Disconnected partitions: mean neg = {mean_d:.4f} (n={len(disconnected)})")
            results[label]['connected_vs_disconnected'] = {
                'connected_mean': round(float(np.mean([x['negativity'] for x in connected])), 6) if connected else None,
                'disconnected_mean': round(float(np.mean([x['negativity'] for x in disconnected])), 6) if disconnected else None,
                'n_connected': len(connected),
                'n_disconnected': len(disconnected),
            }

    # Comparison: how much more varied is 2D?
    negs_1d = [x['negativity'] for x in results['cluster_1d_6']['bipartitions']]
    negs_2d = [x['negativity'] for x in results['cluster_2d_2x3']['bipartitions']]
    print(f"\n  === Comparison ===")
    print(f"    1D spectrum: range [{min(negs_1d):.4f}, {max(negs_1d):.4f}], "
          f"std={np.std(negs_1d):.4f}")
    print(f"    2D spectrum: range [{min(negs_2d):.4f}, {max(negs_2d):.4f}], "
          f"std={np.std(negs_2d):.4f}")
    print(f"    2D spectrum is {'MORE' if np.std(negs_2d) > np.std(negs_1d) else 'LESS'} "
          f"varied than 1D")

    results['comparison'] = {
        '1d_std': round(float(np.std(negs_1d)), 6),
        '2d_std': round(float(np.std(negs_2d)), 6),
        '1d_range': [round(float(min(negs_1d)), 6), round(float(max(negs_1d)), 6)],
        '2d_range': [round(float(min(negs_2d)), 6), round(float(max(negs_2d)), 6)],
    }

    with open('results/sprint_007b_spectrum_2d_vs_1d.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nResults saved to results/sprint_007b_spectrum_2d_vs_1d.json")

if __name__ == '__main__':
    t0 = time.time()
    run()
    print(f"\nTotal time: {time.time()-t0:.1f}s")
