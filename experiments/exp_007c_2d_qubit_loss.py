"""Sprint 007c: 2D cluster state under qubit loss.

2D cluster on 2x3 grid:
  0 - 1 - 2
  |   |   |
  3 - 4 - 5

Topological positions:
- Corner: qubits 0, 2, 3, 5 (degree 2)
- Edge: qubits 1, 4 (degree 3)

Compare to 1D chain (6 qubits):
  0 - 1 - 2 - 3 - 4 - 5
- Edge: 0, 5 (degree 1)
- Interior: 1, 2, 3, 4 (degree 2)

Key question: In 1D, losing qubits near the half-cut boundary caused
entropy jumps (Sprint 003). Does 2D topology create different loss
patterns? Does the higher connectivity protect entanglement?

Also test 2-qubit loss: corner+corner, corner+edge, edge+edge.
"""

import numpy as np
import json
import time
from itertools import combinations
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy, DensityMatrix

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

def entropy_after_loss(sv, qubits_to_lose, n_total):
    """Trace out qubits, return half-cut entropy of remainder."""
    dm = partial_trace(sv, qubits_to_lose)
    remaining = n_total - len(qubits_to_lose)
    if remaining < 2:
        return 0.0
    half = remaining // 2
    # Trace out second half to get rho_A
    rho_A = partial_trace(dm, list(range(half, remaining)))
    return float(entropy(rho_A, base=2))

def grid_label(idx, cols):
    return (idx // cols, idx % cols)

def grid_degree(idx, rows, cols):
    """Number of neighbors on the grid."""
    r, c = grid_label(idx, cols)
    deg = 0
    if r > 0: deg += 1
    if r < rows - 1: deg += 1
    if c > 0: deg += 1
    if c < cols - 1: deg += 1
    return deg

def position_type_2d(idx, rows, cols):
    deg = grid_degree(idx, rows, cols)
    if deg == 2:
        return "corner"
    elif deg == 3:
        return "edge"
    else:
        return "interior"  # degree 4, not present in 2x3

def position_type_1d(idx, n):
    if idx == 0 or idx == n - 1:
        return "endpoint"
    return "interior"

def run():
    n = 6
    rows, cols = 2, 3

    sv_1d = make_cluster_1d(n)
    sv_2d = make_cluster_2d(rows, cols)

    results = {'n_qubits': n, 'grid': '2x3'}

    # Baseline: no loss
    baseline_1d = entropy_after_loss(sv_1d, [], n)  # won't work — let's compute directly
    # Actually for no-loss: just compute half-cut entropy
    rho_half_1d = partial_trace(sv_1d, list(range(n // 2, n)))
    baseline_1d = float(entropy(rho_half_1d, base=2))
    rho_half_2d = partial_trace(sv_2d, list(range(n // 2, n)))
    baseline_2d = float(entropy(rho_half_2d, base=2))
    print(f"Baseline half-cut entropy: 1D={baseline_1d:.4f}, 2D={baseline_2d:.4f}")

    # === Single qubit loss ===
    print("\n=== Single qubit loss ===")

    for label, sv, pos_fn in [
        ('cluster_1d', sv_1d, lambda i: position_type_1d(i, n)),
        ('cluster_2d', sv_2d, lambda i: position_type_2d(i, rows, cols)),
    ]:
        baseline = baseline_1d if '1d' in label else baseline_2d
        print(f"\n  --- {label} ---")
        single_loss = []
        for q in range(n):
            ent = entropy_after_loss(sv, [q], n)
            pos = pos_fn(q)
            delta = ent - baseline
            single_loss.append({
                'qubit': q,
                'position': pos,
                'entropy': round(ent, 6),
                'delta': round(delta, 6),
                'grid_pos': list(grid_label(q, cols)) if '2d' in label else None,
            })
            print(f"    Lose qubit {q} ({pos:8s}): entropy={ent:.4f} (Δ={delta:+.4f})")

        results[f'{label}_single_loss'] = single_loss

        # Summary by position type
        by_type = {}
        for x in single_loss:
            t = x['position']
            if t not in by_type:
                by_type[t] = []
            by_type[t].append(x['entropy'])
        for t, vals in by_type.items():
            print(f"    {t}: mean entropy = {np.mean(vals):.4f}")

    # === Two-qubit loss ===
    print("\n=== Two-qubit loss ===")

    for label, sv, pos_fn in [
        ('cluster_1d', sv_1d, lambda i: position_type_1d(i, n)),
        ('cluster_2d', sv_2d, lambda i: position_type_2d(i, rows, cols)),
    ]:
        baseline = baseline_1d if '1d' in label else baseline_2d
        print(f"\n  --- {label} ---")
        two_loss = []
        for q1, q2 in combinations(range(n), 2):
            ent = entropy_after_loss(sv, [q1, q2], n)
            pos1, pos2 = pos_fn(q1), pos_fn(q2)
            pair_type = f"{pos1}+{pos2}"
            delta = ent - baseline
            two_loss.append({
                'qubits': [q1, q2],
                'pair_type': pair_type,
                'entropy': round(ent, 6),
                'delta': round(delta, 6),
            })
            if abs(delta) > 0.01:
                print(f"    Lose ({q1},{q2}) [{pair_type:20s}]: "
                      f"entropy={ent:.4f} (Δ={delta:+.4f})")

        results[f'{label}_two_loss'] = two_loss

        # Summary by pair type
        by_type = {}
        for x in two_loss:
            t = x['pair_type']
            if t not in by_type:
                by_type[t] = []
            by_type[t].append(x['entropy'])
        print(f"    Summary by pair type:")
        for t in sorted(by_type.keys()):
            vals = by_type[t]
            print(f"      {t:20s}: mean={np.mean(vals):.4f}, "
                  f"range=[{np.min(vals):.4f}, {np.max(vals):.4f}]")

    # === Key comparison: max entropy under loss ===
    max_1d_1 = max(x['entropy'] for x in results['cluster_1d_single_loss'])
    max_2d_1 = max(x['entropy'] for x in results['cluster_2d_single_loss'])
    max_1d_2 = max(x['entropy'] for x in results['cluster_1d_two_loss'])
    max_2d_2 = max(x['entropy'] for x in results['cluster_2d_two_loss'])

    print(f"\n=== Max entropy under loss ===")
    print(f"  1-qubit loss: 1D max={max_1d_1:.4f}, 2D max={max_2d_1:.4f}")
    print(f"  2-qubit loss: 1D max={max_1d_2:.4f}, 2D max={max_2d_2:.4f}")

    results['comparison'] = {
        'baseline_1d': round(baseline_1d, 6),
        'baseline_2d': round(baseline_2d, 6),
        'max_entropy_1loss_1d': round(max_1d_1, 6),
        'max_entropy_1loss_2d': round(max_2d_1, 6),
        'max_entropy_2loss_1d': round(max_1d_2, 6),
        'max_entropy_2loss_2d': round(max_2d_2, 6),
    }

    with open('results/sprint_007c_2d_qubit_loss.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nResults saved to results/sprint_007c_2d_qubit_loss.json")

if __name__ == '__main__':
    t0 = time.time()
    run()
    print(f"\nTotal time: {time.time()-t0:.1f}s")
