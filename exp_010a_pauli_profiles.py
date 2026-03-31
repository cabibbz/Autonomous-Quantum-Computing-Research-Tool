"""Sprint 010a: Pauli expectation value profiles for entanglement archetypes.

Measure all 1-qubit and 2-qubit Pauli observables for GHZ, W, Cluster (1D, 2D).
These profiles are the "fingerprint" visible from local measurements only.
"""

import numpy as np
import json
from itertools import product
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, SparsePauliOp

def make_ghz(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(n-1):
        qc.cx(i, i+1)
    return Statevector.from_instruction(qc)

def make_w(n):
    # W = equal superposition of single-excitation states
    import numpy as np
    state = np.zeros(2**n, dtype=complex)
    for i in range(n):
        idx = 1 << i
        state[idx] = 1.0 / np.sqrt(n)
    return Statevector(state)

def make_cluster_1d(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(n-1):
        qc.cz(i, i+1)
    return Statevector.from_instruction(qc)

def make_cluster_2d(rows, cols):
    n = rows * cols
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for r in range(rows):
        for c in range(cols):
            idx = r * cols + c
            if c + 1 < cols:
                qc.cz(idx, idx + 1)
            if r + 1 < rows:
                qc.cz(idx, idx + cols)
    return Statevector.from_instruction(qc)

def pauli_expectation(sv, pauli_str):
    """Compute <psi|P|psi> for a Pauli string like 'XIZI'."""
    op = SparsePauliOp.from_list([(pauli_str, 1.0)])
    return sv.expectation_value(op).real

def get_pauli_profile(sv, n):
    """Compute all 1-qubit and 2-qubit Pauli expectations."""
    paulis_1q = ['X', 'Y', 'Z']
    paulis_2q = ['XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ']

    results = {'1q': {}, '2q': {}}

    # 1-qubit expectations
    for q in range(n):
        for p in paulis_1q:
            # Build n-qubit Pauli string with p on qubit q, I elsewhere
            # Qiskit uses little-endian: rightmost = qubit 0
            label = ['I'] * n
            label[n - 1 - q] = p
            pstr = ''.join(label)
            val = pauli_expectation(sv, pstr)
            key = f"q{q}_{p}"
            results['1q'][key] = round(val, 6)

    # 2-qubit expectations (all pairs)
    for q1 in range(n):
        for q2 in range(q1+1, n):
            for p in paulis_2q:
                label = ['I'] * n
                label[n - 1 - q1] = p[0]
                label[n - 1 - q2] = p[1]
                pstr = ''.join(label)
                val = pauli_expectation(sv, pstr)
                key = f"q{q1}q{q2}_{p}"
                results['2q'][key] = round(val, 6)

    return results

def summarize_profile(profile, n):
    """Get statistics of Pauli expectations."""
    # 1-qubit summary by Pauli type
    summary = {}
    for p in ['X', 'Y', 'Z']:
        vals = [v for k, v in profile['1q'].items() if k.endswith(f'_{p}')]
        summary[f'1q_{p}_mean'] = round(np.mean(vals), 6)
        summary[f'1q_{p}_max'] = round(np.max(np.abs(vals)), 6)
        summary[f'1q_{p}_std'] = round(np.std(vals), 6)

    # 2-qubit summary by Pauli type
    for p in ['XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ']:
        vals = [v for k, v in profile['2q'].items() if k.endswith(f'_{p}')]
        if vals:
            summary[f'2q_{p}_mean'] = round(np.mean(vals), 6)
            summary[f'2q_{p}_max'] = round(np.max(np.abs(vals)), 6)
            summary[f'2q_{p}_std'] = round(np.std(vals), 6)

    # Count nonzero expectations (threshold 0.01)
    n_nonzero_1q = sum(1 for v in profile['1q'].values() if abs(v) > 0.01)
    n_nonzero_2q = sum(1 for v in profile['2q'].values() if abs(v) > 0.01)
    summary['n_nonzero_1q'] = n_nonzero_1q
    summary['n_nonzero_2q'] = n_nonzero_2q
    summary['total_1q'] = len(profile['1q'])
    summary['total_2q'] = len(profile['2q'])
    summary['sparsity_1q'] = round(1 - n_nonzero_1q / len(profile['1q']), 4)
    summary['sparsity_2q'] = round(1 - n_nonzero_2q / len(profile['2q']), 4)

    return summary

n = 6  # 6 qubits — tractable for all operations

states = {
    'GHZ': make_ghz(n),
    'W': make_w(n),
    'Cluster_1D': make_cluster_1d(n),
    'Cluster_2D': make_cluster_2d(2, 3),  # 2x3 = 6 qubits
}

results = {}
for name, sv in states.items():
    print(f"\n=== {name} ===")
    profile = get_pauli_profile(sv, n)
    summary = summarize_profile(profile, n)
    results[name] = {'profile': profile, 'summary': summary}

    print(f"  1-qubit nonzero: {summary['n_nonzero_1q']}/{summary['total_1q']} (sparsity {summary['sparsity_1q']})")
    print(f"  2-qubit nonzero: {summary['n_nonzero_2q']}/{summary['total_2q']} (sparsity {summary['sparsity_2q']})")
    print(f"  1q X mean={summary['1q_X_mean']:.4f}, Y mean={summary['1q_Y_mean']:.4f}, Z mean={summary['1q_Z_mean']:.4f}")
    print(f"  2q XX mean={summary['2q_XX_mean']:.4f}, YY mean={summary['2q_YY_mean']:.4f}, ZZ mean={summary['2q_ZZ_mean']:.4f}")

# Key discriminating observables
print("\n\n=== KEY DISCRIMINATING OBSERVABLES ===")
print(f"{'Observable':<20} {'GHZ':>8} {'W':>8} {'Cluster1D':>8} {'Cluster2D':>8}")
print("-" * 56)

# Find observables with highest variance across states
all_keys_1q = set()
all_keys_2q = set()
for name in states:
    all_keys_1q.update(results[name]['profile']['1q'].keys())
    all_keys_2q.update(results[name]['profile']['2q'].keys())

# Check which 2-qubit correlators discriminate best
best_discriminators = []
for key in sorted(all_keys_2q):
    vals = [results[name]['profile']['2q'].get(key, 0) for name in states]
    variance = np.var(vals)
    if variance > 0.01:
        best_discriminators.append((key, variance, vals))

best_discriminators.sort(key=lambda x: -x[1])
for key, var, vals in best_discriminators[:15]:
    print(f"{key:<20} {vals[0]:>8.4f} {vals[1]:>8.4f} {vals[2]:>8.4f} {vals[3]:>8.4f}")

# Save
with open('results/sprint_010a_pauli_profiles.json', 'w') as f:
    json.dump({
        'n_qubits': n,
        'states': list(states.keys()),
        'summaries': {name: results[name]['summary'] for name in states},
        'discriminators': [(k, v, vals) for k, v, vals in best_discriminators[:20]],
    }, f, indent=2)

print("\nResults saved to results/sprint_010a_pauli_profiles.json")
