"""Sprint 010b: Fidelity-based GME witnesses.

For each archetype, construct a witness W = alpha*I - |psi><psi|
where alpha is the max overlap with any biseparable state.
If <W> < 0, the state is genuine multipartite entangled.

Also test: stabilizer witnesses for cluster states using stabilizer generators.
"""

import numpy as np
import json
from itertools import combinations
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, DensityMatrix, partial_trace, state_fidelity

def make_ghz(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(n-1):
        qc.cx(i, i+1)
    return Statevector.from_instruction(qc)

def make_w(n):
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

def biseparable_bound(sv, n):
    """Compute max fidelity of |psi> with any biseparable state.

    For a biseparable state across cut A|B:
    max_bisep fidelity = max eigenvalue of rho_A tensor max eigenvalue of rho_B
    Actually, max fidelity with product states across A|B = max singular value of
    the coefficient matrix squared = max Schmidt coefficient squared.

    For biseparable (not necessarily product), it's the max Schmidt coefficient.
    """
    max_fid = 0
    # Check all bipartitions
    for k in range(1, n // 2 + 1):
        for subsys_a in combinations(range(n), k):
            subsys_b = [q for q in range(n) if q not in subsys_a]
            # Get reduced density matrix of subsystem A
            rho = DensityMatrix(sv)
            rho_a = partial_trace(rho, subsys_b)
            # Max eigenvalue of rho_a = max Schmidt coefficient squared
            eigs = np.real(np.linalg.eigvalsh(rho_a.data))
            max_schmidt_sq = np.max(eigs)
            if max_schmidt_sq > max_fid:
                max_fid = max_schmidt_sq
    return max_fid

def apply_depolarizing(sv, p, n):
    """Apply global depolarizing noise."""
    rho = DensityMatrix(sv)
    noisy = (1 - p) * rho.data + p * np.eye(2**n) / (2**n)
    return DensityMatrix(noisy)

n = 6
states = {
    'GHZ': make_ghz(n),
    'W': make_w(n),
    'Cluster_1D': make_cluster_1d(n),
    'Cluster_2D': make_cluster_2d(2, 3),
}

# Part 1: Compute biseparable bounds and witness values
print("=== FIDELITY-BASED GME WITNESSES ===\n")
print(f"{'State':<15} {'Bisep Bound':>12} {'Fidelity':>10} {'Witness':>10} {'GME?':>6}")
print("-" * 55)

results = {}
for name, sv in states.items():
    alpha = biseparable_bound(sv, n)
    fid = 1.0  # fidelity with itself
    witness_val = alpha - fid  # W = alpha*I - |psi><psi|, <W> = alpha - F
    is_gme = witness_val < 0

    results[name] = {
        'bisep_bound': round(float(alpha), 6),
        'fidelity': fid,
        'witness_value': round(float(witness_val), 6),
        'is_gme': bool(is_gme),
        'gme_margin': round(float(-witness_val), 6),  # how far below 0
    }
    print(f"{name:<15} {alpha:>12.6f} {fid:>10.4f} {witness_val:>10.6f} {'YES' if is_gme else 'NO':>6}")

# Part 2: Witness value under depolarizing noise — when does GME detection fail?
print("\n\n=== GME WITNESS UNDER NOISE ===")
noise_levels = np.arange(0, 1.01, 0.05)

noise_results = {}
for name, sv in states.items():
    alpha = results[name]['bisep_bound']
    witness_vals = []
    gme_threshold = None

    for p in noise_levels:
        noisy_rho = apply_depolarizing(sv, p, n)
        # Fidelity with the pure target state
        fid = state_fidelity(noisy_rho, sv)
        w = alpha - fid
        witness_vals.append({'p': round(float(p), 3), 'fidelity': round(float(fid), 6),
                           'witness': round(float(w), 6)})
        if w >= 0 and gme_threshold is None:
            gme_threshold = round(float(p), 3)

    noise_results[name] = {
        'witness_values': witness_vals,
        'gme_death_threshold': gme_threshold,
    }

print(f"\n{'State':<15} {'GME Death p':>12} {'Bisep Bound':>12}")
print("-" * 41)
for name in states:
    thresh = noise_results[name]['gme_death_threshold']
    bound = results[name]['bisep_bound']
    print(f"{name:<15} {thresh if thresh else '>1.0':>12} {bound:>12.6f}")

# Part 3: Cross-fidelity — does each state's witness detect the others?
print("\n\n=== CROSS-FIDELITY MATRIX ===")
print("(Fidelity of row-state with column-target)")
print(f"{'':>15}", end='')
for name in states:
    print(f" {name:>12}", end='')
print()

cross_fid = {}
for name1, sv1 in states.items():
    cross_fid[name1] = {}
    print(f"{name1:>15}", end='')
    for name2, sv2 in states.items():
        f = state_fidelity(sv1, sv2)
        cross_fid[name1][name2] = round(float(f), 6)
        print(f" {f:>12.6f}", end='')
    print()

# Part 4: Can a witness for state A detect GME in state B?
print("\n\n=== WITNESS CROSS-DETECTION ===")
print("(Witness_A applied to state_B: <0 means detected)")
print(f"{'Witness\\State':>15}", end='')
for name in states:
    print(f" {name:>12}", end='')
print()

cross_witness = {}
for name_w in states:
    alpha = results[name_w]['bisep_bound']
    cross_witness[name_w] = {}
    print(f"{name_w:>15}", end='')
    for name_s, sv_s in states.items():
        fid = state_fidelity(sv_s, states[name_w])
        w = alpha - fid
        cross_witness[name_w][name_s] = round(float(w), 6)
        detected = "✓" if w < -0.001 else "✗"
        print(f" {w:>10.4f}{detected}", end='')
    print()

# Save
save_data = {
    'n_qubits': n,
    'gme_witnesses': results,
    'noise_sensitivity': {name: {'gme_death_threshold': noise_results[name]['gme_death_threshold']}
                         for name in states},
    'cross_fidelity': cross_fid,
    'cross_witness': cross_witness,
}

with open('results/sprint_010b_gme_witnesses.json', 'w') as f:
    json.dump(save_data, f, indent=2)

print("\nResults saved to results/sprint_010b_gme_witnesses.json")
