"""Sprint 010c: Minimal Pauli classifier for entanglement archetypes.

From 10a we know:
- GHZ: ZZ=1 everywhere, zero 1-body
- W: Z≠0, XX=YY=ZZ=1/3 everywhere
- Cluster 1D: only neighbor XZ≠0
- Cluster 2D: ALL 1-body and 2-body = 0

Question: What is the MINIMUM number of Pauli measurements needed to
classify a state into {GHZ, W, Cluster_1D, Cluster_2D, Product, Other}?

Also test: robustness under noise (can classifier work on noisy states?).
"""

import numpy as np
import json
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, DensityMatrix, SparsePauliOp

def make_ghz(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(n-1):
        qc.cx(i, i+1)
    return Statevector.from_instruction(qc)

def make_w(n):
    state = np.zeros(2**n, dtype=complex)
    for i in range(n):
        state[1 << i] = 1.0 / np.sqrt(n)
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

def make_product(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    return Statevector.from_instruction(qc)

def pauli_exp(sv_or_rho, pauli_str):
    """Compute expectation value of a Pauli string."""
    op = SparsePauliOp.from_list([(pauli_str, 1.0)])
    if isinstance(sv_or_rho, DensityMatrix):
        return np.real(np.trace(sv_or_rho.data @ op.to_matrix()))
    return sv_or_rho.expectation_value(op).real

def build_pauli_label(n, qubit_pauli_dict):
    """Build n-qubit Pauli label with specified qubits."""
    label = ['I'] * n
    for q, p in qubit_pauli_dict.items():
        label[n - 1 - q] = p  # Qiskit little-endian
    return ''.join(label)

def classify(measurements, n):
    """Classify state from minimal Pauli measurements.

    Decision tree:
    1. Measure Z on qubit 0: nonzero → W (only state with <Z>≠0)
    2. Measure ZZ on qubits 0,1: =1 → GHZ (perfect ZZ correlations)
    3. Measure XZ on qubits 0,1: =±1 → Cluster_1D (stabilizer signature)
    4. If all above ≈0 → Cluster_2D or Product
    5. Measure XZZ on qubits 0,1,2 (3-body): nonzero → Cluster_2D
    """
    z0 = measurements['Z0']
    zz01 = measurements['ZZ01']
    xz01 = measurements['XZ01']
    xzz = measurements.get('XZZ_2d', 0.0)

    scores = {
        'W': 0, 'GHZ': 0, 'Cluster_1D': 0, 'Cluster_2D': 0, 'Product': 0
    }

    # Rule 1: W has nonzero <Z>
    if abs(z0) > 0.1:
        scores['W'] += 3

    # Rule 2: GHZ has <ZZ> = 1
    if abs(zz01) > 0.5:
        scores['GHZ'] += 3

    # Rule 3: Cluster 1D has <XZ> = ±1
    if abs(xz01) > 0.5:
        scores['Cluster_1D'] += 3

    # Rule 4: Cluster 2D has nonzero 3-body correlations
    if abs(xzz) > 0.5 and abs(z0) < 0.1 and abs(zz01) < 0.1 and abs(xz01) < 0.1:
        scores['Cluster_2D'] += 3

    # Rule 5: Product has all ≈ 0 at 2-body and 3-body
    if abs(z0) < 0.1 and abs(zz01) < 0.1 and abs(xz01) < 0.1 and abs(xzz) < 0.1:
        scores['Product'] += 2

    best = max(scores, key=scores.get)
    if scores[best] == 0:
        best = 'Unknown'
    return best, scores

def measure_features(state, n):
    """Extract the minimal feature set for classification."""
    feats = {}
    feats['Z0'] = pauli_exp(state, build_pauli_label(n, {0: 'Z'}))
    feats['ZZ01'] = pauli_exp(state, build_pauli_label(n, {0: 'Z', 1: 'Z'}))
    feats['XZ01'] = pauli_exp(state, build_pauli_label(n, {0: 'X', 1: 'Z'}))
    # 3-body: X on q0, Z on neighbors of q0
    # For 1D chain: neighbors of 0 are {1}  → XZ (already captured)
    # For 2D grid (2x3): neighbors of 0 are {1, 3} → X_0 Z_1 Z_3
    feats['XZZ_2d'] = pauli_exp(state, build_pauli_label(n, {0: 'X', 1: 'Z', 3: 'Z'}))
    return feats

n = 6
states = {
    'GHZ': make_ghz(n),
    'W': make_w(n),
    'Cluster_1D': make_cluster_1d(n),
    'Cluster_2D': make_cluster_2d(2, 3),
    'Product': make_product(n),
}

# Part 1: Classification of pure states
print("=== MINIMAL PAULI CLASSIFIER (4 measurements) ===\n")
print(f"{'State':<15} {'Z0':>8} {'ZZ01':>8} {'XZ01':>8} {'XZZ_2d':>8} {'Class':>15}")
print("-" * 65)

pure_results = {}
for name, sv in states.items():
    feats = measure_features(sv, n)
    pred, scores = classify(feats, n)
    correct = pred == name
    pure_results[name] = {
        'features': {k: round(float(v), 6) for k, v in feats.items()},
        'prediction': pred,
        'correct': correct,
    }
    mark = "✓" if correct else f"✗ ({pred})"
    print(f"{name:<15} {feats['Z0']:>8.4f} {feats['ZZ01']:>8.4f} {feats['XZ01']:>8.4f} {feats['XZZ_2d']:>8.4f} {mark:>15}")

accuracy = sum(1 for r in pure_results.values() if r['correct']) / len(pure_results)
print(f"\nPure state accuracy: {accuracy*100:.0f}% ({sum(1 for r in pure_results.values() if r['correct'])}/{len(pure_results)})")

# Part 2: Classification under noise
print("\n\n=== CLASSIFIER ROBUSTNESS UNDER DEPOLARIZING NOISE ===\n")
noise_levels = [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]

noise_results = {}
print(f"{'p':>6}", end='')
for name in states:
    print(f" {name:>12}", end='')
print(f" {'Accuracy':>10}")
print("-" * 80)

for p in noise_levels:
    correct = 0
    noise_results[str(p)] = {}
    print(f"{p:>6.2f}", end='')
    for name, sv in states.items():
        # Apply depolarizing noise
        rho = DensityMatrix(sv)
        noisy = DensityMatrix((1 - p) * rho.data + p * np.eye(2**n) / (2**n))
        feats = measure_features(noisy, n)
        pred, scores = classify(feats, n)
        is_correct = pred == name
        if is_correct:
            correct += 1
        noise_results[str(p)][name] = {
            'prediction': pred,
            'correct': is_correct,
            'features': {k: round(float(v), 6) for k, v in feats.items()},
        }
        mark = "✓" if is_correct else f"→{pred[:3]}"
        print(f" {mark:>12}", end='')
    acc = correct / len(states)
    print(f" {acc*100:>9.0f}%")

# Part 3: How many measurements do we REALLY need?
print("\n\n=== MINIMUM MEASUREMENTS NEEDED ===")
print("\nDecision tree (pure states):")
print("  1. Measure <Z> on any qubit:")
print("     - |<Z>| > 0.1 → W state")
print("  2. Measure <ZZ> on any pair:")
print("     - |<ZZ>| > 0.5 → GHZ state")
print("  3. Measure <XZ> on adjacent pair:")
print("     - |<XZ>| > 0.5 → Cluster 1D")
print("  4. Measure <XZZ> on triple:")
print("     - |<XZZ>| > 0.5 → Cluster 2D")
print("  5. All ~ 0 → Product state")
print("\n  → WORST CASE: 4 measurements (4 Pauli settings)")
print("     vs full tomography: 3^6 - 1 = 728 settings")
print(f"     Compression ratio: {728/4:.0f}x")

# Part 4: Stabilizer-based detection for cluster states
print("\n\n=== STABILIZER WITNESS FOR CLUSTER STATES ===")
# Cluster 1D stabilizers: K_i = X_i * prod(Z_j for j in neighbors of i)
print("\n1D Cluster stabilizer expectations:")
sv_c1d = states['Cluster_1D']
stab_1d = []
for i in range(n):
    paulis = {i: 'X'}
    if i > 0:
        paulis[i-1] = 'Z'
    if i < n-1:
        paulis[i+1] = 'Z'
    label = build_pauli_label(n, paulis)
    val = pauli_exp(sv_c1d, label)
    stab_1d.append(round(float(val), 6))
    print(f"  K_{i} = {label}: <K_{i}> = {val:.4f}")

print(f"  Sum of stabilizers: {sum(stab_1d):.4f} (should be {n} for perfect cluster)")

# Test stabilizer sum on OTHER states
print("\n  Stabilizer sum on other states:")
for name, sv in states.items():
    if name == 'Cluster_1D':
        continue
    stab_sum = 0
    for i in range(n):
        paulis = {i: 'X'}
        if i > 0:
            paulis[i-1] = 'Z'
        if i < n-1:
            paulis[i+1] = 'Z'
        label = build_pauli_label(n, paulis)
        stab_sum += pauli_exp(sv, label)
    print(f"  {name}: {stab_sum:.4f}")

# Save all results
save_data = {
    'n_qubits': n,
    'n_measurements_needed': 4,
    'full_tomography_settings': 728,
    'compression_ratio': 182,
    'pure_state_results': pure_results,
    'noise_robustness': noise_results,
    'stabilizer_sums': {
        'Cluster_1D_on_self': sum(stab_1d),
    },
}

with open('results/sprint_010c_pauli_classifier.json', 'w') as f:
    json.dump(save_data, f, indent=2)

print("\n\nResults saved to results/sprint_010c_pauli_classifier.json")
