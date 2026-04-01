"""
Sprint 017a: Concatenated bit-flip code threshold emergence.
Compare uncoded, 1-level (3-qubit), and 2-level (9-qubit) concatenated
bit-flip codes under bit-flip noise. Use classical majority voting
(measure all qubits, post-process).

Threshold for repetition code: p_th = 0.5.
F_L ~ 1 - C*(p/p_th)^(2^L)
"""
import numpy as np
import json
from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, pauli_error

def make_bit_flip_noise(p):
    """Bit-flip on 'id' gates only (storage noise)."""
    noise_model = NoiseModel()
    error_1q = pauli_error([('X', p), ('I', 1 - p)])
    noise_model.add_all_qubit_quantum_error(error_1q, ['id'])
    return noise_model

def majority(bits):
    """Classical majority vote."""
    return 1 if sum(bits) > len(bits) / 2 else 0

def run_experiment(n_qubits, init_state, p, shots=50000):
    """Run repetition code: encode, apply noise, measure all, classical decode."""
    sim = AerSimulator()
    qc = QuantumCircuit(n_qubits, n_qubits)

    if init_state == 1:
        qc.x(0)

    # Encode: CNOT from qubit 0 to all others
    for i in range(1, n_qubits):
        qc.cx(0, i)

    # Storage noise (id on each qubit)
    for i in range(n_qubits):
        qc.id(i)

    # Measure all
    qc.measure(range(n_qubits), range(n_qubits))

    if p > 0:
        noise = make_bit_flip_noise(p)
        result = sim.run(qc, noise_model=noise, shots=shots).result()
    else:
        result = sim.run(qc, shots=shots).result()

    return result.get_counts()

def decode_1level(counts, init_state, n=3):
    """Majority vote on n qubits."""
    correct = 0
    total = 0
    for bitstring, count in counts.items():
        bits = [int(b) for b in bitstring[::-1]]  # reverse for qubit ordering
        decoded = majority(bits[:n])
        if decoded == init_state:
            correct += count
        total += count
    return correct / total

def decode_2level(counts, init_state):
    """2-level concatenation: majority on each block of 3, then majority of block results."""
    correct = 0
    total = 0
    for bitstring, count in counts.items():
        bits = [int(b) for b in bitstring[::-1]]
        # Inner decode: blocks [0,1,2], [3,4,5], [6,7,8]
        block0 = majority(bits[0:3])
        block1 = majority(bits[3:6])
        block2 = majority(bits[6:9])
        # Outer decode
        decoded = majority([block0, block1, block2])
        if decoded == init_state:
            correct += count
        total += count
    return correct / total

# Sweep noise
p_values = [0.0, 0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
shots = 100000

results_data = {"p_values": p_values, "uncoded": {}, "level1_3qubit": {}, "level2_9qubit": {}}

for label, init in [("zero", 0), ("one", 1)]:
    uncoded_fids = []
    level1_fids = []
    level2_fids = []

    for p in p_values:
        print(f"  p={p:.2f}, |{label}>...")

        # Uncoded: single qubit
        counts_1 = run_experiment(1, init, p, shots)
        fid_uncoded = counts_1.get(str(init), 0) / shots
        uncoded_fids.append(round(fid_uncoded, 6))

        # 1-level: 3 qubits
        counts_3 = run_experiment(3, init, p, shots)
        fid_1 = decode_1level(counts_3, init, 3)
        level1_fids.append(round(fid_1, 6))

        # 2-level: 9 qubits
        counts_9 = run_experiment(9, init, p, shots)
        fid_2 = decode_2level(counts_9, init)
        level2_fids.append(round(fid_2, 6))

    results_data["uncoded"][label] = uncoded_fids
    results_data["level1_3qubit"][label] = level1_fids
    results_data["level2_9qubit"][label] = level2_fids

# Theory
p_arr = np.array(p_values)
theory_uncoded = (1 - p_arr).tolist()
theory_1level = (1 - 3*p_arr**2 + 2*p_arr**3).tolist()
p_eff = 3*p_arr**2 - 2*p_arr**3
theory_2level = (1 - 3*p_eff**2 + 2*p_eff**3).tolist()
results_data["theory"] = {"uncoded": theory_uncoded, "level1": theory_1level, "level2": theory_2level}
results_data["theory_threshold"] = 0.5

# Print summary
print("\n=== CONCATENATED BIT-FLIP CODE THRESHOLD ===")
print(f"{'p':>6} | {'Uncoded':>8} | {'1-level':>8} | {'2-level':>8} | {'Th-Unc':>8} | {'Th-L1':>8} | {'Th-L2':>8}")
print("-" * 72)
for i, p in enumerate(p_values):
    avg_u = (results_data["uncoded"]["zero"][i] + results_data["uncoded"]["one"][i]) / 2
    avg_1 = (results_data["level1_3qubit"]["zero"][i] + results_data["level1_3qubit"]["one"][i]) / 2
    avg_2 = (results_data["level2_9qubit"]["zero"][i] + results_data["level2_9qubit"]["one"][i]) / 2
    print(f"{p:6.2f} | {avg_u:8.4f} | {avg_1:8.4f} | {avg_2:8.4f} | {theory_uncoded[i]:8.4f} | {theory_1level[i]:8.4f} | {theory_2level[i]:8.4f}")

# Improvement ratios
print("\n=== IMPROVEMENT OF CODING ===")
print(f"{'p':>6} | {'1L vs unc':>10} | {'2L vs 1L':>10} | {'2L vs unc':>10}")
print("-" * 50)
for i, p in enumerate(p_values):
    avg_u = (results_data["uncoded"]["zero"][i] + results_data["uncoded"]["one"][i]) / 2
    avg_1 = (results_data["level1_3qubit"]["zero"][i] + results_data["level1_3qubit"]["one"][i]) / 2
    avg_2 = (results_data["level2_9qubit"]["zero"][i] + results_data["level2_9qubit"]["one"][i]) / 2
    imp_1u = avg_1 - avg_u
    imp_21 = avg_2 - avg_1
    imp_2u = avg_2 - avg_u
    print(f"{p:6.2f} | {imp_1u:+10.4f} | {imp_21:+10.4f} | {imp_2u:+10.4f}")

with open("results/sprint_017a_concatenated_threshold.json", "w") as f:
    json.dump(results_data, f, indent=2)
print("\nResults saved.")
