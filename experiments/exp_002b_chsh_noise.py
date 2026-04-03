"""
Experiment 2b: CHSH Under Depolarizing Noise
How much noise does it take to kill quantum advantage?
2 qubits, statevector + noise model. Should run in <30s.
"""

from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
import numpy as np
import json, os, time

RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

t0 = time.time()

shots = 20000  # enough for good statistics

def measure_chsh_noisy(noise_model, shots=20000):
    """Measure CHSH S at optimal angles with given noise model."""
    sim = AerSimulator(noise_model=noise_model) if noise_model else AerSimulator()
    a_angles = [0, np.pi/2]
    b_angles = [np.pi/4, 3*np.pi/4]
    correlations = {}

    for i, a in enumerate(a_angles):
        for j, b in enumerate(b_angles):
            qc = QuantumCircuit(2, 2)
            qc.h(0)
            qc.cx(0, 1)
            qc.ry(-a, 0)
            qc.ry(-b, 1)
            qc.measure([0, 1], [0, 1])
            counts = sim.run(qc, shots=shots).result().get_counts()
            E = sum((1-2*int(bs[1]))*(1-2*int(bs[0]))*c for bs, c in counts.items()) / shots
            correlations[f"a{i}_b{j}"] = E

    return correlations["a0_b0"] - correlations["a0_b1"] + correlations["a1_b0"] + correlations["a1_b1"]

# Sweep depolarizing noise from 0% to 30%
noise_rates = np.linspace(0, 0.30, 31)
results_data = []

for p in noise_rates:
    if p == 0:
        S = measure_chsh_noisy(None, shots)
    else:
        noise_model = NoiseModel()
        # Add depolarizing error to all single-qubit gates
        error_1q = depolarizing_error(p, 1)
        error_2q = depolarizing_error(p, 2)
        noise_model.add_all_qubit_quantum_error(error_1q, ['h', 'ry'])
        noise_model.add_all_qubit_quantum_error(error_2q, ['cx'])
        S = measure_chsh_noisy(noise_model, shots)

    results_data.append({"noise_rate": float(p), "S": float(S)})
    marker = " <<< CLASSICAL BOUNDARY" if abs(S) <= 2.0 and len(results_data) > 1 and abs(results_data[-2]["S"]) > 2.0 else ""
    print(f"  p={p:.3f}: S={S:.4f}  |S|={'>' if abs(S)>2 else '<='} 2.0{marker}")

# Find the critical noise rate where |S| crosses 2.0
S_values = [d["S"] for d in results_data]
critical_p = None
for i in range(1, len(S_values)):
    if abs(S_values[i-1]) > 2.0 and abs(S_values[i]) <= 2.0:
        # Linear interpolation
        p1, p2 = noise_rates[i-1], noise_rates[i]
        s1, s2 = abs(S_values[i-1]), abs(S_values[i])
        critical_p = p1 + (2.0 - s1) * (p2 - p1) / (s2 - s1)
        break

# Theoretical: for depolarizing channel on both qubits,
# S_noisy = S_ideal * (1-p)^k where k depends on circuit depth
# For our circuit: H, CX, RY, RY = 4 gates, so roughly S ~ 2sqrt(2) * (1-p)^3

dt = time.time() - t0

results = {
    "experiment": "chsh_noise_sensitivity",
    "noise_type": "depolarizing",
    "noise_on": "all gates (1q and 2q)",
    "shots": shots,
    "sweep": results_data,
    "critical_noise_rate": float(critical_p) if critical_p else None,
    "S_at_zero_noise": float(S_values[0]),
    "S_at_max_noise": float(S_values[-1]),
    "runtime_seconds": dt,
}

with open(os.path.join(RESULTS_DIR, "sprint_002b_chsh_noise.json"), "w") as f:
    json.dump(results, f, indent=2)

print(f"\nRuntime: {dt:.1f}s")
print(f"S at zero noise: {S_values[0]:.4f}")
print(f"Critical noise rate: {critical_p:.4f}" if critical_p else "Critical point not found in range")
print(f"S at max noise (p=0.30): {S_values[-1]:.4f}")
print(f"Saved to results/sprint_002b_chsh_noise.json")
