"""Sprint 123c: QPU-ready TFIM circuit — adiabatic state prep + energy measurement.

Strategy:
- Prepare |+...+> (paramagnetic ground state at h>>J)
- Trotterize adiabatic evolution from h=5,J=0 to h=1,J=1 (critical point)
- Measure <ZZ> and <X> to get energy components
- Test on Aer with noise model first

Key insight: Don't try to measure chi_F on QPU (requires overlap).
Instead, measure E(g)/N across the transition and show the QPU can
resolve the phase transition location.
"""
import numpy as np
import json, time
from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error

# === Build circuits ===

def tfim_trotter_step(qc, n, J, h, dt):
    """One Trotter step of exp(-i*H*dt) for H = -J*ZZ - h*X."""
    # ZZ interaction: exp(i*J*dt*Z_i*Z_{i+1}) for each pair
    for i in range(n-1):
        qc.cx(i, i+1)
        qc.rz(2*J*dt, i+1)
        qc.cx(i, i+1)
    # Periodic boundary
    qc.cx(n-1, 0)
    qc.rz(2*J*dt, 0)
    qc.cx(n-1, 0)
    # X field: exp(i*h*dt*X_i) = Rx(2*h*dt)
    for i in range(n):
        qc.rx(2*h*dt, i)

def build_adiabatic_circuit(n, n_steps, J_final=1.0, h_final=1.0):
    """Adiabatic prep from paramagnetic (J=0,h=5) to target (J_final, h_final)."""
    qc = QuantumCircuit(n)
    # Start in |+...+>
    for i in range(n):
        qc.h(i)

    dt = 1.0 / n_steps  # total time = 1
    for step in range(n_steps):
        s = (step + 1) / n_steps  # ramp from 0 to 1
        J = J_final * s
        h = 5.0 * (1 - s) + h_final * s  # ramp from h=5 to h_final
        tfim_trotter_step(qc, n, J, h, dt)

    return qc

def measure_zz(qc_base, n):
    """Measure in Z basis to get <ZZ> correlators."""
    qc = qc_base.copy()
    qc.measure_all()
    return qc

def measure_x(qc_base, n):
    """Measure in X basis to get <X> magnetization."""
    qc = qc_base.copy()
    for i in range(n):
        qc.h(i)
    qc.measure_all()
    return qc

def extract_zz_from_counts(counts, n, shots):
    """Extract <ZZ>_nn from Z-basis measurement counts."""
    zz_sum = 0
    total = 0
    for bitstring, count in counts.items():
        bits = [int(b) for b in bitstring[::-1]]  # reverse for qiskit convention
        for i in range(n):
            j = (i + 1) % n
            # ZZ = +1 if same, -1 if different
            zz_sum += count * (1 if bits[i] == bits[j] else -1)
        total += count
    return zz_sum / (total * n)

def extract_x_from_counts(counts, n, shots):
    """Extract <X>/n from X-basis measurement counts."""
    x_sum = 0
    total = 0
    for bitstring, count in counts.items():
        bits = [int(b) for b in bitstring[::-1]]
        for i in range(n):
            x_sum += count * (1 - 2*bits[i])  # X eigenvalue: 0->+1, 1->-1
        total += count
    return x_sum / (total * n)

# === Test parameters ===
n_qubits = 4
shots = 8192

print("SPRINT 123c: TFIM QPU CIRCUIT DESIGN")
print("=" * 60)

# Test different Trotter step counts
print(f"\n--- Noiseless simulator: n={n_qubits} ---")
sim_ideal = AerSimulator(method='statevector')

# Exact values for reference (from 123a)
exact_zz = {4: 0.653281, 6: 0.643951}
exact_x = {4: 0.653281, 6: 0.643951}

results = {'circuits': []}

for n_steps in [5, 10, 20, 30]:
    qc = build_adiabatic_circuit(n_qubits, n_steps, J_final=1.0, h_final=1.0)
    n_cnots = qc.count_ops().get('cx', 0)
    depth = qc.depth()

    qc_z = measure_zz(qc, n_qubits)
    qc_x = measure_x(qc, n_qubits)

    # Run ideal
    job_z = sim_ideal.run(qc_z, shots=shots)
    job_x = sim_ideal.run(qc_x, shots=shots)
    counts_z = job_z.result().get_counts()
    counts_x = job_x.result().get_counts()

    zz_meas = extract_zz_from_counts(counts_z, n_qubits, shots)
    x_meas = extract_x_from_counts(counts_x, n_qubits, shots)
    energy = -zz_meas - x_meas  # E/N = -J*<ZZ>/N - h*<X>/N at J=h=1

    zz_err = abs(zz_meas - exact_zz[n_qubits])
    x_err = abs(x_meas - exact_x[n_qubits])

    print(f"\n  steps={n_steps:3d}: CNOTs={n_cnots:4d}, depth={depth:4d}")
    print(f"    <ZZ>_nn = {zz_meas:.4f} (exact: {exact_zz[n_qubits]:.4f}, err: {zz_err:.4f})")
    print(f"    <X>/n   = {x_meas:.4f} (exact: {exact_x[n_qubits]:.4f}, err: {x_err:.4f})")
    print(f"    E/N     = {energy:.4f} (exact: {-exact_zz[n_qubits]-exact_x[n_qubits]:.4f})")

    results['circuits'].append({
        'n': n_qubits, 'steps': n_steps, 'cnots': n_cnots, 'depth': depth,
        'zz_ideal': float(zz_meas), 'x_ideal': float(x_meas),
        'zz_err': float(zz_err), 'x_err': float(x_err),
    })

# Now test with noise model
print(f"\n{'='*60}")
print("NOISY SIMULATION (depolarizing noise)")
print("="*60)

# Typical IBM hardware: ~0.3% single-qubit error, ~1% two-qubit error
noise = NoiseModel()
noise.add_all_qubit_quantum_error(depolarizing_error(0.003, 1), ['rx', 'ry', 'rz', 'h'])
noise.add_all_qubit_quantum_error(depolarizing_error(0.01, 2), ['cx'])

sim_noisy = AerSimulator(noise_model=noise)

for n_steps in [5, 10, 20]:
    qc = build_adiabatic_circuit(n_qubits, n_steps, J_final=1.0, h_final=1.0)
    n_cnots = qc.count_ops().get('cx', 0)

    qc_z = measure_zz(qc, n_qubits)
    qc_x = measure_x(qc, n_qubits)

    job_z = sim_noisy.run(qc_z, shots=shots)
    job_x = sim_noisy.run(qc_x, shots=shots)
    counts_z = job_z.result().get_counts()
    counts_x = job_x.result().get_counts()

    zz_meas = extract_zz_from_counts(counts_z, n_qubits, shots)
    x_meas = extract_x_from_counts(counts_x, n_qubits, shots)
    energy = -zz_meas - x_meas

    # Signal survival: fraction of ideal signal retained
    zz_signal = zz_meas / exact_zz[n_qubits]
    x_signal = x_meas / exact_x[n_qubits]

    print(f"\n  steps={n_steps:3d}: CNOTs={n_cnots:4d}")
    print(f"    <ZZ>_nn = {zz_meas:.4f} (signal: {zz_signal:.1%})")
    print(f"    <X>/n   = {x_meas:.4f} (signal: {x_signal:.1%})")
    print(f"    E/N     = {energy:.4f}")

# Phase transition scan (with noise)
print(f"\n{'='*60}")
print("NOISY PHASE TRANSITION SCAN: n=4, 10 Trotter steps")
print("="*60)

h_values = [0.5, 0.7, 0.85, 1.0, 1.15, 1.3, 1.5, 2.0]
n_steps_scan = 10

print(f"\n  {'h':>5s}  {'<ZZ>':>8s}  {'<X>/n':>8s}  {'E/N':>8s}")
scan_results = []
for h in h_values:
    qc = build_adiabatic_circuit(n_qubits, n_steps_scan, J_final=1.0, h_final=h)

    qc_z = measure_zz(qc, n_qubits)
    qc_x = measure_x(qc, n_qubits)

    job_z = sim_noisy.run(qc_z, shots=shots)
    job_x = sim_noisy.run(qc_x, shots=shots)
    counts_z = job_z.result().get_counts()
    counts_x = job_x.result().get_counts()

    zz = extract_zz_from_counts(counts_z, n_qubits, shots)
    x = extract_x_from_counts(counts_x, n_qubits, shots)
    energy = -zz - h * x

    print(f"  {h:5.2f}  {zz:8.4f}  {x:8.4f}  {energy:8.4f}")
    scan_results.append({'h': h, 'zz': float(zz), 'x': float(x), 'energy': float(energy)})

results['phase_scan_noisy'] = scan_results

# Summary and QPU recommendation
print(f"\n{'='*60}")
print("QPU RECOMMENDATION")
print("="*60)
print("""
RECOMMENDED QPU EXPERIMENT:
  - Model: TFIM, n=4 qubits, periodic BC
  - Circuit: Adiabatic prep with 10 Trotter steps (~60 CNOTs)
  - Measurements: <ZZ> (Z basis) and <X> (X basis) at 8 field values
  - Shots: 8192 per circuit
  - Total circuits: 16 (8 h-values x 2 bases)
  - Estimated QPU time: ~30s (well within 580s budget)

  GOAL: Resolve the quantum phase transition location from hardware data.
  Compare QPU energy curve E(h) with exact curve.

  STRETCH: Repeat at n=6 (10 Trotter steps ~ 100 CNOTs, marginal fidelity).
""")

with open('../results/sprint_123c_qpu_circuit.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("Results saved.")
