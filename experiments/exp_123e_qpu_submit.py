"""Sprint 123e: Submit TFIM critical ground state measurement to QPU.

Circuit: exact state preparation (11 CX gates for n=4).
Measurements: <ZZ>_nn and <X>/n at criticality (J=h=1).
Also measure at h=0.5 (ordered) and h=2.0 (paramagnetic) as controls.

QPU justification:
1. Prediction: <ZZ>_nn = 0.653, <X>/n = 0.653 at criticality for n=4 TFIM
2. Surprising if: hardware value deviates >20% from exact (beyond noise model prediction)
3. Action: If <10% deviation, validates quantum critical state prep on hardware.
   If >20%, quantifies the noise floor for critical-point measurements.
"""
import numpy as np
import json, time
from scipy.sparse import kron, eye, csr_matrix
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from gpu_utils import eigsh

from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

def pauli_z():
    return csr_matrix(np.array([[1, 0], [0, -1]], dtype=float))

def pauli_x():
    return csr_matrix(np.array([[0, 1], [1, 0]], dtype=float))

def tfim_hamiltonian(n, J, h):
    dim = 2**n
    H = csr_matrix((dim, dim), dtype=float)
    I2 = eye(2, format='csr')
    for i in range(n):
        j = (i + 1) % n
        ops_zz = [I2] * n
        ops_zz[i] = pauli_z()
        ops_zz[j] = pauli_z()
        term = ops_zz[0]
        for k in range(1, n):
            term = kron(term, ops_zz[k], format='csr')
        H -= J * term
        ops_x = [I2] * n
        ops_x[i] = pauli_x()
        term = ops_x[0]
        for k in range(1, n):
            term = kron(term, ops_x[k], format='csr')
        H -= h * term
    return H

def extract_observables(counts, n):
    """Extract <ZZ>_nn and <X>/n from counts (Z-basis measurement)."""
    zz_sum = 0
    total = 0
    for bitstring, count in counts.items():
        bits = [int(b) for b in bitstring.replace(' ', '')[::-1]]
        if len(bits) < n:
            bits = bits + [0] * (n - len(bits))
        for i in range(n):
            j = (i + 1) % n
            zz_sum += count * (1 if bits[i] == bits[j] else -1)
        total += count
    return zz_sum / (total * n)

def extract_x(counts, n):
    """Extract <X>/n from X-basis counts."""
    x_sum = 0
    total = 0
    for bitstring, count in counts.items():
        bits = [int(b) for b in bitstring.replace(' ', '')[::-1]]
        if len(bits) < n:
            bits = bits + [0] * (n - len(bits))
        for i in range(n):
            x_sum += count * (1 - 2*bits[i])
        total += count
    return x_sum / (total * n)

n = 4
shots = 8192

print("SPRINT 123e: QPU SUBMISSION — TFIM CRITICAL STATE")
print("=" * 60)

# Prepare ground states at 3 coupling values
h_values = [0.5, 1.0, 2.0]  # ordered, critical, paramagnetic
labels = ['ordered', 'critical', 'paramagnetic']

circuits_z = []  # Z-basis measurement circuits
circuits_x = []  # X-basis measurement circuits
exact_vals = []

for h, label in zip(h_values, labels):
    H = tfim_hamiltonian(n, 1.0, h)
    evals, evecs = eigsh(H, k=1, which='SA')
    psi0 = evecs[:, 0]
    E0 = evals[0]

    # Build circuit
    qc = QuantumCircuit(n)
    qc.initialize(psi0, range(n))

    # Decompose
    pm = generate_preset_pass_manager(optimization_level=3, basis_gates=['cx', 'rz', 'sx', 'x', 'id'])
    qc_dec = pm.run(qc)
    n_cx = qc_dec.count_ops().get('cx', 0)

    # Z-basis measurement
    qc_z = qc_dec.copy()
    qc_z.measure_all()
    circuits_z.append(qc_z)

    # X-basis measurement
    qc_x = qc_dec.copy()
    for i in range(n):
        qc_x.h(i)
    qc_x.measure_all()
    circuits_x.append(qc_x)

    # Exact observables
    I2 = eye(2, format='csr')
    zz_exact = 0
    for i in range(n):
        j = (i + 1) % n
        ops = [I2] * n
        ops[i] = pauli_z()
        ops[j] = pauli_z()
        term = ops[0]
        for k in range(1, n):
            term = kron(term, ops[k], format='csr')
        zz_exact += psi0 @ term.dot(psi0)
    zz_exact /= n

    x_exact = 0
    for i in range(n):
        ops = [I2] * n
        ops[i] = pauli_x()
        term = ops[0]
        for k in range(1, n):
            term = kron(term, ops[k], format='csr')
        x_exact += psi0 @ term.dot(psi0)
    x_exact /= n

    exact_vals.append({'h': h, 'label': label, 'E0_per_n': float(E0/n),
                       'zz_exact': float(zz_exact), 'x_exact': float(x_exact),
                       'n_cx': n_cx})
    print(f"\n  h={h} ({label}): E0/n={E0/n:.6f}, <ZZ>_nn={zz_exact:.6f}, <X>/n={x_exact:.6f}, CX={n_cx}")

# Run on local simulator first
print(f"\n{'='*60}")
print("LOCAL SIMULATION (ideal + noisy)")
print("="*60)

sim_ideal = AerSimulator(method='statevector')
from qiskit_aer.noise import NoiseModel, depolarizing_error
noise = NoiseModel()
noise.add_all_qubit_quantum_error(depolarizing_error(0.003, 1), ['sx', 'x', 'rz', 'h'])
noise.add_all_qubit_quantum_error(depolarizing_error(0.01, 2), ['cx'])
sim_noisy = AerSimulator(noise_model=noise)

results = {'exact': exact_vals, 'ideal_sim': [], 'noisy_sim': []}

for i, (h, label) in enumerate(zip(h_values, labels)):
    # Ideal
    job_z = sim_ideal.run(circuits_z[i], shots=shots)
    job_x = sim_ideal.run(circuits_x[i], shots=shots)
    zz_id = extract_observables(job_z.result().get_counts(), n)
    x_id = extract_x(job_x.result().get_counts(), n)

    # Noisy
    job_z = sim_noisy.run(circuits_z[i], shots=shots)
    job_x = sim_noisy.run(circuits_x[i], shots=shots)
    zz_ny = extract_observables(job_z.result().get_counts(), n)
    x_ny = extract_x(job_x.result().get_counts(), n)

    zz_ex = exact_vals[i]['zz_exact']
    x_ex = exact_vals[i]['x_exact']
    zz_sig = zz_ny / zz_ex if zz_ex != 0 else 0
    x_sig = x_ny / x_ex if x_ex != 0 else 0

    print(f"\n  h={h} ({label}):")
    print(f"    Ideal:  <ZZ>={zz_id:.4f}, <X>/n={x_id:.4f}")
    print(f"    Noisy:  <ZZ>={zz_ny:.4f} ({zz_sig:.0%}), <X>/n={x_ny:.4f} ({x_sig:.0%})")
    print(f"    Exact:  <ZZ>={zz_ex:.4f}, <X>/n={x_ex:.4f}")

    results['ideal_sim'].append({'h': h, 'zz': float(zz_id), 'x': float(x_id)})
    results['noisy_sim'].append({'h': h, 'zz': float(zz_ny), 'x': float(x_ny),
                                  'zz_signal': float(zz_sig), 'x_signal': float(x_sig)})

# QPU submission
print(f"\n{'='*60}")
print("QPU SUBMISSION")
print("="*60)

try:
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler

    service = QiskitRuntimeService()
    backend = service.least_busy(operational=True, simulator=False, min_num_qubits=4)
    print(f"  Backend: {backend.name}")
    print(f"  Qubits: {backend.num_qubits}")

    # Transpile for hardware
    pm_hw = generate_preset_pass_manager(backend=backend, optimization_level=3)
    all_circuits = circuits_z + circuits_x
    all_transpiled = pm_hw.run(all_circuits)

    for i, tc in enumerate(all_transpiled):
        basis = 'Z' if i < 3 else 'X'
        h_idx = i % 3
        n_cx_hw = tc.count_ops().get('cx', 0) + tc.count_ops().get('ecr', 0) + tc.count_ops().get('cz', 0)
        depth_hw = tc.depth()
        print(f"  Circuit {i} (h={h_values[h_idx]}, {basis}-basis): {n_cx_hw} 2Q gates, depth={depth_hw}")

    # Submit
    sampler = Sampler(mode=backend)
    print(f"\n  Submitting {len(all_transpiled)} circuits, {shots} shots each...")
    job = sampler.run(all_transpiled, shots=shots)
    print(f"  Job ID: {job.job_id()}")
    print(f"  Waiting for results...")

    t0 = time.time()
    result = job.result()
    dt = time.time() - t0
    print(f"  Results received in {dt:.1f}s")

    # Extract metrics
    try:
        metrics = job.metrics()
        qpu_seconds = metrics.get('usage', {}).get('quantum_seconds', 'unknown')
        print(f"  QPU time: {qpu_seconds}s")
        results['qpu_time'] = qpu_seconds
    except Exception as e:
        print(f"  Could not get metrics: {e}")

    # Process results
    results['qpu'] = []
    print(f"\n  QPU RESULTS:")
    for i in range(6):
        counts = result[i].data.meas.get_counts()
        basis = 'Z' if i < 3 else 'X'
        h_idx = i % 3
        h = h_values[h_idx]
        label = labels[h_idx]

        if basis == 'Z':
            val = extract_observables(counts, n)
            obs_name = '<ZZ>_nn'
            exact = exact_vals[h_idx]['zz_exact']
        else:
            val = extract_x(counts, n)
            obs_name = '<X>/n'
            exact = exact_vals[h_idx]['x_exact']

        signal = val / exact if exact != 0 else 0
        print(f"    h={h} ({label}), {obs_name}: QPU={val:.4f}, exact={exact:.4f}, signal={signal:.0%}")
        results['qpu'].append({'h': h, 'basis': basis, 'observable': obs_name,
                               'qpu_value': float(val), 'exact_value': float(exact),
                               'signal': float(signal)})

except ImportError:
    print("  qiskit-ibm-runtime not installed. Skipping QPU submission.")
    print("  Install with: pip install qiskit-ibm-runtime")
except Exception as e:
    print(f"  QPU submission failed: {e}")
    print("  Results from noisy simulation saved as fallback.")

with open('../results/sprint_123e_qpu_submit.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved.")
