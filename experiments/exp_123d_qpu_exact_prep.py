"""Sprint 123d: QPU test using exact ground state preparation.

For n=4: use Qiskit's initialize to prepare the exact TFIM ground state,
then measure observables. This gives the minimum-depth circuit for
exact state preparation.

Key question: Is the circuit depth after decomposition feasible on hardware?
"""
import numpy as np
import json, time
from scipy.sparse import kron, eye, csr_matrix
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from gpu_utils import eigsh

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
from qiskit.transpiler import generate_preset_pass_manager

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

print("SPRINT 123d: EXACT STATE PREP QPU TEST")
print("=" * 60)

results = {}

for n in [4]:
    print(f"\n--- n={n} ---")
    H = tfim_hamiltonian(n, 1.0, 1.0)
    evals, evecs = eigsh(H, k=2, which='SA')
    idx0 = np.argmin(evals)
    E0 = evals[idx0]
    psi0 = evecs[:, idx0]

    print(f"  Exact E0/n = {E0/n:.6f}")

    # Prepare circuit with exact state
    qc = QuantumCircuit(n)
    qc.initialize(psi0, range(n))

    # Decompose into basic gates
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    # Use a generic backend with CX basis
    pm = generate_preset_pass_manager(optimization_level=3, basis_gates=['cx', 'rz', 'sx', 'x', 'id'])
    qc_decomposed = pm.run(qc)

    n_cx = qc_decomposed.count_ops().get('cx', 0)
    depth = qc_decomposed.depth()
    print(f"  Decomposed circuit: {n_cx} CX gates, depth={depth}")
    print(f"  Gate counts: {dict(qc_decomposed.count_ops())}")

    # Test with ideal simulator
    sim_ideal = AerSimulator(method='statevector')
    shots = 8192

    # Z basis measurement
    qc_z = qc_decomposed.copy()
    qc_z.measure_all()
    job_z = sim_ideal.run(qc_z, shots=shots)
    counts_z = job_z.result().get_counts()

    zz_sum = 0
    total = 0
    for bitstring, count in counts_z.items():
        bits = [int(b) for b in bitstring.replace(' ', '')[::-1]]
        if len(bits) < n:
            bits = bits + [0] * (n - len(bits))
        for i in range(n):
            j = (i + 1) % n
            zz_sum += count * (1 if bits[i] == bits[j] else -1)
        total += count
    zz_ideal = zz_sum / (total * n)

    # X basis measurement
    qc_x = qc_decomposed.copy()
    for i in range(n):
        qc_x.h(i)
    qc_x.measure_all()
    job_x = sim_ideal.run(qc_x, shots=shots)
    counts_x = job_x.result().get_counts()

    x_sum = 0
    total = 0
    for bitstring, count in counts_x.items():
        bits = [int(b) for b in bitstring.replace(' ', '')[::-1]]
        if len(bits) < n:
            bits = bits + [0] * (n - len(bits))
        for i in range(n):
            x_sum += count * (1 - 2*bits[i])
        total += count
    x_ideal = x_sum / (total * n)

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

    print(f"\n  Ideal sim: <ZZ>={zz_ideal:.4f} (exact {zz_exact:.4f}), <X>/n={x_ideal:.4f} (exact {x_exact:.4f})")

    # Test with noise
    print(f"\n  Noisy simulation (1% CX error, 0.3% 1Q error):")
    noise = NoiseModel()
    noise.add_all_qubit_quantum_error(depolarizing_error(0.003, 1), ['sx', 'x', 'rz', 'h'])
    noise.add_all_qubit_quantum_error(depolarizing_error(0.01, 2), ['cx'])
    sim_noisy = AerSimulator(noise_model=noise)

    job_z = sim_noisy.run(qc_z, shots=shots)
    counts_z = job_z.result().get_counts()
    zz_sum = 0; total = 0
    for bitstring, count in counts_z.items():
        bits = [int(b) for b in bitstring.replace(' ', '')[::-1]]
        if len(bits) < n:
            bits = bits + [0] * (n - len(bits))
        for i in range(n):
            j = (i + 1) % n
            zz_sum += count * (1 if bits[i] == bits[j] else -1)
        total += count
    zz_noisy = zz_sum / (total * n)

    job_x = sim_noisy.run(qc_x, shots=shots)
    counts_x = job_x.result().get_counts()
    x_sum = 0; total = 0
    for bitstring, count in counts_x.items():
        bits = [int(b) for b in bitstring.replace(' ', '')[::-1]]
        if len(bits) < n:
            bits = bits + [0] * (n - len(bits))
        for i in range(n):
            x_sum += count * (1 - 2*bits[i])
        total += count
    x_noisy = x_sum / (total * n)

    zz_signal = zz_noisy / zz_exact
    x_signal = x_noisy / x_exact
    energy_noisy = -zz_noisy - x_noisy
    energy_exact = -zz_exact - x_exact

    print(f"  <ZZ>={zz_noisy:.4f} (signal: {zz_signal:.1%})")
    print(f"  <X>/n={x_noisy:.4f} (signal: {x_signal:.1%})")
    print(f"  E/N={energy_noisy:.4f} (exact: {energy_exact:.4f})")

    # Estimate with 2% CX error (more realistic for some QPUs)
    print(f"\n  Noisy simulation (2% CX error):")
    noise2 = NoiseModel()
    noise2.add_all_qubit_quantum_error(depolarizing_error(0.005, 1), ['sx', 'x', 'rz', 'h'])
    noise2.add_all_qubit_quantum_error(depolarizing_error(0.02, 2), ['cx'])
    sim_noisy2 = AerSimulator(noise_model=noise2)

    job_z = sim_noisy2.run(qc_z, shots=shots)
    counts_z = job_z.result().get_counts()
    zz_sum = 0; total = 0
    for bitstring, count in counts_z.items():
        bits = [int(b) for b in bitstring.replace(' ', '')[::-1]]
        if len(bits) < n:
            bits = bits + [0] * (n - len(bits))
        for i in range(n):
            j = (i + 1) % n
            zz_sum += count * (1 if bits[i] == bits[j] else -1)
        total += count
    zz_noisy2 = zz_sum / (total * n)

    job_x = sim_noisy2.run(qc_x, shots=shots)
    counts_x = job_x.result().get_counts()
    x_sum = 0; total = 0
    for bitstring, count in counts_x.items():
        bits = [int(b) for b in bitstring.replace(' ', '')[::-1]]
        if len(bits) < n:
            bits = bits + [0] * (n - len(bits))
        for i in range(n):
            x_sum += count * (1 - 2*bits[i])
        total += count
    x_noisy2 = x_sum / (total * n)

    print(f"  <ZZ>={zz_noisy2:.4f} (signal: {zz_noisy2/zz_exact:.1%})")
    print(f"  <X>/n={x_noisy2:.4f} (signal: {x_noisy2/x_exact:.1%})")

    results[f'n{n}'] = {
        'n_cx': n_cx, 'depth': depth,
        'zz_exact': float(zz_exact), 'x_exact': float(x_exact),
        'zz_ideal': float(zz_ideal), 'x_ideal': float(x_ideal),
        'zz_noisy_1pct': float(zz_noisy), 'x_noisy_1pct': float(x_noisy),
        'zz_signal_1pct': float(zz_signal), 'x_signal_1pct': float(x_signal),
        'zz_noisy_2pct': float(zz_noisy2), 'x_noisy_2pct': float(x_noisy2),
    }

# Decision
print(f"\n{'='*60}")
print("QPU DECISION")
print("="*60)
n_cx_4 = results['n4']['n_cx']
sig_1 = results['n4']['zz_signal_1pct']
print(f"""
Circuit depth for n=4: {n_cx_4} CX gates
Signal at 1% CX error: {sig_1:.1%} of ideal <ZZ>

Assessment:
- If CX count < 30: QPU measurement is feasible with error mitigation
- If CX count 30-100: Marginal, need zero-noise extrapolation
- If CX count > 100: Not feasible on current hardware

The initialize decomposition for n=4 requires {n_cx_4} CX gates.
""")

if n_cx_4 <= 50:
    print("VERDICT: PROCEED with QPU submission")
    print("  Circuit is shallow enough for meaningful hardware measurement.")
else:
    print("VERDICT: QPU not productive for this experiment")
    print("  Circuit too deep. Consider:")
    print("  - Free-fermion circuit (JW transform, Givens rotations)")
    print("  - Shallower approximate state prep")
    print("  - Different observable (not requiring full GS prep)")

with open('../results/sprint_123d_qpu_exact_prep.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved.")
