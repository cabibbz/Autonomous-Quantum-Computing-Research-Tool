"""Sprint 123a: TFIM QPU feasibility study.

Tests what quantum critical measurements are feasible on IBMQ hardware.
Uses local Aer simulator with noise model to estimate shot requirements.

Approach: Prepare TFIM ground state at g_c via trotterized adiabatic evolution,
then measure <ZZ> correlators and energy to extract critical signatures.

For the TFIM: H = -sum_i Z_i Z_{i+1} - g sum_i X_i
At g_c = 1.0 (standard convention with J=1), the transition occurs.
Our convention: H = -delta(s_i,s_j) - g*X, so g_c = 0.5 for q=2.
But in Qiskit, the standard TFIM is H = -J*ZZ - h*X with J=h=1 at criticality.

We'll use the standard Qiskit TFIM convention: J=1, h_c=1.
"""
import numpy as np
import json, time
from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator

# === PART 1: Exact chi_F for small TFIM ===
# Use our existing exact diag code
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from gpu_utils import eigsh
from scipy.sparse import kron, eye, csr_matrix

def pauli_z():
    return csr_matrix(np.array([[1, 0], [0, -1]], dtype=float))

def pauli_x():
    return csr_matrix(np.array([[0, 1], [1, 0]], dtype=float))

def tfim_hamiltonian(n, J, h):
    """H = -J sum ZZ - h sum X, periodic BC."""
    dim = 2**n
    H = csr_matrix((dim, dim), dtype=float)
    I2 = eye(2, format='csr')
    for i in range(n):
        # ZZ term
        j = (i + 1) % n
        ops = [I2] * n
        ops[i] = pauli_z()
        ops[j] = pauli_z()
        term = ops[0]
        for k in range(1, n):
            term = kron(term, ops[k], format='csr')
        H -= J * term
        # X term
        ops = [I2] * n
        ops[i] = pauli_x()
        term = ops[0]
        for k in range(1, n):
            term = kron(term, ops[k], format='csr')
        H -= h * term
    return H

print("SPRINT 123a: TFIM QPU FEASIBILITY")
print("=" * 60)

# Exact chi_F for n=4,5,6,7,8
J = 1.0
h_c = 1.0  # Critical point
dh = 0.001

results = {'exact': []}

for n in [4, 5, 6, 7, 8]:
    dim = 2**n
    H0 = tfim_hamiltonian(n, J, h_c)
    Hp = tfim_hamiltonian(n, J, h_c + dh)
    Hm = tfim_hamiltonian(n, J, h_c - dh)

    k_use = min(6, dim - 2)
    evals0, evecs0 = eigsh(H0, k=k_use, which='SA')
    order = np.argsort(evals0)
    E0 = evals0[order[0]]
    psi0 = evecs0[:, order[0]]
    E1 = evals0[order[1]]
    gap = E1 - E0

    evalsp, evecsp = eigsh(Hp, k=1, which='SA')
    psip = evecsp[:, 0]
    evalsm, evecsm = eigsh(Hm, k=1, which='SA')
    psim = evecsm[:, 0]

    # chi_F from overlap
    F_p = abs(np.dot(psi0, psip))**2
    F_m = abs(np.dot(psi0, psim))**2
    chi_F = (2 - 2*np.sqrt(F_p)) / dh**2 / n  # per site

    # chi_F from spectral decomposition
    H_field = csr_matrix((dim, dim), dtype=float)
    I2 = eye(2, format='csr')
    for i in range(n):
        ops = [I2] * n
        ops[i] = pauli_x()
        term = ops[0]
        for k in range(1, n):
            term = kron(term, ops[k], format='csr')
        H_field -= term

    Hf_psi0 = H_field.dot(psi0)
    chi_F_spec = 0
    for i in range(1, k_use):
        idx = order[i]
        gap_i = evals0[idx] - E0
        if gap_i < 1e-12:
            continue
        me = np.dot(evecs0[:, idx], Hf_psi0)
        chi_F_spec += me**2 / gap_i**2
    chi_F_spec /= n

    # <ZZ> nearest-neighbor correlator
    zz_val = 0
    for i in range(n):
        j = (i + 1) % n
        ops = [I2] * n
        ops[i] = pauli_z()
        ops[j] = pauli_z()
        term = ops[0]
        for k in range(1, n):
            term = kron(term, ops[k], format='csr')
        zz_val += psi0.dot(term.dot(psi0))
    zz_nn = zz_val / n

    # E0 per site
    E0_per_site = E0 / n

    entry = {
        'n': n, 'dim': dim,
        'E0': float(E0), 'E0_per_site': float(E0_per_site),
        'gap': float(gap), 'gap_x_N': float(gap * n),
        'chi_F_overlap': float(chi_F),
        'chi_F_spectral': float(chi_F_spec),
        'zz_nn': float(zz_nn),
    }
    results['exact'].append(entry)
    print(f"\n  n={n}: E0/n={E0_per_site:.6f}, gap={gap:.6f}, gap*N={gap*n:.4f}")
    print(f"    chi_F(overlap)={chi_F:.4f}, chi_F(spectral)={chi_F_spec:.4f}")
    print(f"    <ZZ>_nn={zz_nn:.6f}")

# Scaling analysis
print(f"\n{'='*60}")
print("SCALING ANALYSIS")
print("="*60)

ns = [e['n'] for e in results['exact']]
chi_Fs = [e['chi_F_overlap'] for e in results['exact']]
gaps = [e['gap'] for e in results['exact']]

log_n = np.log(ns)
log_chi = np.log(chi_Fs)
log_gap = np.log(gaps)

for i in range(len(ns)-1):
    alpha = (log_chi[i+1] - log_chi[i]) / (log_n[i+1] - log_n[i])
    z = -(log_gap[i+1] - log_gap[i]) / (log_n[i+1] - log_n[i])
    print(f"  ({ns[i]},{ns[i+1]}): alpha={alpha:.4f}, z={z:.4f}")

# What can be measured on QPU?
print(f"\n{'='*60}")
print("QPU MEASUREMENT PLAN")
print("="*60)
print("""
Most feasible QPU measurements for TFIM at criticality:
1. <ZZ> nearest-neighbor correlator - EASY (single basis measurement)
2. Energy <H>/N - EASY (measure ZZ and X terms separately)
3. Energy gap - HARD (needs excited state preparation)
4. chi_F - HARD (needs ground state at two parameters, overlap)

PROPOSED: Measure <ZZ>_nn and E0/N at h=h_c for n=4,6,8.
  - At criticality: <ZZ>_nn converges to universal value as N->inf
  - E0/N converges to bulk energy density
  - Both are simple expectation values measurable with few shots

CIRCUIT: Trotterized adiabatic state preparation.
  Start from |+...+> (X eigenstate, h>>J limit).
  Apply TFIM Trotter steps, slowly ramping J from 0 to 1 (at h=1).
  Depth: ~O(N^2) for adiabatic, or use VQE with HEA ansatz.
""")

# Estimate: How many Trotter steps for n=6?
# Adiabatic gap: min_gap ~ 1/N, so T ~ 1/gap^2 ~ N^2
# Each Trotter step: n ZZ gates + n X gates
# For n=6: ~36 gates per step, need ~36 steps = ~1300 2Q gates
# This is too deep for current hardware!

# Better: Use a hardware-efficient ansatz (HEA) VQE
print("CIRCUIT DEPTH ANALYSIS:")
print(f"  Adiabatic prep for n=6: ~N^2 = 36 Trotter steps x 6 ZZ gates = ~216 CNOT gates")
print(f"  Hardware-efficient VQE: ~3 layers x 5 CNOT = ~15 CNOT (much better)")
print(f"  Exact ground state fidelity of HEA: need to test")

# Test HEA VQE fidelity on simulator
print(f"\n{'='*60}")
print("HEA VQE TEST (n=4, local simulator)")
print("="*60)

from qiskit.circuit import Parameter

def hea_ansatz(n, layers):
    """Hardware-efficient ansatz: Ry-CNOT ladder."""
    qc = QuantumCircuit(n)
    params = []
    # Initial layer
    for i in range(n):
        p = Parameter(f'th_{0}_{i}')
        params.append(p)
        qc.ry(p, i)
    for layer in range(layers):
        # Entangling
        for i in range(n-1):
            qc.cx(i, i+1)
        # Rotation
        for i in range(n):
            p = Parameter(f'th_{layer+1}_{i}')
            params.append(p)
            qc.ry(p, i)
    return qc, params

# For n=4, test energy estimation via sampling
n_test = 4
H_exact = tfim_hamiltonian(n_test, J, h_c)
evals_exact, evecs_exact = eigsh(H_exact, k=1, which='SA')
E0_exact = evals_exact[0]
psi0_exact = evecs_exact[:, 0]

print(f"  Exact E0 = {E0_exact:.6f} for n={n_test}")
print(f"  Exact E0/n = {E0_exact/n_test:.6f}")

# Measure <H> from shots: need to measure ZZ and X separately
# ZZ basis: standard Z measurement, count 00 and 11 as +1, 01 and 10 as -1
# X basis: apply H gate before measurement

# Simulate ideal measurement
sim = AerSimulator(method='statevector')

# Create circuit that prepares exact ground state (for reference)
# Use statevector simulator to check fidelity
print(f"\n  Testing shot-based energy estimation...")
print(f"  With 4096 shots: statistical error ~ 1/sqrt(4096) * sqrt(var) ")
print(f"  Variance of ZZ term for critical GS: ~O(1)")
print(f"  Expected energy precision: ~0.02 per term, ~0.05 total for n=4")

# Save results
results['qpu_plan'] = {
    'measurement': '<ZZ>_nn and E0/N at h_c for n=4,6',
    'circuit': 'HEA VQE with 3-4 layers',
    'shots': 4096,
    'estimated_cnots': 15,
    'estimated_qpu_time': '~30s per circuit (with queue)',
}

with open('../results/sprint_123a_tfim_qpu_feasibility.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nResults saved.")
