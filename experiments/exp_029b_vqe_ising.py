"""
Sprint 029b: VQE for TFIM Ground State — Entanglement During Optimization
Track how entanglement measures evolve as VQE converges to ground state.
Compare h/J = 0.5 (ordered), 1.0 (critical), 2.0 (disordered).
"""

import numpy as np
from scipy.sparse import kron, eye, csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.optimize import minimize
from qiskit.circuit import QuantumCircuit, Parameter
from qiskit_aer import AerSimulator
import json, time

# Pauli matrices
I2 = eye(2, format='csr')
X_mat = csr_matrix(np.array([[0, 1], [1, 0]], dtype=complex))
Z_mat = csr_matrix(np.array([[1, 0], [0, -1]], dtype=complex))

def pauli_op(op, qubit, n):
    ops = [I2] * n
    ops[qubit] = op
    result = ops[0]
    for o in ops[1:]:
        result = kron(result, o, format='csr')
    return result

def tfim_hamiltonian(n, h_over_J):
    dim = 2**n
    H = csr_matrix((dim, dim), dtype=complex)
    for i in range(n - 1):
        H -= pauli_op(Z_mat, i, n) @ pauli_op(Z_mat, i + 1, n)
    for i in range(n):
        H -= h_over_J * pauli_op(X_mat, i, n)
    return H

def partial_trace_simple(state_vec, keep, n):
    keep = sorted(keep)
    trace_out = sorted(set(range(n)) - set(keep))
    rho = np.outer(state_vec, state_vec.conj()).reshape([2]*n + [2]*n)
    offset = 0
    for q in sorted(trace_out, reverse=True):
        rho = np.trace(rho, axis1=q, axis2=q + n - offset)
        offset += 1
    k = len(keep)
    return rho.reshape(2**k, 2**k)

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return -np.sum(evals * np.log2(evals))

def mutual_information(state_vec, i, j, n):
    S_i = von_neumann_entropy(partial_trace_simple(state_vec, [i], n))
    S_j = von_neumann_entropy(partial_trace_simple(state_vec, [j], n))
    S_ij = von_neumann_entropy(partial_trace_simple(state_vec, [i, j], n))
    return S_i + S_j - S_ij

def tripartite_info(state_vec, i, j, k, n):
    S_i = von_neumann_entropy(partial_trace_simple(state_vec, [i], n))
    S_j = von_neumann_entropy(partial_trace_simple(state_vec, [j], n))
    S_k = von_neumann_entropy(partial_trace_simple(state_vec, [k], n))
    S_ij = von_neumann_entropy(partial_trace_simple(state_vec, [i, j], n))
    S_ik = von_neumann_entropy(partial_trace_simple(state_vec, [i, k], n))
    S_jk = von_neumann_entropy(partial_trace_simple(state_vec, [j, k], n))
    S_ijk = von_neumann_entropy(partial_trace_simple(state_vec, [i, j, k], n))
    return S_i + S_j + S_k - S_ij - S_ik - S_jk + S_ijk

def half_cut_entropy(state_vec, n):
    keep = list(range(n // 2))
    return von_neumann_entropy(partial_trace_simple(state_vec, keep, n))

# VQE with Hamiltonian Variational Ansatz (HVA)
# Each layer: RZZ on nearest-neighbors + RX on all qubits
def hva_circuit(n, n_layers, params):
    """Hamiltonian Variational Ansatz for TFIM."""
    qc = QuantumCircuit(n)
    # Initial state: |+>^n (good for disordered phase, neutral overall)
    for i in range(n):
        qc.h(i)

    idx = 0
    for layer in range(n_layers):
        # ZZ interaction layer
        for i in range(n - 1):
            qc.cx(i, i + 1)
            qc.rz(params[idx], i + 1)
            qc.cx(i, i + 1)
            idx += 1
        # X field layer
        for i in range(n):
            qc.rx(params[idx], i)
            idx += 1
    return qc

def get_statevector(qc):
    """Get statevector from circuit using Aer."""
    sim = AerSimulator(method='statevector')
    qc_copy = qc.copy()
    qc_copy.save_statevector()
    result = sim.run(qc_copy).result()
    return np.array(result.get_statevector())

def energy_expectation(state_vec, H):
    """Compute <psi|H|psi>."""
    return np.real(state_vec.conj() @ H @ state_vec)

# Main experiment
n = 6  # 6 qubits
n_layers = 4  # 4 HVA layers
n_params_per_layer = (n - 1) + n  # 5 ZZ + 6 RX = 11
n_params = n_layers * n_params_per_layer

h_values = [0.5, 1.0, 2.0]  # ordered, critical, disordered

results = {}

for h in h_values:
    print(f"\n=== h/J = {h} ===")
    t0 = time.time()

    H = tfim_hamiltonian(n, h)
    # Exact ground state energy
    E_exact = eigsh(H, k=1, which='SA', return_eigenvectors=False)[0]
    print(f"  Exact E = {E_exact:.6f}")

    # VQE optimization with trajectory tracking
    trajectory = {
        'iteration': [],
        'energy': [],
        'energy_error': [],
        'half_cut_entropy': [],
        'avg_MI': [],
        'avg_I3': [],
        'min_I3': [],
    }

    eval_count = [0]

    def objective(params):
        qc = hva_circuit(n, n_layers, params)
        psi = get_statevector(qc)
        E = energy_expectation(psi, H.toarray())

        # Track every 5th evaluation
        eval_count[0] += 1
        if eval_count[0] % 5 == 1 or eval_count[0] <= 3:
            S_half = half_cut_entropy(psi, n)

            mi_vals = []
            for i in range(n):
                for j in range(i+1, n):
                    mi_vals.append(mutual_information(psi, i, j, n))

            i3_vals = []
            for i in range(n):
                for j in range(i+1, n):
                    for k in range(j+1, n):
                        i3_vals.append(tripartite_info(psi, i, j, k, n))

            trajectory['iteration'].append(eval_count[0])
            trajectory['energy'].append(float(E))
            trajectory['energy_error'].append(float(abs(E - E_exact)))
            trajectory['half_cut_entropy'].append(float(S_half))
            trajectory['avg_MI'].append(float(np.mean(mi_vals)))
            trajectory['avg_I3'].append(float(np.mean(i3_vals)))
            trajectory['min_I3'].append(float(np.min(i3_vals)))

        return E

    # Random initial parameters, small values
    np.random.seed(42)
    x0 = np.random.uniform(-0.1, 0.1, n_params)

    # Optimize with COBYLA (gradient-free, robust)
    result = minimize(objective, x0, method='COBYLA',
                     options={'maxiter': 200, 'rhobeg': 0.5})

    # Final state analysis
    qc_final = hva_circuit(n, n_layers, result.x)
    psi_final = get_statevector(qc_final)
    E_final = energy_expectation(psi_final, H.toarray())
    S_final = half_cut_entropy(psi_final, n)

    mi_final = []
    for i in range(n):
        for j in range(i+1, n):
            mi_final.append(mutual_information(psi_final, i, j, n))

    i3_final = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                i3_final.append(tripartite_info(psi_final, i, j, k, n))

    elapsed = time.time() - t0

    results[f'h={h}'] = {
        'h_over_J': h,
        'exact_energy': float(E_exact),
        'vqe_energy': float(E_final),
        'energy_error': float(abs(E_final - E_exact)),
        'relative_error': float(abs(E_final - E_exact) / abs(E_exact)),
        'final_entropy': float(S_final),
        'final_avg_MI': float(np.mean(mi_final)),
        'final_avg_I3': float(np.mean(i3_final)),
        'final_min_I3': float(np.min(i3_final)),
        'n_evals': eval_count[0],
        'trajectory': trajectory,
        'elapsed_s': elapsed,
    }

    print(f"  VQE E = {E_final:.6f} (error: {abs(E_final - E_exact):.6f}, "
          f"rel: {abs(E_final - E_exact)/abs(E_exact)*100:.3f}%)")
    print(f"  Entropy: {S_final:.4f}, avg MI: {np.mean(mi_final):.4f}, "
          f"avg I3: {np.mean(i3_final):.4f}, min I3: {np.min(i3_final):.4f}")
    print(f"  {eval_count[0]} evaluations in {elapsed:.1f}s")

    # Reset counter for next run
    eval_count[0] = 0

# Save results
with open('results/sprint_029b_vqe_ising.json', 'w') as f:
    json.dump(results, f, indent=2)

# Analysis
print("\n=== ANALYSIS ===")
print("\nComparison with exact (from 29a for n=6):")
# Exact values from 29a
exact_ref = {
    0.5: {'S': 1.003, 'MI': 0.812, 'I3': 0.729},
    1.0: {'S': 0.467, 'MI': 0.244, 'I3': 0.029},
    2.0: {'S': 0.127, 'MI': 0.068, 'I3': -0.007},
}

for h in h_values:
    r = results[f'h={h}']
    ex = exact_ref[h]
    print(f"\nh/J = {h}:")
    print(f"  Energy: VQE={r['vqe_energy']:.4f}, exact={r['exact_energy']:.4f}, "
          f"error={r['relative_error']*100:.3f}%")
    print(f"  Entropy: VQE={r['final_entropy']:.4f}, exact≈{ex['S']:.3f}")
    print(f"  avg MI: VQE={r['final_avg_MI']:.4f}, exact≈{ex['MI']:.3f}")
    print(f"  avg I3: VQE={r['final_avg_I3']:.4f}, exact≈{ex['I3']:.3f}")

    # How many trajectory points?
    traj = r['trajectory']
    if traj['iteration']:
        print(f"  Trajectory: {len(traj['iteration'])} points, "
              f"entropy {traj['half_cut_entropy'][0]:.4f} → {traj['half_cut_entropy'][-1]:.4f}")

print("\nResults saved to results/sprint_029b_vqe_ising.json")
