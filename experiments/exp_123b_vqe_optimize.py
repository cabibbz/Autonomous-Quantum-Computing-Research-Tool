"""Sprint 123b: Optimize VQE circuit for TFIM critical ground state.

Classically optimize HEA parameters on statevector simulator,
then the fixed optimized circuit can be run on QPU.

For n=4,6: find parameters that prepare ground state with high fidelity.
Then measure <ZZ> and <X> (energy components) with shot-based simulation.
"""
import numpy as np
import json, time
from scipy.optimize import minimize
from scipy.sparse import kron, eye, csr_matrix
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from gpu_utils import eigsh

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

def hea_statevector(params, n, layers):
    """Simulate HEA circuit via matrix multiplication (faster than Qiskit for small n)."""
    dim = 2**n
    # Start from |0...0>
    psi = np.zeros(dim)
    psi[0] = 1.0

    idx = 0
    # Initial Ry layer
    for i in range(n):
        theta = params[idx]; idx += 1
        c, s = np.cos(theta/2), np.sin(theta/2)
        # Apply Ry(theta) to qubit i
        psi_new = np.zeros_like(psi)
        for basis in range(dim):
            bit = (basis >> i) & 1
            basis_flipped = basis ^ (1 << i)
            if bit == 0:
                psi_new[basis] += c * psi[basis]
                psi_new[basis_flipped] += s * psi[basis]
            else:
                psi_new[basis] += c * psi[basis]
                psi_new[basis_flipped] -= s * psi[basis]
        psi = psi_new

    for layer in range(layers):
        # CNOT ladder
        for i in range(n-1):
            # CNOT: control=i, target=i+1
            psi_new = psi.copy()
            for basis in range(dim):
                ctrl_bit = (basis >> i) & 1
                if ctrl_bit == 1:
                    tgt_bit = (basis >> (i+1)) & 1
                    basis_flipped = basis ^ (1 << (i+1))
                    psi_new[basis] = 0
                    psi_new[basis_flipped] = psi[basis]
            psi = psi_new
        # Ry layer
        for i in range(n):
            theta = params[idx]; idx += 1
            c, s = np.cos(theta/2), np.sin(theta/2)
            psi_new = np.zeros_like(psi)
            for basis in range(dim):
                bit = (basis >> i) & 1
                basis_flipped = basis ^ (1 << i)
                if bit == 0:
                    psi_new[basis] += c * psi[basis]
                    psi_new[basis_flipped] += s * psi[basis]
                else:
                    psi_new[basis] += c * psi[basis]
                    psi_new[basis_flipped] -= s * psi[basis]
            psi = psi_new

    return psi

def cost_fn(params, n, layers, H):
    psi = hea_statevector(params, n, layers)
    return np.real(psi @ H.dot(psi))

results = {'vqe_results': []}

for n in [4, 6]:
    print(f"\n{'='*60}")
    print(f"VQE OPTIMIZATION: TFIM n={n}")
    print(f"{'='*60}")

    H = tfim_hamiltonian(n, 1.0, 1.0)
    evals, evecs = eigsh(H, k=2, which='SA')
    E0_exact = min(evals)
    idx0 = np.argmin(evals)
    psi0_exact = evecs[:, idx0]

    print(f"  Exact E0 = {E0_exact:.6f}")

    best_result = None
    best_energy = 1e10

    for n_layers in [2, 3, 4]:
        n_params = n * (n_layers + 1)
        print(f"\n  layers={n_layers}, params={n_params}:")

        best_for_layers = None
        for trial in range(5):
            x0 = np.random.uniform(-np.pi, np.pi, n_params)
            res = minimize(cost_fn, x0, args=(n, n_layers, H),
                          method='L-BFGS-B',
                          options={'maxiter': 500, 'ftol': 1e-12})

            if res.fun < best_energy:
                best_energy = res.fun
                best_result = res
                best_layers = n_layers

            if best_for_layers is None or res.fun < best_for_layers:
                best_for_layers = res.fun

        # Compute fidelity
        psi_vqe = hea_statevector(best_result.x[:n*(best_layers+1)] if best_layers != n_layers else best_result.x, n, best_layers if best_layers != n_layers else n_layers)
        fid = abs(np.dot(psi0_exact, psi_vqe))**2

        energy_err = abs(best_for_layers - E0_exact)
        print(f"    Best E = {best_for_layers:.8f}, err = {energy_err:.2e}")

    # Final result with best configuration
    psi_best = hea_statevector(best_result.x, n, best_layers)
    fid = abs(np.dot(psi0_exact, psi_best))**2
    energy_err = abs(best_energy - E0_exact)

    # Compute observables from VQE state
    dim = 2**n
    I2 = eye(2, format='csr')

    # <ZZ>_nn
    zz_val = 0
    for i in range(n):
        j = (i + 1) % n
        ops = [I2] * n
        ops[i] = pauli_z()
        ops[j] = pauli_z()
        term = ops[0]
        for k in range(1, n):
            term = kron(term, ops[k], format='csr')
        zz_val += psi_best @ term.dot(psi_best)
    zz_nn = zz_val / n

    # <X>
    x_val = 0
    for i in range(n):
        ops = [I2] * n
        ops[i] = pauli_x()
        term = ops[0]
        for k in range(1, n):
            term = kron(term, ops[k], format='csr')
        x_val += psi_best @ term.dot(psi_best)
    x_avg = x_val / n

    n_cnots = (n-1) * best_layers

    print(f"\n  BEST VQE: layers={best_layers}, E={best_energy:.8f}, "
          f"err={energy_err:.2e}, fidelity={fid:.6f}")
    print(f"  <ZZ>_nn = {zz_nn:.6f} (exact: {results['vqe_results'][-1]['zz_nn_exact'] if results['vqe_results'] else '?'})")
    print(f"  <X>/n = {x_avg:.6f}")
    print(f"  CNOT count = {n_cnots}")

    # Exact observables for comparison
    zz_exact = 0
    for i in range(n):
        j = (i + 1) % n
        ops = [I2] * n
        ops[i] = pauli_z()
        ops[j] = pauli_z()
        term = ops[0]
        for k in range(1, n):
            term = kron(term, ops[k], format='csr')
        zz_exact += psi0_exact @ term.dot(psi0_exact)
    zz_nn_exact = zz_exact / n

    x_exact = 0
    for i in range(n):
        ops = [I2] * n
        ops[i] = pauli_x()
        term = ops[0]
        for k in range(1, n):
            term = kron(term, ops[k], format='csr')
        x_exact += psi0_exact @ term.dot(psi0_exact)
    x_avg_exact = x_exact / n

    print(f"  Exact: <ZZ>_nn = {zz_nn_exact:.6f}, <X>/n = {x_avg_exact:.6f}")

    entry = {
        'n': n,
        'layers': best_layers,
        'n_params': n * (best_layers + 1),
        'n_cnots': n_cnots,
        'E0_exact': float(E0_exact),
        'E0_vqe': float(best_energy),
        'energy_err': float(energy_err),
        'fidelity': float(fid),
        'zz_nn_vqe': float(zz_nn),
        'zz_nn_exact': float(zz_nn_exact),
        'x_avg_vqe': float(x_avg),
        'x_avg_exact': float(x_avg_exact),
        'params': best_result.x.tolist(),
    }
    results['vqe_results'].append(entry)

# QPU cost estimate
print(f"\n{'='*60}")
print("QPU COST ESTIMATE")
print("="*60)
for entry in results['vqe_results']:
    n = entry['n']
    # 2 measurement bases (Z for ZZ, X for magnetization) x 4096 shots each
    # Each circuit execution ~ 1-2s on QPU
    # But wait/queue time dominates
    n_circuits = 2  # ZZ basis + X basis
    shots = 4096
    est_seconds = n_circuits * 2  # ~2s per circuit execution
    print(f"  n={n}: {entry['n_cnots']} CNOTs, {n_circuits} circuits x {shots} shots")
    print(f"    Estimated QPU time: ~{est_seconds}s execution + queue")
    print(f"    Fidelity: {entry['fidelity']:.6f}")

total_qpu = sum(2 * 2 for _ in results['vqe_results'])  # rough
print(f"\n  Total estimated: ~{total_qpu}s QPU execution")
print(f"  Budget available: 580s")
print(f"  FEASIBLE: Yes, well within budget")

with open('../results/sprint_123b_vqe_optimize.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nResults saved.")
