"""
Sprint 029a: TFIM Ground State Entanglement Across Phase Transition
Exact diagonalization of H = -J Σ Z_i Z_{i+1} - h Σ X_i
Sweep h/J from 0 to 3, track entropy, MI, I3
"""

import numpy as np
from scipy.sparse import kron, eye, csr_matrix
from scipy.sparse.linalg import eigsh
import json, time

# Pauli matrices (sparse)
I2 = eye(2, format='csr')
X = csr_matrix(np.array([[0, 1], [1, 0]], dtype=complex))
Z = csr_matrix(np.array([[1, 0], [0, -1]], dtype=complex))

def pauli_op(op, qubit, n):
    """Single Pauli operator on qubit `qubit` in n-qubit system."""
    ops = [I2] * n
    ops[qubit] = op
    result = ops[0]
    for o in ops[1:]:
        result = kron(result, o, format='csr')
    return result

def tfim_hamiltonian(n, h_over_J):
    """Build TFIM Hamiltonian: H = -Σ Z_i Z_{i+1} - h Σ X_i (open boundary)."""
    dim = 2**n
    H = csr_matrix((dim, dim), dtype=complex)
    # ZZ interaction (J=1)
    for i in range(n - 1):
        H -= pauli_op(Z, i, n) @ pauli_op(Z, i + 1, n)
    # Transverse field
    for i in range(n):
        H -= h_over_J * pauli_op(X, i, n)
    return H

def partial_trace(rho, keep, n):
    """Partial trace of n-qubit density matrix, keeping qubits in `keep`."""
    d = 2
    rho_full = rho.reshape([d] * 2 * n)
    trace_out = sorted(set(range(n)) - set(keep))
    # Move traced-out indices to the end
    for idx in sorted(trace_out, reverse=True):
        rho_full = np.trace(rho_full, axis1=idx, axis2=idx + n - (n - len(trace_out) - (len(trace_out) - len([t for t in trace_out if t > idx]))))
    # Simpler approach: use einsum-like logic
    # Actually let's use the standard approach
    keep = sorted(keep)
    trace_out = sorted(set(range(n)) - set(keep))

    rho_t = rho.reshape([2]*n + [2]*n)
    # Contract traced-out indices
    for i, q in enumerate(sorted(trace_out, reverse=True)):
        rho_t = np.trace(rho_t, axis1=q, axis2=q + n - i)
        n -= 1  # effectively reduce n for subsequent traces

    k = len(keep)
    return rho_t.reshape(2**k, 2**k)

def partial_trace_simple(state_vec, keep, n):
    """Partial trace from state vector, keeping qubits in `keep`."""
    keep = sorted(keep)
    trace_out = sorted(set(range(n)) - set(keep))

    psi = state_vec.reshape([2]*n)
    rho = np.outer(state_vec, state_vec.conj()).reshape([2]*n + [2]*n)

    # Trace out qubits one by one from the highest index
    offset = 0
    for q in sorted(trace_out, reverse=True):
        rho = np.trace(rho, axis1=q, axis2=q + n - offset)
        offset += 1

    k = len(keep)
    return rho.reshape(2**k, 2**k)

def von_neumann_entropy(rho):
    """Von Neumann entropy of density matrix."""
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return -np.sum(evals * np.log2(evals))

def mutual_information(state_vec, i, j, n):
    """MI between qubits i and j."""
    S_i = von_neumann_entropy(partial_trace_simple(state_vec, [i], n))
    S_j = von_neumann_entropy(partial_trace_simple(state_vec, [j], n))
    S_ij = von_neumann_entropy(partial_trace_simple(state_vec, [i, j], n))
    return S_i + S_j - S_ij

def tripartite_info(state_vec, i, j, k, n):
    """I3(i:j:k) = I(i:j) + I(i:k) - I(i:jk)."""
    S_i = von_neumann_entropy(partial_trace_simple(state_vec, [i], n))
    S_j = von_neumann_entropy(partial_trace_simple(state_vec, [j], n))
    S_k = von_neumann_entropy(partial_trace_simple(state_vec, [k], n))
    S_ij = von_neumann_entropy(partial_trace_simple(state_vec, [i, j], n))
    S_ik = von_neumann_entropy(partial_trace_simple(state_vec, [i, k], n))
    S_jk = von_neumann_entropy(partial_trace_simple(state_vec, [j, k], n))
    S_ijk = von_neumann_entropy(partial_trace_simple(state_vec, [i, j, k], n))

    I3 = (S_i + S_j + S_k - S_ij - S_ik - S_jk + S_ijk)
    return I3

def half_cut_entropy(state_vec, n):
    """Entropy of half-cut bipartition."""
    keep = list(range(n // 2))
    rho_half = partial_trace_simple(state_vec, keep, n)
    return von_neumann_entropy(rho_half)

# Main sweep
results = {}
h_values = np.linspace(0.01, 3.0, 25)  # avoid h=0 degeneracy

for n in [4, 6, 8]:
    print(f"\n=== n={n} qubits ===")
    t0 = time.time()

    n_results = {
        'h_values': [],
        'half_cut_entropy': [],
        'energy': [],
        'energy_gap': [],
        'avg_MI': [],
        'max_MI': [],
        'nearest_MI': [],
        'avg_I3': [],
        'min_I3': [],
        'consecutive_I3': [],
    }

    for h in h_values:
        H = tfim_hamiltonian(n, h)
        # Get ground state and first excited state
        evals, evecs = eigsh(H, k=2, which='SA')
        gs_energy = evals[0]
        gap = evals[1] - evals[0]
        psi = evecs[:, 0]

        # Half-cut entropy
        S_half = half_cut_entropy(psi, n)

        # Pairwise MI (all pairs)
        mi_values = []
        nn_mi = []  # nearest-neighbor only
        for i in range(n):
            for j in range(i+1, n):
                mi = mutual_information(psi, i, j, n)
                mi_values.append(mi)
                if j == i + 1:
                    nn_mi.append(mi)

        # I3 for selected triples (consecutive and some non-consecutive)
        i3_values = []
        consec_i3 = []
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if n <= 6 or (k - i <= 3):  # limit for n=8
                        i3 = tripartite_info(psi, i, j, k, n)
                        i3_values.append(i3)
                        if j == i+1 and k == j+1:
                            consec_i3.append(i3)

        n_results['h_values'].append(float(h))
        n_results['half_cut_entropy'].append(float(S_half))
        n_results['energy'].append(float(gs_energy))
        n_results['energy_gap'].append(float(gap))
        n_results['avg_MI'].append(float(np.mean(mi_values)))
        n_results['max_MI'].append(float(np.max(mi_values)))
        n_results['nearest_MI'].append(float(np.mean(nn_mi)))
        n_results['avg_I3'].append(float(np.mean(i3_values)))
        n_results['min_I3'].append(float(np.min(i3_values)))
        n_results['consecutive_I3'].append(float(np.mean(consec_i3)) if consec_i3 else 0.0)

        print(f"  h/J={h:.2f}: E={gs_energy:.4f}, gap={gap:.4f}, S_half={S_half:.4f}, "
              f"avg_MI={np.mean(mi_values):.4f}, min_I3={np.min(i3_values):.4f}")

    elapsed = time.time() - t0
    print(f"  n={n} completed in {elapsed:.1f}s")
    results[f'n={n}'] = n_results

# Save results
with open('results/sprint_029a_tfim_ground_state.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\n=== ANALYSIS ===")
for n_key in results:
    data = results[n_key]
    h_vals = data['h_values']
    S_vals = data['half_cut_entropy']
    gap_vals = data['energy_gap']

    # Find entropy peak
    peak_idx = np.argmax(S_vals)
    print(f"\n{n_key}:")
    print(f"  Entropy peak at h/J = {h_vals[peak_idx]:.3f}, S = {S_vals[peak_idx]:.4f}")
    print(f"  Gap minimum at h/J = {h_vals[np.argmin(gap_vals)]:.3f}, gap = {min(gap_vals):.4f}")
    print(f"  Ordered phase (h=0.01): MI={data['avg_MI'][0]:.4f}, I3={data['avg_I3'][0]:.4f}")
    print(f"  Critical (peak): MI={data['avg_MI'][peak_idx]:.4f}, I3={data['avg_I3'][peak_idx]:.4f}")
    print(f"  Disordered (h=3.0): MI={data['avg_MI'][-1]:.4f}, I3={data['avg_I3'][-1]:.4f}")

print("\nResults saved to results/sprint_029a_tfim_ground_state.json")
