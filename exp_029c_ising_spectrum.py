"""
Sprint 029c: Entanglement Spectrum at Criticality
Full negativity spectrum across all bipartitions.
Compare TFIM ground states to archetype states (GHZ, W, Cluster).
Test: does the critical ground state match any known archetype?
"""

import numpy as np
from scipy.sparse import kron, eye, csr_matrix
from scipy.sparse.linalg import eigsh
from itertools import combinations
import json, time

# Pauli matrices
I2 = eye(2, format='csr')
X_mat = csr_matrix(np.array([[0, 1], [1, 0]], dtype=complex))
Z_mat = csr_matrix(np.array([[1, 0], [0, -1]], dtype=complex))
H_gate = np.array([[1, 1], [1, -1]]) / np.sqrt(2)

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

def negativity(state_vec, subsystem_A, n):
    """Logarithmic negativity of bipartition A vs complement."""
    rho = np.outer(state_vec, state_vec.conj()).reshape([2]*n + [2]*n)
    # Partial transpose over A
    A = sorted(subsystem_A)
    rho_pt = rho.copy()
    for q in A:
        rho_pt = np.swapaxes(rho_pt, q, q + n)
    rho_pt = rho_pt.reshape(2**n, 2**n)
    evals = np.linalg.eigvalsh(rho_pt)
    neg = (np.sum(np.abs(evals)) - 1.0) / 2.0
    return neg

def make_ghz(n):
    """GHZ state |00...0⟩ + |11...1⟩ normalized."""
    psi = np.zeros(2**n, dtype=complex)
    psi[0] = 1.0 / np.sqrt(2)
    psi[-1] = 1.0 / np.sqrt(2)
    return psi

def make_w(n):
    """W state: equal superposition of all single-excitation states."""
    psi = np.zeros(2**n, dtype=complex)
    for i in range(n):
        idx = 1 << (n - 1 - i)
        psi[idx] = 1.0 / np.sqrt(n)
    return psi

def make_cluster_1d(n):
    """1D cluster state: |+⟩^n with CZ between neighbors."""
    # Start with |+⟩^n
    psi = np.ones(2**n, dtype=complex) / np.sqrt(2**n)
    # Apply CZ gates
    for i in range(n - 1):
        for idx in range(2**n):
            bits = [(idx >> (n-1-q)) & 1 for q in range(n)]
            if bits[i] == 1 and bits[i+1] == 1:
                psi[idx] *= -1
    return psi

def entanglement_spectrum(state_vec, n):
    """Compute negativity for ALL bipartitions up to half-size."""
    spectrum = {}
    for size in range(1, n // 2 + 1):
        negs = []
        for subset in combinations(range(n), size):
            neg = negativity(state_vec, list(subset), n)
            negs.append(neg)
        spectrum[f'|A|={size}'] = {
            'values': [float(v) for v in negs],
            'mean': float(np.mean(negs)),
            'std': float(np.std(negs)),
            'min': float(np.min(negs)),
            'max': float(np.max(negs)),
        }
    return spectrum

def pairwise_mi_matrix(state_vec, n):
    """Full pairwise MI matrix."""
    mi = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            S_i = von_neumann_entropy(partial_trace_simple(state_vec, [i], n))
            S_j = von_neumann_entropy(partial_trace_simple(state_vec, [j], n))
            S_ij = von_neumann_entropy(partial_trace_simple(state_vec, [i, j], n))
            mi[i, j] = S_i + S_j - S_ij
            mi[j, i] = mi[i, j]
    return mi

# Main experiment
n = 6
results = {}

# TFIM ground states at three phase points
h_values = {'ordered': 0.3, 'critical': 0.8, 'disordered': 2.0}

print("=== TFIM Ground States ===")
for label, h in h_values.items():
    print(f"\n--- {label} (h/J={h}) ---")
    t0 = time.time()

    H = tfim_hamiltonian(n, h)
    evals, evecs = eigsh(H, k=1, which='SA')
    psi = evecs[:, 0]

    spec = entanglement_spectrum(psi, n)
    mi_mat = pairwise_mi_matrix(psi, n)
    S_half = von_neumann_entropy(partial_trace_simple(psi, list(range(n//2)), n))

    results[f'tfim_{label}'] = {
        'h_over_J': h,
        'energy': float(evals[0]),
        'half_cut_entropy': float(S_half),
        'spectrum': spec,
        'mi_matrix': mi_mat.tolist(),
        'total_MI': float(np.sum(mi_mat) / 2),
        'mi_uniformity': float(np.std(mi_mat[np.triu_indices(n, 1)]) /
                               (np.mean(mi_mat[np.triu_indices(n, 1)]) + 1e-10)),
    }

    for size_key, data in spec.items():
        print(f"  {size_key}: mean={data['mean']:.4f}, std={data['std']:.4f}, "
              f"range=[{data['min']:.4f}, {data['max']:.4f}]")
    print(f"  Half-cut entropy: {S_half:.4f}")
    print(f"  Total MI: {np.sum(mi_mat)/2:.4f}, MI uniformity (CV): "
          f"{np.std(mi_mat[np.triu_indices(n, 1)]) / (np.mean(mi_mat[np.triu_indices(n, 1)]) + 1e-10):.4f}")
    print(f"  Time: {time.time()-t0:.1f}s")

# Archetype states for comparison
print("\n=== Archetype States ===")
archetypes = {
    'GHZ': make_ghz(n),
    'W': make_w(n),
    'Cluster_1D': make_cluster_1d(n),
}

for label, psi in archetypes.items():
    print(f"\n--- {label} ---")
    t0 = time.time()

    spec = entanglement_spectrum(psi, n)
    mi_mat = pairwise_mi_matrix(psi, n)
    S_half = von_neumann_entropy(partial_trace_simple(psi, list(range(n//2)), n))

    results[f'archetype_{label}'] = {
        'half_cut_entropy': float(S_half),
        'spectrum': spec,
        'mi_matrix': mi_mat.tolist(),
        'total_MI': float(np.sum(mi_mat) / 2),
        'mi_uniformity': float(np.std(mi_mat[np.triu_indices(n, 1)]) /
                               (np.mean(mi_mat[np.triu_indices(n, 1)]) + 1e-10)),
    }

    for size_key, data in spec.items():
        print(f"  {size_key}: mean={data['mean']:.4f}, std={data['std']:.4f}, "
              f"range=[{data['min']:.4f}, {data['max']:.4f}]")
    print(f"  Half-cut entropy: {S_half:.4f}")
    print(f"  Total MI: {np.sum(mi_mat)/2:.4f}, MI uniformity (CV): "
          f"{np.std(mi_mat[np.triu_indices(n, 1)]) / (np.mean(mi_mat[np.triu_indices(n, 1)]) + 1e-10):.4f}")
    print(f"  Time: {time.time()-t0:.1f}s")

# Similarity analysis
print("\n=== SIMILARITY ANALYSIS ===")
print("\nComparing TFIM ground states to archetypes (negativity spectrum cosine similarity):")

def spectrum_vector(spec):
    """Flatten spectrum into vector for comparison."""
    vec = []
    for size_key in sorted(spec.keys()):
        vec.extend(spec[size_key]['values'])
    return np.array(vec)

def cosine_sim(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-10 or nb < 1e-10:
        return 0.0
    return float(np.dot(a, b) / (na * nb))

for tfim_label in ['ordered', 'critical', 'disordered']:
    tfim_vec = spectrum_vector(results[f'tfim_{tfim_label}']['spectrum'])
    print(f"\n  TFIM {tfim_label} (h/J={h_values[tfim_label]}):")
    for arch_label in ['GHZ', 'W', 'Cluster_1D']:
        arch_vec = spectrum_vector(results[f'archetype_{arch_label}']['spectrum'])
        sim = cosine_sim(tfim_vec, arch_vec)
        print(f"    vs {arch_label}: {sim:.4f}")

# Save
with open('results/sprint_029c_ising_spectrum.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nResults saved to results/sprint_029c_ising_spectrum.json")
