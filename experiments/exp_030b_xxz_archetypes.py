"""
Sprint 030b: XXZ Archetype Classification
Map each phase to the 5 known archetypes: Democratic (GHZ), Distributed (W),
Geometric (Cluster), Topological (Toric), Scale-Free (Critical).

Key improvements over 30a:
- Use Sz=0 sector for FM phase (unique ground state = proper GHZ superposition)
- Compute archetype distances using MI pattern, I3 sign, negativity uniformity
- Compare to TFIM reference values from Sprint 029
"""

import numpy as np
from scipy.sparse import kron, eye, csr_matrix
from scipy.sparse.linalg import eigsh
import json, time

# Pauli matrices
I2 = eye(2, format='csr')
Sp = csr_matrix(np.array([[0, 1], [0, 0]], dtype=complex))
Sm = csr_matrix(np.array([[0, 0], [1, 0]], dtype=complex))
Sz_op = csr_matrix(np.array([[0.5, 0], [0, -0.5]], dtype=complex))

def chain_op(op, qubit, n):
    ops = [I2] * n
    ops[qubit] = op
    result = ops[0]
    for o in ops[1:]:
        result = kron(result, o, format='csr')
    return result

def xxz_hamiltonian(n, delta):
    dim = 2**n
    H = csr_matrix((dim, dim), dtype=complex)
    for i in range(n - 1):
        H += 0.5 * (chain_op(Sp, i, n) @ chain_op(Sm, i+1, n) +
                     chain_op(Sm, i, n) @ chain_op(Sp, i+1, n))
        H += delta * chain_op(Sz_op, i, n) @ chain_op(Sz_op, i+1, n)
    return H

def sz_total_operator(n):
    """Total Sz operator."""
    Sz_total = csr_matrix((2**n, 2**n), dtype=complex)
    for i in range(n):
        Sz_total += chain_op(Sz_op, i, n)
    return Sz_total

def project_to_sz_sector(H, n, sz_target):
    """Project Hamiltonian to Sz=sz_target sector.
    Returns projected H and the projector basis vectors."""
    dim = 2**n
    # Find basis states with Sz = sz_target
    basis_indices = []
    for i in range(dim):
        bits = format(i, f'0{n}b')
        sz = sum(0.5 if b == '0' else -0.5 for b in bits)
        if abs(sz - sz_target) < 1e-10:
            basis_indices.append(i)

    # Build projector
    P = np.zeros((dim, len(basis_indices)))
    for j, idx in enumerate(basis_indices):
        P[idx, j] = 1.0

    # Project H
    H_dense = H.toarray()
    H_sector = P.T @ H_dense @ P
    return H_sector, P, basis_indices

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

def negativity(state_vec, subsys_a, n):
    rho = np.outer(state_vec, state_vec.conj())
    rho_reshaped = rho.reshape([2]*n + [2]*n)
    perm = list(range(2*n))
    for q in subsys_a:
        perm[q], perm[q+n] = perm[q+n], perm[q]
    rho_pt = np.transpose(rho_reshaped, perm).reshape(2**n, 2**n)
    evals = np.linalg.eigvalsh(rho_pt)
    neg = (np.sum(np.abs(evals)) - 1.0) / 2.0
    return max(0.0, neg)

def compute_archetype_features(psi, n):
    """Compute the feature vector for archetype classification."""
    # All pairwise MI
    mi_values = []
    nn_mi = []
    for i in range(n):
        for j in range(i+1, n):
            mi = mutual_information(psi, i, j, n)
            mi_values.append(mi)
            if j == i + 1:
                nn_mi.append(mi)

    # All I3 (n=6: 20 triples)
    i3_values = []
    consec_i3 = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                i3 = tripartite_info(psi, i, j, k, n)
                i3_values.append(i3)
                if j == i+1 and k == j+1:
                    consec_i3.append(i3)

    # Negativity spectrum (all bipartitions of size 1,2,3)
    neg_values = []
    for size in [1, 2, 3]:
        from itertools import combinations
        for subset in combinations(range(n), size):
            neg = negativity(psi, list(subset), n)
            neg_values.append(neg)

    mi_arr = np.array(mi_values)
    i3_arr = np.array(i3_values)
    neg_arr = np.array(neg_values)

    features = {
        'avg_MI': float(np.mean(mi_arr)),
        'mi_uniformity': float(np.std(mi_arr) / np.mean(mi_arr)) if np.mean(mi_arr) > 1e-10 else float('inf'),
        'nn_mi_fraction': float(np.sum(nn_mi) / np.sum(mi_arr)) if np.sum(mi_arr) > 1e-10 else 0.0,
        'avg_I3': float(np.mean(i3_arr)),
        'i3_sign': 'positive' if np.mean(i3_arr) > 0.05 else ('negative' if np.mean(i3_arr) < -0.05 else 'mixed'),
        'fraction_negative_I3': float(np.mean(i3_arr < -0.01)),
        'avg_negativity': float(np.mean(neg_arr)),
        'neg_uniformity': float(np.std(neg_arr) / np.mean(neg_arr)) if np.mean(neg_arr) > 1e-10 else float('inf'),
        'half_cut_entropy': float(von_neumann_entropy(partial_trace_simple(psi, list(range(n//2)), n))),
    }
    return features

# Known archetype reference signatures (from Sprints 004-010, 029)
# Format: (avg_MI, MI_CV, I3_sign, nn_mi_fraction)
ARCHETYPES = {
    'Democratic (GHZ)':   {'avg_MI': 1.0, 'mi_cv': 0.0, 'avg_I3': 1.0, 'nn_frac': 0.33},
    'Distributed (W)':    {'avg_MI': 0.38, 'mi_cv': 0.0, 'avg_I3': 0.20, 'nn_frac': 0.33},
    'Geometric (Cluster)':{'avg_MI': 0.13, 'mi_cv': 1.5, 'avg_I3': -0.4, 'nn_frac': 1.0},
    'Scale-Free (Crit)':  {'avg_MI': 0.45, 'mi_cv': 0.4, 'avg_I3': 0.0, 'nn_frac': 0.5},
}

def archetype_distance(features, arch):
    """Simple normalized distance to archetype."""
    d = 0
    d += ((features['avg_MI'] - arch['avg_MI']) / max(arch['avg_MI'], 0.1))**2
    d += ((features['mi_uniformity'] - arch['mi_cv']) / max(arch['mi_cv'], 0.1))**2
    d += ((features['avg_I3'] - arch['avg_I3']) / max(abs(arch['avg_I3']), 0.1))**2
    d += ((features['nn_mi_fraction'] - arch['nn_frac']) / max(arch['nn_frac'], 0.1))**2
    return np.sqrt(d)

# ---- Main computation ----
n = 6
print(f"=== XXZ Archetype Classification, n={n} ===\n")

# Representative Δ values for each phase + transitions
test_points = {
    'FM deep (Δ=-2)':     -2.0,
    'FM near (Δ=-1.2)':   -1.2,
    'FM transition (Δ=-1)': -1.0,
    'XY near-FM (Δ=-0.5)': -0.5,
    'XY isotropic (Δ=0)': 0.0,
    'XY near-BKT (Δ=0.5)': 0.5,
    'XY at Δ=0.9':        0.9,
    'BKT transition (Δ=1)': 1.0,
    'Néel near (Δ=1.5)':  1.5,
    'Néel deep (Δ=3)':    3.0,
}

results = {}

for label, delta in test_points.items():
    print(f"--- {label} ---")
    t0 = time.time()

    H = xxz_hamiltonian(n, delta)

    # For FM phase, use Sz=0 sector to get unique ground state
    if delta <= -1.0:
        print("  Using Sz=0 sector projection")
        H_sector, P, basis_idx = project_to_sz_sector(H, n, sz_target=0.0)
        evals, evecs = np.linalg.eigh(H_sector)
        gs_energy = evals[0]
        gap = evals[1] - evals[0] if len(evals) > 1 else 0.0
        # Reconstruct full-space state vector
        psi = P @ evecs[:, 0]
    else:
        evals, evecs = eigsh(H, k=3, which='SA')
        gs_energy = evals[0]
        gap = evals[1] - evals[0]
        psi = evecs[:, 0]

    # Make psi real if possible (avoid phase issues)
    if np.max(np.abs(psi.imag)) < 1e-10:
        psi = psi.real
    psi = psi / np.linalg.norm(psi)

    print(f"  E = {gs_energy:.6f}, gap = {gap:.6f}")

    # Compute archetype features
    features = compute_archetype_features(psi, n)
    features['energy'] = float(gs_energy)
    features['gap'] = float(gap)
    features['delta'] = float(delta)

    print(f"  S_half = {features['half_cut_entropy']:.4f}")
    print(f"  avg_MI = {features['avg_MI']:.4f}, MI_CV = {features['mi_uniformity']:.4f}")
    print(f"  avg_I3 = {features['avg_I3']:.4f} ({features['i3_sign']})")
    print(f"  nn_MI_frac = {features['nn_mi_fraction']:.4f}")
    print(f"  neg = {features['avg_negativity']:.4f}, neg_CV = {features['neg_uniformity']:.4f}")

    # Classify
    distances = {}
    for arch_name, arch_ref in ARCHETYPES.items():
        dist = archetype_distance(features, arch_ref)
        distances[arch_name] = float(dist)

    best = min(distances, key=distances.get)
    features['archetype_distances'] = distances
    features['best_archetype'] = best

    print(f"  Distances: {', '.join(f'{k.split()[0]}={v:.3f}' for k,v in distances.items())}")
    print(f"  >>> BEST MATCH: {best} (d={distances[best]:.3f})")

    results[label] = features
    print(f"  ({time.time()-t0:.1f}s)")

# Save
with open('results/sprint_030b_xxz_archetypes.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

# Summary table
print("\n\n=== ARCHETYPE MAP ===")
print(f"{'Phase':<30} {'S_half':>6} {'avg_MI':>7} {'MI_CV':>7} {'avg_I3':>7} {'Best Archetype':<25}")
print("-" * 90)
for label, feat in results.items():
    print(f"{label:<30} {feat['half_cut_entropy']:>6.3f} {feat['avg_MI']:>7.4f} "
          f"{feat['mi_uniformity']:>7.3f} {feat['avg_I3']:>7.4f} {feat['best_archetype']:<25}")

print("\n=== COMPARISON WITH TFIM (Sprint 029) ===")
print("TFIM ordered (h→0):  GHZ-like, MI~1.0, I3~+1.0, CV~0.03")
print("TFIM critical (h=1): Scale-Free, MI~0.45, I3~0.0, CV~0.39")
print("TFIM disordered (h→∞): Product-like, MI→0, I3→0, CV~1.16")
print()
print("XXZ FM:    Should match TFIM ordered (both ferromagnetic)")
print("XXZ XY:    Critical throughout — should be Scale-Free?")
print("XXZ Néel:  Staggered order — novel archetype or variant?")

print("\nResults saved to results/sprint_030b_xxz_archetypes.json")
