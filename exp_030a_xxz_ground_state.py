"""
Sprint 030a: XXZ Ground State Entanglement Across Phase Diagram
Exact diagonalization of H = Σ (Sx_i Sx_{i+1} + Sy_i Sy_{i+1} + Δ Sz_i Sz_{i+1})
Sweep Δ from -2.0 to 3.0, track entropy, MI, I3, negativity
Three phases: FM (Δ<-1), XY/Luttinger (|Δ|<1), Néel (Δ>1)
Two transitions: first-order (Δ=-1), BKT (Δ=1)
"""

import numpy as np
from scipy.sparse import kron, eye, csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.linalg import sqrtm
import json, time

# Pauli matrices (sparse)
I2 = eye(2, format='csr')
Sx = csr_matrix(np.array([[0, 0.5], [0.5, 0]], dtype=complex))
Sy = csr_matrix(np.array([[0, -0.5j], [0.5j, 0]], dtype=complex))
Sz = csr_matrix(np.array([[0.5, 0], [0, -0.5]], dtype=complex))

# Raising/lowering for efficiency: Sx Sx + Sy Sy = 0.5*(S+ S- + S- S+)
Sp = csr_matrix(np.array([[0, 1], [0, 0]], dtype=complex))
Sm = csr_matrix(np.array([[0, 0], [1, 0]], dtype=complex))

def chain_op(op, qubit, n):
    """Single operator on qubit `qubit` in n-qubit system."""
    ops = [I2] * n
    ops[qubit] = op
    result = ops[0]
    for o in ops[1:]:
        result = kron(result, o, format='csr')
    return result

def xxz_hamiltonian(n, delta):
    """Build XXZ Hamiltonian: H = Σ (Sx_i Sx_{i+1} + Sy_i Sy_{i+1} + Δ Sz_i Sz_{i+1})
    Using: Sx Sx + Sy Sy = 0.5*(S+ S- + S- S+)
    Open boundary conditions.
    """
    dim = 2**n
    H = csr_matrix((dim, dim), dtype=complex)
    for i in range(n - 1):
        # XX + YY = 0.5 * (S+ S- + S- S+)
        H += 0.5 * (chain_op(Sp, i, n) @ chain_op(Sm, i+1, n) +
                     chain_op(Sm, i, n) @ chain_op(Sp, i+1, n))
        # Δ * ZZ
        H += delta * chain_op(Sz, i, n) @ chain_op(Sz, i+1, n)
    return H

def partial_trace_simple(state_vec, keep, n):
    """Partial trace from state vector, keeping qubits in `keep`."""
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
    """I3(i:j:k) = S_i + S_j + S_k - S_ij - S_ik - S_jk + S_ijk."""
    S_i = von_neumann_entropy(partial_trace_simple(state_vec, [i], n))
    S_j = von_neumann_entropy(partial_trace_simple(state_vec, [j], n))
    S_k = von_neumann_entropy(partial_trace_simple(state_vec, [k], n))
    S_ij = von_neumann_entropy(partial_trace_simple(state_vec, [i, j], n))
    S_ik = von_neumann_entropy(partial_trace_simple(state_vec, [i, k], n))
    S_jk = von_neumann_entropy(partial_trace_simple(state_vec, [j, k], n))
    S_ijk = von_neumann_entropy(partial_trace_simple(state_vec, [i, j, k], n))
    return S_i + S_j + S_k - S_ij - S_ik - S_jk + S_ijk

def negativity(state_vec, subsys_a, n):
    """Logarithmic negativity for bipartition."""
    rho = np.outer(state_vec, state_vec.conj())
    # Partial transpose w.r.t. subsys_a
    rho_reshaped = rho.reshape([2]*n + [2]*n)
    # Swap bra/ket indices for subsystem A
    perm = list(range(2*n))
    for q in subsys_a:
        perm[q], perm[q+n] = perm[q+n], perm[q]
    rho_pt = np.transpose(rho_reshaped, perm).reshape(2**n, 2**n)
    evals = np.linalg.eigvalsh(rho_pt)
    neg = (np.sum(np.abs(evals)) - 1.0) / 2.0
    return max(0.0, neg)

def half_cut_entropy(state_vec, n):
    """Entropy of half-cut bipartition."""
    keep = list(range(n // 2))
    rho_half = partial_trace_simple(state_vec, keep, n)
    return von_neumann_entropy(rho_half)

# ---- Main sweep ----
n = 6  # Safe for density matrix ops
print(f"=== XXZ Ground State Sweep, n={n} qubits ===")

# Dense sweep: 50 points from -2.0 to 3.0
# Extra points near transitions at Δ=-1 and Δ=1
delta_coarse = np.linspace(-2.0, 3.0, 40)
delta_fine_fm = np.linspace(-1.3, -0.7, 15)  # near FM transition
delta_fine_bkt = np.linspace(0.7, 1.3, 15)   # near BKT transition
delta_values = np.unique(np.sort(np.concatenate([delta_coarse, delta_fine_fm, delta_fine_bkt])))

results = {
    'n': n,
    'delta_values': [],
    'half_cut_entropy': [],
    'energy': [],
    'energy_gap': [],
    'single_qubit_entropy': [],
    'mi_matrix': [],  # all pairwise MI
    'avg_MI': [],
    'max_MI': [],
    'nearest_MI': [],
    'avg_I3': [],
    'min_I3': [],
    'max_I3': [],
    'consecutive_I3': [],
    'negativity_half': [],
}

t0 = time.time()

for delta in delta_values:
    H = xxz_hamiltonian(n, delta)

    # Ground state — need 2 eigenvalues for gap
    # For FM phase (Δ << -1), ground state is degenerate; use k=3 to be safe
    try:
        evals, evecs = eigsh(H, k=3, which='SA')
        gs_energy = evals[0]
        gap = evals[1] - evals[0]
        psi = evecs[:, 0]
    except Exception as e:
        print(f"  Δ={delta:.3f}: eigsh failed ({e}), skipping")
        continue

    # Half-cut entropy
    S_half = half_cut_entropy(psi, n)

    # Single-qubit entropy (avg)
    sq_entropies = []
    for i in range(n):
        sq_entropies.append(von_neumann_entropy(partial_trace_simple(psi, [i], n)))
    avg_sq = np.mean(sq_entropies)

    # Pairwise MI (all pairs)
    mi_mat = np.zeros((n, n))
    mi_values = []
    nn_mi = []
    for i in range(n):
        for j in range(i+1, n):
            mi = mutual_information(psi, i, j, n)
            mi_mat[i, j] = mi
            mi_mat[j, i] = mi
            mi_values.append(mi)
            if j == i + 1:
                nn_mi.append(mi)

    # I3 for all triples (feasible at n=6: C(6,3)=20)
    i3_values = []
    consec_i3 = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                i3 = tripartite_info(psi, i, j, k, n)
                i3_values.append(i3)
                if j == i+1 and k == j+1:
                    consec_i3.append(i3)

    # Half-cut negativity
    neg_half = negativity(psi, list(range(n//2)), n)

    results['delta_values'].append(float(delta))
    results['half_cut_entropy'].append(float(S_half))
    results['energy'].append(float(gs_energy))
    results['energy_gap'].append(float(gap))
    results['single_qubit_entropy'].append(float(avg_sq))
    results['mi_matrix'].append(mi_mat.tolist())
    results['avg_MI'].append(float(np.mean(mi_values)))
    results['max_MI'].append(float(np.max(mi_values)))
    results['nearest_MI'].append(float(np.mean(nn_mi)))
    results['avg_I3'].append(float(np.mean(i3_values)))
    results['min_I3'].append(float(np.min(i3_values)))
    results['max_I3'].append(float(np.max(i3_values)))
    results['consecutive_I3'].append(float(np.mean(consec_i3)) if consec_i3 else 0.0)
    results['negativity_half'].append(float(neg_half))

    print(f"  Δ={delta:+.3f}: E={gs_energy:.4f}, gap={gap:.4f}, S_half={S_half:.4f}, "
          f"avg_MI={np.mean(mi_values):.4f}, avg_I3={np.mean(i3_values):.4f}, neg={neg_half:.4f}")

elapsed = time.time() - t0
print(f"\nTotal time: {elapsed:.1f}s")

# Save results
with open('results/sprint_030a_xxz_ground_state.json', 'w') as f:
    json.dump(results, f, indent=2)

# Analysis
print("\n=== ANALYSIS ===")
d = np.array(results['delta_values'])
S = np.array(results['half_cut_entropy'])
MI = np.array(results['avg_MI'])
I3 = np.array(results['avg_I3'])
gap = np.array(results['energy_gap'])
neg = np.array(results['negativity_half'])

# Identify phases
fm_mask = d < -1.0
xy_mask = (d > -1.0) & (d < 1.0)
neel_mask = d > 1.0

print(f"\nFM phase (Δ<-1):")
if np.any(fm_mask):
    print(f"  Entropy: {np.mean(S[fm_mask]):.4f} ± {np.std(S[fm_mask]):.4f}")
    print(f"  Avg MI:  {np.mean(MI[fm_mask]):.4f} ± {np.std(MI[fm_mask]):.4f}")
    print(f"  Avg I3:  {np.mean(I3[fm_mask]):.4f} ± {np.std(I3[fm_mask]):.4f}")
    print(f"  Gap:     {np.mean(gap[fm_mask]):.4f}")
    print(f"  Negativity: {np.mean(neg[fm_mask]):.4f}")

print(f"\nXY phase (-1<Δ<1):")
if np.any(xy_mask):
    print(f"  Entropy: {np.mean(S[xy_mask]):.4f} ± {np.std(S[xy_mask]):.4f}")
    print(f"  Avg MI:  {np.mean(MI[xy_mask]):.4f} ± {np.std(MI[xy_mask]):.4f}")
    print(f"  Avg I3:  {np.mean(I3[xy_mask]):.4f} ± {np.std(I3[xy_mask]):.4f}")
    print(f"  Gap:     {np.mean(gap[xy_mask]):.4f}")
    print(f"  Negativity: {np.mean(neg[xy_mask]):.4f}")

print(f"\nNéel phase (Δ>1):")
if np.any(neel_mask):
    print(f"  Entropy: {np.mean(S[neel_mask]):.4f} ± {np.std(S[neel_mask]):.4f}")
    print(f"  Avg MI:  {np.mean(MI[neel_mask]):.4f} ± {np.std(MI[neel_mask]):.4f}")
    print(f"  Avg I3:  {np.mean(I3[neel_mask]):.4f} ± {np.std(I3[neel_mask]):.4f}")
    print(f"  Gap:     {np.mean(gap[neel_mask]):.4f}")
    print(f"  Negativity: {np.mean(neg[neel_mask]):.4f}")

# Find entropy peak
peak_idx = np.argmax(S)
print(f"\nEntropy peak at Δ = {d[peak_idx]:.3f}, S = {S[peak_idx]:.4f}")
print(f"Gap minimum at Δ = {d[np.argmin(gap)]:.3f}, gap = {min(gap):.6f}")

# MI uniformity CV (Sprint 029's new order parameter)
mi_cv_values = []
for mi_mat in results['mi_matrix']:
    mi_arr = np.array(mi_mat)
    upper = mi_arr[np.triu_indices(n, k=1)]
    if np.mean(upper) > 1e-10:
        cv = np.std(upper) / np.mean(upper)
    else:
        cv = float('inf')
    mi_cv_values.append(cv)

results['mi_uniformity_cv'] = mi_cv_values
print(f"\nMI uniformity CV:")
for phase, mask in [("FM", fm_mask), ("XY", xy_mask), ("Néel", neel_mask)]:
    cv_arr = np.array(mi_cv_values)[mask]
    if len(cv_arr) > 0:
        finite = cv_arr[np.isfinite(cv_arr)]
        if len(finite) > 0:
            print(f"  {phase}: {np.mean(finite):.4f} ± {np.std(finite):.4f}")
        else:
            print(f"  {phase}: inf (zero MI)")

# Re-save with CV
with open('results/sprint_030a_xxz_ground_state.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nResults saved to results/sprint_030a_xxz_ground_state.json")
