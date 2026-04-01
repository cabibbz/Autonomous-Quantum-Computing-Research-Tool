"""
Sprint 033b: 3-State Potts Entanglement Spectrum at Criticality

Compare entanglement spectrum structure across the Potts phase transition.
The critical Potts model has c=4/5 CFT with W₃ algebra, which should produce
different multiplet structure than TFIM's c=1/2 Ising CFT (Sprint 031b).

Key questions:
1. Does the Potts spectrum show triplet structure from Z₃ symmetry (vs doublet from Z₂)?
2. What is the entanglement gap behavior across the transition?
3. Does the Schmidt rank show a different pattern than TFIM?

System: n=8 (dim=6561), half-cut A = left 4 sites (dim_A = 81).
"""

import numpy as np
from scipy import linalg as la
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import eigsh
import json, time

t0 = time.time()

# ---- Sparse Potts Hamiltonian (reuse from 033a) ----
def sparse_potts_hamiltonian(n, h_over_J, boundary='open'):
    d = 3
    dim = d**n
    P = csr_matrix(np.array([[0,0,1],[1,0,0],[0,1,0]], dtype=float))
    Pd = P.T.tocsr()
    I_d = sp_eye(d, format='csr')
    delta = csr_matrix((d*d, d*d))
    for a in range(d):
        idx = a * d + a
        row = csr_matrix(([1.0], ([idx], [idx])), shape=(d*d, d*d))
        delta = delta + row
    H = csr_matrix((dim, dim))
    n_bonds = n - 1 if boundary == 'open' else n
    for bond in range(n_bonds):
        j = (bond + 1) % n
        if j == bond + 1:
            prefix = sp_eye(d**bond, format='csr') if bond > 0 else sp_eye(1, format='csr')
            suffix = sp_eye(d**(n - bond - 2), format='csr') if bond + 2 < n else sp_eye(1, format='csr')
            H = H - sp_kron(sp_kron(prefix, delta), suffix, format='csr')
        else:
            for a in range(d):
                proj = csr_matrix(([1.0], ([a], [a])), shape=(d, d))
                middle = sp_eye(d**(n-2), format='csr')
                H = H - sp_kron(sp_kron(proj, middle), proj, format='csr')
    PPd = P + Pd
    for site in range(n):
        prefix = sp_eye(d**site, format='csr') if site > 0 else sp_eye(1, format='csr')
        suffix = sp_eye(d**(n - site - 1), format='csr') if site < n - 1 else sp_eye(1, format='csr')
        H = H - h_over_J * sp_kron(sp_kron(prefix, PPd), suffix, format='csr')
    return H

def partial_trace_general(rho, dims, keep):
    n = len(dims)
    trace_out = [i for i in range(n) if i not in keep]
    shape = list(dims) + list(dims)
    rho_tensor = rho.reshape(shape)
    for q in sorted(trace_out, reverse=True):
        n_remaining = rho_tensor.ndim // 2
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q + n_remaining)
    dim_keep = int(np.prod([dims[k] for k in keep]))
    return rho_tensor.reshape(dim_keep, dim_keep)

# ---- Main sweep: entanglement spectrum across transition ----
n = 8
d = 3
dims = [d] * n
n_A = n // 2
dim_A = d**n_A  # = 81

# Sweep h/J through the transition (based on 033a: critical region ≈ 0.15-0.35)
h_values = np.concatenate([
    [0.01, 0.05],              # deep ordered
    np.linspace(0.10, 0.35, 11),  # transition region
    [0.40, 0.50, 0.75, 1.0, 2.0]  # disordered
])

print(f"=== Sprint 033b: Potts Entanglement Spectrum (n={n}, dim_A={dim_A}) ===")
print(f"Sweeping {len(h_values)} h/J values")

all_spectra = []
for h_J in h_values:
    t1 = time.time()
    H = sparse_potts_hamiltonian(n, h_J, 'open')
    eigvals_H, eigvecs_H = eigsh(H, k=4, which='SA')
    sort_idx = np.argsort(eigvals_H)
    psi_gs = eigvecs_H[:, sort_idx[0]]

    rho = np.outer(psi_gs, psi_gs.conj())
    rho_A = partial_trace_general(rho, dims, list(range(n_A)))

    # Entanglement spectrum: ξ_i = -ln(λ_i)
    eigvals_rho = la.eigvalsh(rho_A)
    eigvals_rho = np.sort(eigvals_rho)[::-1]  # descending

    # Significant eigenvalues
    sig = eigvals_rho[eigvals_rho > 1e-12]
    ent_energies = -np.log(sig)

    # Normalize: shift so ground state ξ₀ = 0
    ent_energies_shifted = ent_energies - ent_energies[0]

    # Entropy
    S = -np.sum(sig * np.log2(sig))

    # Schmidt rank
    rank = len(sig)

    # Entanglement gap
    if len(ent_energies_shifted) > 1:
        ent_gap = ent_energies_shifted[1]
    else:
        ent_gap = float('inf')

    # Degeneracy analysis: count near-degenerate groups
    tol_degen = 0.05  # tolerance for degeneracy
    groups = []
    current_group = [ent_energies_shifted[0]]
    for k in range(1, len(ent_energies_shifted)):
        if ent_energies_shifted[k] - ent_energies_shifted[k-1] < tol_degen:
            current_group.append(ent_energies_shifted[k])
        else:
            groups.append(len(current_group))
            current_group = [ent_energies_shifted[k]]
    groups.append(len(current_group))

    # Normalized participation ratio
    npr = np.sum(sig**2)  # Tr(ρ²), lower = more spread
    npr_inv = 1.0 / (rank * npr) if npr > 0 else 0  # normalized, 1 = flat

    result = {
        'h_J': float(h_J),
        'entropy_bits': float(S),
        'schmidt_rank': int(rank),
        'ent_gap': float(ent_gap),
        'ent_energies': [float(x) for x in ent_energies_shifted[:20]],
        'top_eigenvalues': [float(x) for x in sig[:20]],
        'degeneracy_groups': groups[:15],
        'npr': float(npr_inv),
        'time_s': float(time.time() - t1)
    }
    all_spectra.append(result)

    # Print summary
    group_str = str(groups[:8])
    print(f"  h/J={h_J:.3f}: S={S:.3f}, rank={rank:2d}, gap={ent_gap:.3f}, "
          f"groups={group_str} [{time.time()-t1:.1f}s]")

# ---- Detailed spectrum analysis at critical point ----
print("\n=== Spectrum at Critical Region (h/J ≈ 0.225) ===")
# Use h/J=0.225 (near steepest entropy drop / critical point)
crit_idx = min(range(len(all_spectra)), key=lambda i: abs(all_spectra[i]['h_J'] - 0.225))
crit = all_spectra[crit_idx]
print(f"h/J = {crit['h_J']:.3f}")
print(f"Entropy = {crit['entropy_bits']:.4f} bits")
print(f"Schmidt rank = {crit['schmidt_rank']}")
print(f"Entanglement gap = {crit['ent_gap']:.4f}")
print(f"\nEntanglement energies (ξ - ξ₀):")
for k, xi in enumerate(crit['ent_energies'][:15]):
    print(f"  ξ_{k} = {xi:.4f}")
print(f"\nDegeneracy pattern: {crit['degeneracy_groups'][:12]}")

# Check for triplet structure (Z₃ → expect degeneracies in multiples of 3?)
print("\n=== Triplet Analysis ===")
for sp in all_spectra:
    groups = sp['degeneracy_groups']
    n_triplets = sum(1 for g in groups if g == 3)
    n_doublets = sum(1 for g in groups if g == 2)
    n_singlets = sum(1 for g in groups if g == 1)
    if sp['h_J'] in [0.01, 0.15, 0.225, 0.30, 0.50, 1.0] or abs(sp['h_J'] - 0.225) < 0.03:
        print(f"  h/J={sp['h_J']:.3f}: singlets={n_singlets}, doublets={n_doublets}, "
              f"triplets={n_triplets}, groups={groups[:10]}")

# ---- Spacing analysis at criticality ----
print("\n=== Entanglement Energy Spacings at Criticality ===")
xi = crit['ent_energies']
spacings = [xi[k+1] - xi[k] for k in range(min(14, len(xi)-1))]
print("Spacings (ξ_{k+1} - ξ_k):")
for k, s in enumerate(spacings):
    print(f"  Δξ_{k}→{k+1} = {s:.4f}")

# Ratio test: for CFT, spacings should be related to scaling dimensions
print("\nSpacing ratios (Δξ_k / Δξ_0):")
if spacings[0] > 1e-10:
    for k, s in enumerate(spacings[:10]):
        print(f"  Δξ_{k}/{k+1} / Δξ_0 = {s/spacings[0]:.4f}")

# ---- Compare to TFIM pattern (Sprint 031b) ----
print("\n=== Comparison to TFIM (Sprint 031b) ===")
print("TFIM c=1/2: doublet pattern (groups [1,1,...] with alternating spacings ~1.0 and ~2.75)")
print("Potts c=4/5: expect different CFT content — larger central charge, W₃ algebra")
print(f"Potts critical spectrum: groups = {crit['degeneracy_groups'][:12]}")

# ---- Evolution of spectrum across transition ----
print("\n=== Spectrum Evolution ===")
print("  h/J    rank  gap    NPR    dominant_group_sizes")
for sp in all_spectra:
    print(f"  {sp['h_J']:5.3f}  {sp['schmidt_rank']:3d}   {sp['ent_gap']:5.3f}  {sp['npr']:.4f}  {sp['degeneracy_groups'][:8]}")

# ---- Save ----
results = {
    'experiment': '033b',
    'description': 'Potts entanglement spectrum across phase transition',
    'parameters': {'n': n, 'local_dim': d, 'n_A': n_A, 'dim_A': dim_A},
    'spectra': all_spectra,
    'critical_analysis': {
        'h_J': crit['h_J'],
        'spacings': spacings,
        'degeneracy_groups': crit['degeneracy_groups'],
    },
    'runtime_seconds': time.time() - t0
}

with open('results/sprint_033b_potts_spectrum.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nSaved. Runtime: {time.time()-t0:.1f}s")
