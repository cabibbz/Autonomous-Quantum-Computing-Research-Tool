#!/usr/bin/env python3
"""Sprint 031a: Entanglement spectrum of 5 entanglement archetypes.

Compute the full eigenvalue spectrum of the half-cut reduced density matrix for:
1. GHZ (Democratic) - n=8
2. W (Distributed) - n=8
3. Cluster 1D (Geometric) - n=8
4. Cluster 2D (Topological) - 2x4=8 qubits
5. Scale-Free (TFIM critical) - n=8

Characterize each spectrum: shape, gap, participation ratio, spectral entropy vs von Neumann.
"""

import numpy as np
import json
from datetime import datetime
from qiskit.circuit import QuantumCircuit
from qiskit_aer import AerSimulator

def get_statevector(qc):
    """Get statevector from circuit."""
    sim = AerSimulator(method='statevector')
    qc_copy = qc.copy()
    qc_copy.save_statevector()
    result = sim.run(qc_copy).result()
    return np.array(result.get_statevector(qc_copy))

def half_cut_rdm(sv, n):
    """Get reduced density matrix for first n//2 qubits."""
    nA = n // 2
    nB = n - nA
    psi = sv.reshape(2**nA, 2**nB)
    rho = psi @ psi.conj().T
    return rho

def entanglement_spectrum(rho):
    """Compute entanglement spectrum (eigenvalues of RDM, sorted descending)."""
    eigenvalues = np.linalg.eigvalsh(rho)
    # Filter out numerical noise
    eigenvalues = np.real(eigenvalues)
    eigenvalues = eigenvalues[eigenvalues > 1e-12]
    eigenvalues = np.sort(eigenvalues)[::-1]
    return eigenvalues

def spectrum_metrics(evals):
    """Compute metrics characterizing the spectrum shape."""
    # Von Neumann entropy
    S_vn = -np.sum(evals * np.log2(evals + 1e-30))

    # Spectral gap: ratio of largest to second-largest eigenvalue
    if len(evals) > 1:
        gap = evals[0] - evals[1]  # absolute gap
        gap_ratio = evals[0] / evals[1] if evals[1] > 1e-12 else np.inf
    else:
        gap = 0
        gap_ratio = np.inf

    # Participation ratio: 1/sum(λ²), measures effective number of terms
    PR = 1.0 / np.sum(evals**2)

    # Normalized participation ratio: PR / dim
    NPR = PR / len(evals)

    # Renyi-2 entropy: -log2(sum(λ²))
    S2 = -np.log2(np.sum(evals**2))

    # Entanglement energy spectrum: ξ_i = -ln(λ_i)
    ent_energies = -np.log(evals + 1e-30)

    # Spectrum flatness: std of eigenvalues / mean
    flatness = np.std(evals) / np.mean(evals) if np.mean(evals) > 0 else 0

    # Schmidt rank (number of nonzero eigenvalues)
    schmidt_rank = len(evals)

    # Max entropy for this rank
    S_max = np.log2(schmidt_rank)

    # Entropy ratio: how close to maximally mixed
    entropy_ratio = S_vn / S_max if S_max > 0 else 0

    return {
        'von_neumann_entropy': float(S_vn),
        'renyi_2_entropy': float(S2),
        'spectral_gap': float(gap),
        'gap_ratio': float(min(gap_ratio, 1e6)),
        'participation_ratio': float(PR),
        'normalized_PR': float(NPR),
        'flatness': float(flatness),
        'schmidt_rank': int(schmidt_rank),
        'max_entropy': float(S_max),
        'entropy_ratio': float(entropy_ratio),
        'entanglement_energies': [float(e) for e in ent_energies[:10]],  # top 10
        'eigenvalues': [float(e) for e in evals[:16]],  # top 16
    }


# ---- State preparation circuits ----

def make_ghz(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(1, n):
        qc.cx(0, i)
    return qc

def make_w(n):
    """W state: equal superposition of single-excitation states."""
    qc = QuantumCircuit(n)
    qc.x(0)
    for i in range(n - 1):
        theta = 2 * np.arccos(np.sqrt(1 / (n - i)))
        qc.ry(theta, i)
        qc.cx(i, i + 1)
        qc.x(i)
    return qc

def make_cluster_1d(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        qc.cz(i, i + 1)
    return qc

def make_cluster_2d(rows, cols):
    """2D cluster state on rows x cols grid."""
    n = rows * cols
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    # Horizontal edges
    for r in range(rows):
        for c in range(cols - 1):
            qc.cz(r * cols + c, r * cols + c + 1)
    # Vertical edges
    for r in range(rows - 1):
        for c in range(cols):
            qc.cz(r * cols + c, (r + 1) * cols + c)
    return qc

def make_tfim_critical(n):
    """TFIM ground state at critical point (h/J = 1) via exact diagonalization."""
    dim = 2**n
    H = np.zeros((dim, dim))

    for i in range(dim):
        # Transverse field: -h * sum(X_i)
        for q in range(n):
            j = i ^ (1 << q)  # flip qubit q
            H[i, j] -= 1.0  # h = J = 1

        # Ising coupling: -J * sum(Z_i Z_{i+1})
        for q in range(n - 1):
            zi = 1 - 2 * ((i >> q) & 1)
            zi1 = 1 - 2 * ((i >> (q + 1)) & 1)
            H[i, i] -= zi * zi1  # J = 1

    eigenvalues, eigenvectors = np.linalg.eigh(H)
    return eigenvectors[:, 0]  # ground state


# ---- Main ----

print("Sprint 031a: Entanglement spectrum of 5 archetypes")
print("=" * 60)

n = 8
results = {}

# 1. GHZ
print("\n1. GHZ (Democratic)...")
sv_ghz = get_statevector(make_ghz(n))
rho_ghz = half_cut_rdm(sv_ghz, n)
evals_ghz = entanglement_spectrum(rho_ghz)
metrics_ghz = spectrum_metrics(evals_ghz)
results['GHZ'] = metrics_ghz
print(f"   Eigenvalues: {evals_ghz}")
print(f"   S_vN={metrics_ghz['von_neumann_entropy']:.4f}, PR={metrics_ghz['participation_ratio']:.2f}, gap={metrics_ghz['spectral_gap']:.4f}")

# 2. W
print("\n2. W (Distributed)...")
sv_w = get_statevector(make_w(n))
rho_w = half_cut_rdm(sv_w, n)
evals_w = entanglement_spectrum(rho_w)
metrics_w = spectrum_metrics(evals_w)
results['W'] = metrics_w
print(f"   Eigenvalues: {evals_w[:8]}")
print(f"   S_vN={metrics_w['von_neumann_entropy']:.4f}, PR={metrics_w['participation_ratio']:.2f}, gap={metrics_w['spectral_gap']:.4f}")

# 3. Cluster 1D
print("\n3. Cluster 1D (Geometric)...")
sv_cl1d = get_statevector(make_cluster_1d(n))
rho_cl1d = half_cut_rdm(sv_cl1d, n)
evals_cl1d = entanglement_spectrum(rho_cl1d)
metrics_cl1d = spectrum_metrics(evals_cl1d)
results['Cluster_1D'] = metrics_cl1d
print(f"   Eigenvalues: {evals_cl1d[:8]}")
print(f"   S_vN={metrics_cl1d['von_neumann_entropy']:.4f}, PR={metrics_cl1d['participation_ratio']:.2f}, gap={metrics_cl1d['spectral_gap']:.4f}")

# 4. Cluster 2D (2x4)
print("\n4. Cluster 2D (Topological) [2x4]...")
sv_cl2d = get_statevector(make_cluster_2d(2, 4))
rho_cl2d = half_cut_rdm(sv_cl2d, n)
evals_cl2d = entanglement_spectrum(rho_cl2d)
metrics_cl2d = spectrum_metrics(evals_cl2d)
results['Cluster_2D'] = metrics_cl2d
print(f"   Eigenvalues: {evals_cl2d[:8]}")
print(f"   S_vN={metrics_cl2d['von_neumann_entropy']:.4f}, PR={metrics_cl2d['participation_ratio']:.2f}, gap={metrics_cl2d['spectral_gap']:.4f}")

# 5. TFIM Critical (Scale-Free)
print("\n5. TFIM Critical (Scale-Free) [h/J=1.0]...")
sv_tfim = make_tfim_critical(n)
rho_tfim = half_cut_rdm(sv_tfim, n)
evals_tfim = entanglement_spectrum(rho_tfim)
metrics_tfim = spectrum_metrics(evals_tfim)
results['Scale_Free'] = metrics_tfim
print(f"   Eigenvalues: {evals_tfim[:8]}")
print(f"   S_vN={metrics_tfim['von_neumann_entropy']:.4f}, PR={metrics_tfim['participation_ratio']:.2f}, gap={metrics_tfim['spectral_gap']:.4f}")

# ---- Comparative analysis ----
print("\n" + "=" * 60)
print("COMPARATIVE ANALYSIS")
print("=" * 60)

print(f"\n{'Archetype':<15} {'S_vN':>6} {'S2':>6} {'PR':>6} {'NPR':>6} {'Gap':>8} {'Rank':>5} {'S/Smax':>6} {'Flat':>6}")
print("-" * 75)
for name in ['GHZ', 'W', 'Cluster_1D', 'Cluster_2D', 'Scale_Free']:
    m = results[name]
    print(f"{name:<15} {m['von_neumann_entropy']:6.3f} {m['renyi_2_entropy']:6.3f} "
          f"{m['participation_ratio']:6.2f} {m['normalized_PR']:6.3f} "
          f"{m['spectral_gap']:8.4f} {m['schmidt_rank']:5d} "
          f"{m['entropy_ratio']:6.3f} {m['flatness']:6.3f}")

print("\nEntanglement energies (ξ = -ln λ):")
for name in ['GHZ', 'W', 'Cluster_1D', 'Cluster_2D', 'Scale_Free']:
    m = results[name]
    energies_str = ', '.join(f'{e:.3f}' for e in m['entanglement_energies'][:6])
    print(f"  {name:<15} [{energies_str}]")

# Spectrum shape classification
print("\nSpectrum Shape Classification:")
for name in ['GHZ', 'W', 'Cluster_1D', 'Cluster_2D', 'Scale_Free']:
    m = results[name]
    evals = np.array(m['eigenvalues'])
    if m['schmidt_rank'] <= 2:
        shape = "BINARY (2-level)"
    elif m['normalized_PR'] > 0.9:
        shape = "FLAT (maximally mixed)"
    elif m['normalized_PR'] > 0.5:
        shape = "BROAD"
    elif m['gap_ratio'] > 10:
        shape = "GAPPED (dominated by largest)"
    else:
        # Check for power-law
        if len(evals) >= 4:
            log_evals = np.log(evals[:min(8, len(evals))])
            log_idx = np.log(np.arange(1, len(log_evals) + 1))
            if len(log_idx) > 2:
                slope, _ = np.polyfit(log_idx, log_evals, 1)
                shape = f"POWER-LAW (α≈{-slope:.2f})" if abs(slope) > 0.3 else "GRADUAL"
            else:
                shape = "SHORT"
        else:
            shape = "SHORT"
    print(f"  {name:<15} → {shape}  (rank={m['schmidt_rank']}, NPR={m['normalized_PR']:.3f})")


# ---- Save results ----
output = {
    'experiment': 'sprint_031a_archetype_entanglement_spectra',
    'timestamp': datetime.now().isoformat(),
    'n_qubits': n,
    'description': 'Entanglement spectrum of 5 archetypes at n=8, half-cut bipartition',
    'results': results,
}

with open('results/sprint_031a_archetype_spectra.json', 'w') as f:
    json.dump(output, f, indent=2)

print("\nResults saved to results/sprint_031a_archetype_spectra.json")
