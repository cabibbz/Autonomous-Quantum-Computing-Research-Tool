#!/usr/bin/env python3
"""Sprint 031b: Entanglement spectrum across TFIM phase diagram.

Track how the entanglement spectrum evolves through the TFIM phase transition:
  ordered (h/J << 1) → critical (h/J = 1) → disordered (h/J >> 1)

Key questions:
1. Does the entanglement gap close at criticality?
2. Does the spectrum shape encode universality class information?
3. Can we see CFT tower structure at the critical point?
"""

import numpy as np
import json
from datetime import datetime

def tfim_ground_state(n, h_over_J):
    """TFIM ground state via exact diagonalization. H = -J sum(ZZ) - h sum(X)."""
    dim = 2**n
    H = np.zeros((dim, dim))
    for i in range(dim):
        for q in range(n):
            j = i ^ (1 << q)
            H[i, j] -= h_over_J  # transverse field
        for q in range(n - 1):
            zi = 1 - 2 * ((i >> q) & 1)
            zi1 = 1 - 2 * ((i >> (q + 1)) & 1)
            H[i, i] -= zi * zi1  # Ising coupling (J=1)
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    return eigenvectors[:, 0], eigenvalues[0], eigenvalues[1] - eigenvalues[0]

def half_cut_spectrum(sv, n):
    """Entanglement spectrum from half-cut."""
    nA = n // 2
    nB = n - nA
    psi = sv.reshape(2**nA, 2**nB)
    rho = psi @ psi.conj().T
    evals = np.linalg.eigvalsh(rho)
    evals = np.real(evals)
    evals = evals[evals > 1e-14]
    return np.sort(evals)[::-1]

def spectrum_metrics(evals):
    """Compute metrics from eigenvalue spectrum."""
    S_vn = -np.sum(evals * np.log2(evals + 1e-30))
    S2 = -np.log2(np.sum(evals**2))

    if len(evals) > 1:
        gap = evals[0] - evals[1]
        gap_ratio = evals[0] / evals[1]
    else:
        gap = evals[0]
        gap_ratio = np.inf

    PR = 1.0 / np.sum(evals**2)
    NPR = PR / len(evals) if len(evals) > 0 else 0
    schmidt_rank = len(evals)

    # Entanglement energies
    ent_energies = -np.log(evals + 1e-30)

    # Entanglement gap (in energy space): ξ₂ - ξ₁
    if len(ent_energies) > 1:
        ent_gap = float(ent_energies[1] - ent_energies[0])
    else:
        ent_gap = float('inf')

    # Level spacing statistics (for CFT tower detection)
    if len(ent_energies) > 2:
        spacings = np.diff(ent_energies)
        spacing_ratio = float(np.min(spacings[:3]) / np.max(spacings[:3])) if len(spacings) >= 3 else 0
    else:
        spacing_ratio = 0

    return {
        'von_neumann_entropy': float(S_vn),
        'renyi_2': float(S2),
        'spectral_gap': float(gap),
        'gap_ratio': float(min(gap_ratio, 1e10)),
        'entanglement_gap': ent_gap,
        'participation_ratio': float(PR),
        'normalized_PR': float(NPR),
        'schmidt_rank': int(schmidt_rank),
        'spacing_ratio': spacing_ratio,
        'eigenvalues': [float(e) for e in evals[:16]],
        'entanglement_energies': [float(e) for e in ent_energies[:10]],
    }


# ---- Main sweep ----
print("Sprint 031b: Entanglement spectrum across TFIM phase diagram")
print("=" * 60)

n = 8
# Dense sweep around critical point
h_values = np.concatenate([
    np.linspace(0.01, 0.5, 8),    # deep ordered
    np.linspace(0.6, 0.9, 6),     # approaching critical
    np.linspace(0.92, 1.08, 9),   # critical region (dense)
    np.linspace(1.1, 1.5, 6),     # leaving critical
    np.linspace(1.6, 3.0, 8),     # deep disordered
])

results = []
print(f"\nSweeping h/J from {h_values[0]:.2f} to {h_values[-1]:.2f} ({len(h_values)} points)")

for i, h in enumerate(h_values):
    sv, E0, energy_gap = tfim_ground_state(n, h)
    evals = half_cut_spectrum(sv, n)
    metrics = spectrum_metrics(evals)
    metrics['h_over_J'] = float(h)
    metrics['ground_energy'] = float(E0)
    metrics['energy_gap'] = float(energy_gap)
    results.append(metrics)

    if i % 10 == 0 or abs(h - 1.0) < 0.05:
        print(f"  h/J={h:.3f}: S={metrics['von_neumann_entropy']:.4f}, "
              f"rank={metrics['schmidt_rank']}, "
              f"ent_gap={metrics['entanglement_gap']:.3f}, "
              f"NPR={metrics['normalized_PR']:.4f}, "
              f"E_gap={energy_gap:.4f}")


# ---- Analysis ----
print("\n" + "=" * 60)
print("ANALYSIS")
print("=" * 60)

# Find critical point indicators
h_vals = np.array([r['h_over_J'] for r in results])
entropies = np.array([r['von_neumann_entropy'] for r in results])
ent_gaps = np.array([r['entanglement_gap'] for r in results])
nprs = np.array([r['normalized_PR'] for r in results])
ranks = np.array([r['schmidt_rank'] for r in results])
energy_gaps = np.array([r['energy_gap'] for r in results])

# Entropy peak
peak_idx = np.argmax(entropies)
print(f"\nEntropy peak: h/J={h_vals[peak_idx]:.3f}, S={entropies[peak_idx]:.4f}")

# Entanglement gap minimum
min_ent_gap_idx = np.argmin(ent_gaps)
print(f"Entanglement gap minimum: h/J={h_vals[min_ent_gap_idx]:.3f}, ξ_gap={ent_gaps[min_ent_gap_idx]:.4f}")

# Energy gap minimum (should be at h/J=1)
min_E_gap_idx = np.argmin(energy_gaps)
print(f"Energy gap minimum: h/J={h_vals[min_E_gap_idx]:.3f}, ΔE={energy_gaps[min_E_gap_idx]:.4f}")

# NPR peak
npr_peak_idx = np.argmax(nprs)
print(f"NPR peak: h/J={h_vals[npr_peak_idx]:.3f}, NPR={nprs[npr_peak_idx]:.4f}")

# Schmidt rank maximum
rank_max_idx = np.argmax(ranks)
print(f"Schmidt rank max: h/J={h_vals[rank_max_idx]:.3f}, rank={ranks[rank_max_idx]}")

# Three phases comparison
ordered = [r for r in results if r['h_over_J'] < 0.5]
critical = [r for r in results if 0.95 < r['h_over_J'] < 1.05]
disordered = [r for r in results if r['h_over_J'] > 2.0]

print(f"\nPhase comparison:")
print(f"  Ordered   (h<0.5):  S={np.mean([r['von_neumann_entropy'] for r in ordered]):.4f}, "
      f"rank={int(np.mean([r['schmidt_rank'] for r in ordered]))}, "
      f"ent_gap={np.mean([r['entanglement_gap'] for r in ordered]):.3f}")
print(f"  Critical  (≈1.0):   S={np.mean([r['von_neumann_entropy'] for r in critical]):.4f}, "
      f"rank={int(np.mean([r['schmidt_rank'] for r in critical]))}, "
      f"ent_gap={np.mean([r['entanglement_gap'] for r in critical]):.3f}")
print(f"  Disordered (h>2.0): S={np.mean([r['von_neumann_entropy'] for r in disordered]):.4f}, "
      f"rank={int(np.mean([r['schmidt_rank'] for r in disordered]))}, "
      f"ent_gap={np.mean([r['entanglement_gap'] for r in disordered]):.3f}")

# CFT tower analysis at critical point
print("\nCFT Tower Analysis at h/J=1.0:")
crit_result = min(results, key=lambda r: abs(r['h_over_J'] - 1.0))
ent_energies = crit_result['entanglement_energies']
print(f"  Entanglement energies: {[f'{e:.3f}' for e in ent_energies[:8]]}")
if len(ent_energies) > 1:
    # Shift to ground state = 0
    shifted = [e - ent_energies[0] for e in ent_energies]
    print(f"  Shifted (ξ - ξ₀):     {[f'{e:.3f}' for e in shifted[:8]]}")
    if shifted[1] > 0:
        normalized = [e / shifted[1] for e in shifted]
        print(f"  Normalized (Δξ/Δξ₁):  {[f'{e:.3f}' for e in normalized[:8]]}")
        print(f"  (CFT would give: 0, 1, integers or simple fractions)")

# Entanglement gap vs energy gap correlation
print("\nEntanglement gap vs Energy gap (selected points):")
for r in results:
    if r['h_over_J'] in [0.01, 0.5, 0.8, 0.96, 1.0, 1.04, 1.2, 1.5, 2.0, 3.0]:
        continue
    if abs(r['h_over_J'] - 0.5) < 0.02 or abs(r['h_over_J'] - 1.0) < 0.02 or abs(r['h_over_J'] - 2.0) < 0.02:
        print(f"  h/J={r['h_over_J']:.2f}: ξ_gap={r['entanglement_gap']:.3f}, ΔE={r['energy_gap']:.4f}")


# ---- Save ----
output = {
    'experiment': 'sprint_031b_tfim_entanglement_spectrum',
    'timestamp': datetime.now().isoformat(),
    'n_qubits': n,
    'description': 'Entanglement spectrum across TFIM phase diagram, n=8',
    'h_values': [float(h) for h in h_values],
    'results': results,
}

with open('results/sprint_031b_tfim_spectrum.json', 'w') as f:
    json.dump(output, f, indent=2)

print("\nResults saved to results/sprint_031b_tfim_spectrum.json")
