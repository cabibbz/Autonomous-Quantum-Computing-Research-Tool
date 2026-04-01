#!/usr/bin/env python3
"""Sprint 031c: Entanglement spectrum across XXZ phase diagram.

Track entanglement spectrum across XXZ model: FM → XY → Néel.
H = -sum(XX + YY + Δ ZZ)

Key questions:
1. Does I3 sign change at Δ≈0.7 correspond to a spectrum transition?
2. Do FM (Δ=-1) and BKT (Δ=1) transitions show different spectral signatures?
3. Is the entanglement spectrum "dissociated" from bulk transitions?
"""

import numpy as np
import json
from datetime import datetime

def xxz_ground_state(n, delta):
    """XXZ ground state via exact diag. H = -sum(XX + YY + Δ*ZZ)."""
    dim = 2**n
    H = np.zeros((dim, dim))

    for i in range(dim):
        for q in range(n - 1):
            # ZZ coupling
            zi = 1 - 2 * ((i >> q) & 1)
            zi1 = 1 - 2 * ((i >> (q + 1)) & 1)
            H[i, i] -= delta * zi * zi1

            # XX + YY = 2(|01><10| + |10><01|) = flip-flop
            bq = (i >> q) & 1
            bq1 = (i >> (q + 1)) & 1
            if bq != bq1:  # one is 0, other is 1
                j = i ^ (1 << q) ^ (1 << (q + 1))  # flip both
                H[i, j] -= 1.0  # XX contribution
                H[i, j] -= 1.0  # YY contribution (real part for flip-flop)

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

def pairwise_mi_and_i3(sv, n):
    """Compute average pairwise MI and sample I3 for archetype comparison."""
    dim = 2**n
    rho_full = np.outer(sv, sv.conj())

    # Single-qubit entropies
    s1 = []
    for q in range(n):
        # Trace out all but qubit q
        rho_q = np.zeros((2, 2), dtype=complex)
        for a in range(2):
            for b in range(2):
                for other in range(dim):
                    i_row = (other & ~(1 << q)) | (a << q)
                    i_col = (other & ~(1 << q)) | (b << q)
                    rho_q[a, b] += rho_full[i_row, i_col]
        evals_q = np.linalg.eigvalsh(rho_q)
        evals_q = evals_q[evals_q > 1e-14]
        s1.append(-np.sum(evals_q * np.log2(evals_q + 1e-30)))

    # Nearest-neighbor 2-qubit RDMs for MI
    mi_values = []
    for q in range(n - 1):
        rho_2 = np.zeros((4, 4), dtype=complex)
        for a in range(4):
            for b in range(4):
                a1, a2 = a >> 1, a & 1
                b1, b2 = b >> 1, b & 1
                for other_bits in range(2**(n-2)):
                    # Build full index with qubits q and q+1 set
                    idx_a = 0
                    idx_b = 0
                    bit_pos = 0
                    for p in range(n):
                        if p == q:
                            idx_a |= (a1 << p)
                            idx_b |= (b1 << p)
                        elif p == q + 1:
                            idx_a |= (a2 << p)
                            idx_b |= (b2 << p)
                        else:
                            bit_val = (other_bits >> bit_pos) & 1
                            idx_a |= (bit_val << p)
                            idx_b |= (bit_val << p)
                            bit_pos += 1
                    rho_2[a, b] += rho_full[idx_a, idx_b]

        evals_2 = np.linalg.eigvalsh(rho_2)
        evals_2 = evals_2[evals_2 > 1e-14]
        s2 = -np.sum(evals_2 * np.log2(evals_2 + 1e-30))
        mi = s1[q] + s1[q+1] - s2
        mi_values.append(float(mi))

    avg_mi = np.mean(mi_values)
    mi_cv = np.std(mi_values) / avg_mi if avg_mi > 1e-6 else 0

    return float(avg_mi), float(mi_cv), mi_values

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
    ent_energies = -np.log(evals + 1e-30)

    ent_gap = float(ent_energies[1] - ent_energies[0]) if len(ent_energies) > 1 else float('inf')

    # Doublet detection: check if consecutive entanglement energies come in near-degenerate pairs
    if len(ent_energies) >= 4:
        shifted = ent_energies - ent_energies[0]
        spacings = np.diff(shifted)
        if len(spacings) >= 3:
            # Alternating pattern: small, large, small, large...
            even_spacings = spacings[0::2][:4]  # 0, 2, 4, ...
            odd_spacings = spacings[1::2][:4]   # 1, 3, 5, ...
            alt_ratio = np.mean(even_spacings) / np.mean(odd_spacings) if np.mean(odd_spacings) > 1e-6 else 0
        else:
            alt_ratio = 0
    else:
        alt_ratio = 0

    # Degeneracy count: number of eigenvalue pairs within 1% of each other
    degeneracies = 0
    for i in range(len(evals) - 1):
        if abs(evals[i] - evals[i+1]) / (evals[i] + 1e-30) < 0.01:
            degeneracies += 1

    return {
        'von_neumann_entropy': float(S_vn),
        'renyi_2': float(S2),
        'spectral_gap': float(gap),
        'gap_ratio': float(min(gap_ratio, 1e10)),
        'entanglement_gap': ent_gap,
        'participation_ratio': float(PR),
        'normalized_PR': float(NPR),
        'schmidt_rank': int(len(evals)),
        'alternating_ratio': float(alt_ratio),
        'degeneracies': int(degeneracies),
        'eigenvalues': [float(e) for e in evals[:16]],
        'entanglement_energies': [float(e) for e in ent_energies[:10]],
    }


# ---- Main sweep ----
print("Sprint 031c: Entanglement spectrum across XXZ phase diagram")
print("=" * 60)

n = 6  # Use n=6 for speed (MI/I3 computation is expensive)

# Dense sweep
delta_values = np.concatenate([
    np.linspace(-2.0, -1.2, 5),   # deep FM
    np.linspace(-1.1, -0.8, 7),   # FM transition region
    np.linspace(-0.7, 0.5, 7),    # XY phase
    np.linspace(0.6, 0.8, 5),     # I3 sign change region (Δ≈0.7)
    np.linspace(0.85, 1.15, 7),   # BKT transition region
    np.linspace(1.2, 3.0, 6),     # Néel phase
])

results = []
print(f"\nSweeping Δ from {delta_values[0]:.2f} to {delta_values[-1]:.2f} ({len(delta_values)} points)")

for i, delta in enumerate(delta_values):
    sv, E0, energy_gap = xxz_ground_state(n, delta)
    evals = half_cut_spectrum(sv, n)
    metrics = spectrum_metrics(evals)
    metrics['delta'] = float(delta)
    metrics['ground_energy'] = float(E0)
    metrics['energy_gap'] = float(energy_gap)

    # Compute MI and MI-CV
    avg_mi, mi_cv, mi_list = pairwise_mi_and_i3(sv, n)
    metrics['avg_mi'] = avg_mi
    metrics['mi_cv'] = mi_cv
    metrics['nn_mi'] = mi_list

    results.append(metrics)

    if i % 8 == 0 or abs(delta - 0.7) < 0.1 or abs(delta + 1.0) < 0.1 or abs(delta - 1.0) < 0.1:
        print(f"  Δ={delta:+.3f}: S={metrics['von_neumann_entropy']:.4f}, "
              f"rank={metrics['schmidt_rank']}, "
              f"ent_gap={metrics['entanglement_gap']:.3f}, "
              f"degen={metrics['degeneracies']}, "
              f"MI_CV={mi_cv:.3f}")


# ---- Analysis ----
print("\n" + "=" * 60)
print("ANALYSIS")
print("=" * 60)

deltas = np.array([r['delta'] for r in results])
entropies = np.array([r['von_neumann_entropy'] for r in results])
ent_gaps = np.array([r['entanglement_gap'] for r in results])
nprs = np.array([r['normalized_PR'] for r in results])
ranks = np.array([r['schmidt_rank'] for r in results])
degens = np.array([r['degeneracies'] for r in results])
mi_cvs = np.array([r['mi_cv'] for r in results])
energy_gaps = np.array([r['energy_gap'] for r in results])

# Key transitions
print("\nKey transition indicators:")
# FM transition (Δ=-1)
fm_region = [(r['delta'], r) for r in results if -1.2 < r['delta'] < -0.8]
print(f"\nFM transition (Δ=-1):")
for d, r in fm_region:
    print(f"  Δ={d:+.2f}: S={r['von_neumann_entropy']:.4f}, ent_gap={r['entanglement_gap']:.3f}, "
          f"degen={r['degeneracies']}, rank={r['schmidt_rank']}, E_gap={r['energy_gap']:.4f}")

# I3 sign change region (Δ≈0.7)
i3_region = [(r['delta'], r) for r in results if 0.5 < r['delta'] < 0.9]
print(f"\nI3 sign change region (Δ≈0.7):")
for d, r in i3_region:
    print(f"  Δ={d:+.2f}: S={r['von_neumann_entropy']:.4f}, ent_gap={r['entanglement_gap']:.3f}, "
          f"degen={r['degeneracies']}, rank={r['schmidt_rank']}, MI_CV={r['mi_cv']:.3f}")

# BKT transition (Δ=1)
bkt_region = [(r['delta'], r) for r in results if 0.85 < r['delta'] < 1.2]
print(f"\nBKT transition (Δ=1):")
for d, r in bkt_region:
    print(f"  Δ={d:+.2f}: S={r['von_neumann_entropy']:.4f}, ent_gap={r['entanglement_gap']:.3f}, "
          f"degen={r['degeneracies']}, rank={r['schmidt_rank']}, E_gap={r['energy_gap']:.4f}")

# Spectrum shape evolution
print("\nSpectrum shape evolution:")
for r in results[::4]:
    evals = r['eigenvalues'][:6]
    evals_str = ', '.join(f'{e:.4f}' for e in evals)
    print(f"  Δ={r['delta']:+.2f}: [{evals_str}]  rank={r['schmidt_rank']}")

# Degeneracy transitions
print(f"\nDegeneracy count across phase diagram:")
for r in results:
    if r['degeneracies'] > 0 or abs(r['delta'] - round(r['delta'])) < 0.1:
        print(f"  Δ={r['delta']:+.2f}: {r['degeneracies']} near-degenerate pairs")

# Entanglement gap extrema
min_gap_idx = np.argmin(ent_gaps)
max_gap_idx = np.argmax(ent_gaps)
print(f"\nEntanglement gap: min at Δ={deltas[min_gap_idx]:.2f} (ξ_gap={ent_gaps[min_gap_idx]:.4f}), "
      f"max at Δ={deltas[max_gap_idx]:.2f} (ξ_gap={ent_gaps[max_gap_idx]:.4f})")

# Does entanglement gap track energy gap?
print("\nCorrelation check: entanglement gap vs energy gap")
valid = energy_gaps > 1e-6
corr = np.corrcoef(ent_gaps[valid], energy_gaps[valid])[0, 1]
print(f"  Pearson correlation (where E_gap > 0): {corr:.4f}")


# ---- Save ----
output = {
    'experiment': 'sprint_031c_xxz_entanglement_spectrum',
    'timestamp': datetime.now().isoformat(),
    'n_qubits': n,
    'description': 'Entanglement spectrum across XXZ phase diagram, n=6',
    'delta_values': [float(d) for d in delta_values],
    'results': results,
}

with open('results/sprint_031c_xxz_spectrum.json', 'w') as f:
    json.dump(output, f, indent=2)

print("\nResults saved to results/sprint_031c_xxz_spectrum.json")
