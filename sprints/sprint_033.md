# Sprint 033 — 3-State Potts Model: S₃ Symmetry Tests BW Locality Prediction

**Date:** 2026-04-01
**Status:** In Progress

## Idea

Sprint 032 discovered that symmetry group controls BW locality: Z₂ (TFIM) → 91%, U(1) (XXZ) → 100%. The 3-state Potts model has S₃ symmetry (6 elements — between Z₂=2 and U(1)=∞). This gives a clean quantitative test of the symmetry → BW accuracy prediction.

Additionally, the critical 3-state Potts model has central charge c=4/5 (vs c=1/2 for Ising), a W-algebra symmetry, and a non-diagonal partition function — all of which should produce distinct entanglement spectrum signatures compared to TFIM (Sprint 031).

## Literature Check

- Giudici et al. (2018), PRB 98, 134403: Tested BW for multiple 1D models including Potts. Confirmed BW works qualitatively. Did NOT compare BW accuracy quantitatively across symmetry groups.
- 3-state Potts: c=4/5, second-order transition, S₃ = Z₃ ⋊ Z₂ symmetry.
- Gap: No systematic comparison of BW accuracy vs symmetry group size at small scale. Our Sprint 032 data + Potts fills this gap.

## Potts Model

H = -J Σ δ(σᵢ, σⱼ) - h Σ (Pᵢ + Pᵢ†)

- σᵢ ∈ {0, 1, 2}, δ = Kronecker delta (Potts interaction)
- P: cyclic permutation |0⟩→|1⟩→|2⟩→|0⟩
- S₃ symmetry: Z₃ (cyclic) + Z₂ (P ↔ P†, charge conjugation)
- h/J = 0: ordered (3 degenerate ground states), h/J → ∞: disordered (unique)
- Critical point at h/J ≈ 0.2-0.3 for this convention (P+P† has max eigenvalue 2, shifting effective field scale)
- h/J → ∞: disordered (unique paramagnetic ground state)

## Experiments

### 33a: Potts Ground State Entanglement Sweep
**Status:** Complete

Swept h/J from 0.01 to 3.0 for n=6 (dim=729) with n=8 (dim=6561) at key points.

**Results (n=6):**
| Phase | h/J | S (bits) | MI_CV | I3(0,1,2) |
|-------|-----|----------|-------|-----------|
| Ordered | 0.01 | 1.585 = log₂(3) | 0.000 | +1.584 |
| Transition | 0.20 | 1.504 | 0.131 | +0.956 |
| Steepest drop | 0.27 | 1.172 | 0.348 | +0.430 |
| I3 sign change | 0.42 | ~0 | ~0.97 | ~0 |
| Disordered | 1.00 | 0.068 | 1.292 | -0.006 |

**Key findings:**
1. Ordered Potts = GHZ-3: S = log₂(3), uniform MI (CV=0), I3 = +1.58 (maximally redundant). Same archetype as TFIM but with q=3 degeneracy.
2. Steepest entropy drop at h/J ≈ 0.27 (transition region).
3. Entropy INCREASES with system size at h/J=0.25 (n=6: 1.315 → n=8: 1.426), confirming critical-like scaling.
4. I3 changes sign at h/J ≈ 0.42 — same pattern as TFIM (redundant → synergistic).
5. Disordered phase has near-zero entropy (product state), MI_CV → 1.35.
6. Pattern matches TFIM exactly: GHZ → Scale-Free-like → Product.

**Surprise:** The critical point is NOT at h/J=1 (self-dual) for this convention. The transverse field P+P† has max eigenvalue 2, making the effective field 2h and shifting the critical point to h/J ≈ 0.2-0.3.

### 33b: Potts Entanglement Spectrum at Criticality
**Status:** Complete

n=8 (dim_A = 81), swept h/J from 0.01 to 2.0.

**Results:**
| h/J | Phase | Rank | Gap | Degeneracy Pattern |
|-----|-------|------|-----|-------------------|
| 0.05 | Ordered | 9 | 0.001 | **[3, 3, 3]** — perfect triplets |
| 0.15 | Near-critical | 15 | 0.089 | [1, 2, 1, 2, 2, ...] |
| 0.225 | Critical | 24 | 0.609 | [1, 2, 1, 2, 2, 1, ...] |
| 0.50 | Disordered | 15 | 3.954 | [1, 2, 2, 1, 2, ...] |
| 1.00 | Deep disordered | 9 | 5.631 | [1, 2, 2, 1, 2, 1] |

**Key findings:**
1. **Z₃ → triplet structure in ordered phase**: GHZ-3 state has perfect [3,3,3] triplet degeneracy. This is the Z₃ fingerprint — NOT present in TFIM (which has Z₂ → doublets).
2. **Triplet → doublet transition**: As h/J increases, spectrum transitions from triplets [3,3,3] to doublets [1,2,1,2,...]. The Z₃ symmetry becomes invisible at the spectral level in disordered phase.
3. **Critical spectrum dominated by 3 levels**: Gap between ξ₂=0.609 and ξ₃=7.43 is huge. The dominant levels are singlet + doublet.
4. **Schmidt rank peaks at 24** near criticality (vs max 81).
5. **NPR minimum at criticality** (0.10) — most spread-out spectrum.

**Comparison to TFIM (Sprint 031b):**
- TFIM: doublets throughout (Z₂). Potts: triplets→doublets (Z₃→broken).
- The symmetry group determines the ordered-phase degeneracy pattern, but both look similar in disordered phase.

### 33c: BW Locality for Potts — S₃ vs Z₂ vs U(1)
**Status:** Complete

n=8, subsystem A = left 4 sites (dim_A = 81). Swept h/J from 0.05 to 2.0.

**Results:**
| h/J | Phase | Potts Locality | BW Variance Captured | Best Envelope |
|-----|-------|---------------|---------------------|--------------|
| 0.05 | Deep ordered | 51.1% | 49.9% | linear |
| 0.15 | Ordered | 65.8% | 65.1% | linear |
| 0.25 | Critical | **76.5%** (peak) | 75.8% | linear |
| 0.50 | Disordered | 59.2% | 58.8% | linear |
| 1.00 | Deep disordered | 45.9% | 44.7% | linear |

**KEY RESULT — Sprint 032's prediction is WRONG:**
| Model | Symmetry | |G| | d | Peak BW Locality |
|-------|----------|-----|---|-----------------|
| XXZ | U(1) | ∞ | 2 | **100.0%** |
| TFIM | Z₂ | 2 | 2 | **91%** |
| Potts | S₃ | 6 | 3 | **76.5%** |

S₃ has order 6 (between Z₂=2 and U(1)=∞), but gives the LOWEST BW locality!

**Why the prediction fails:**
1. BW locality = fraction of H_E captured by physical Hamiltonian terms
2. The relevant quantity is NOT |G| but the ratio of Hamiltonian terms to symmetry-allowed terms
3. For d=2 (qubit): operator space has 4^{n_A} = 256 terms. Z₂/U(1) strongly constrain which operators appear.
4. For d=3 (qutrit): operator space has 9^{n_A} = 6561 terms. Even S₃ (order 6) allows many more non-Hamiltonian operators.
5. S₃ allows additional 2-body invariants beyond δ and P+P†: e.g., P_i P_j + P_i† P_j†, P_i P_j† + P_i† P_j — these are S₃-invariant but NOT in the physical Hamiltonian.
6. U(1) is special: for d=2, it forces ALL 2-body terms to be XX+YY+ZZ type, which ARE the Hamiltonian terms.

**Corrected principle:** BW locality is controlled by the ratio dim(Hamiltonian operator space) / dim(symmetry-allowed operator space). Group size alone is insufficient; local Hilbert space dimension matters crucially.

**Coupling profiles confirm BW envelope:**
- δ coefficients decrease from far-from-cut to near-cut: 10.1 → 7.3 → 3.5 (h/J=0.10)
- P+P† coefficients similarly decrease: 0.80 → 0.58 → 0.54 → 0.16
- Entanglement temperature gradient IS Unruh-like (hot at cut, cold in bulk)
- BW picture qualitatively correct even though quantitatively weaker (76.5% vs 91-100%)
