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
[pending]

### 33c: BW Locality for Potts — S₃ vs Z₂ vs U(1)
[pending]
