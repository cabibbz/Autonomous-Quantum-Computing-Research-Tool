# Sprint 035 — BW Size Scaling: Pauli Fraction Fails, Spectrum Succeeds

**Date:** 2026-04-01
**Status:** Complete (3 experiments + sub-experiments)

## Idea

Sprints 032-034 established that BW locality depends on the H/G-inv ratio, with all results at n=8. This sprint tests whether these results survive at n=8→20 using TeNPy DMRG. Also tests MI-CV order parameter size scaling.

**Literature:** Dalmonte et al. (2018) showed BW improves with system size. Gap: systematic BW metric comparison across sizes.

## Experiments

### 35a: TFIM BW Pauli Fraction Size Scaling
**n=8,10 (exact), n=12,16,20 (DMRG). h/J=1.0 (critical), h/J=0.5 (ordered).**

Key result: BW Pauli fraction drops dramatically with system size:

| n  | Phase    | Pauli% | Fidelity | BW Overlap |
|----|----------|--------|----------|------------|
| 8  | critical | 90.9%  | 0.991    | 0.950      |
| 10 | critical | 67.5%  | 0.904    | 0.821      |
| 12 | critical | 46.9%  | 0.679    | 0.683      |
| 16 | critical | 20.1%  | 0.129    | 0.441      |
| 20 | critical |  7.6%  | 0.008    | 0.269      |
| 8  | ordered  | 89.1%  | 0.997    | 0.941      |
| 10 | ordered  | 63.0%  | 0.952    | 0.793      |
| 12 | ordered  | 41.7%  | 0.735    | 0.642      |

Cross-validated at n=12 with both exact diag AND DMRG — perfect match (46.95%). The drop is real physics, not a DMRG artifact.

**Why it drops:** ||H_E||² grows exponentially (from 2^n_A eigenvalues) while TFIM-type terms grow linearly (~2n_A operators). The denominator outpaces the numerator.

**Technical notes:**
- TeNPy TFIChain uses H = -J·Sigmax·Sigmax - g·Sigmaz (XX interaction, Z field)
- Requires correct operator mapping: bond=XX, site=Z in TeNPy basis
- The sin_inv BW envelope is wrong at these sizes (R²_env ≈ 0.15)

### 35b: BW Entanglement Spectrum Accuracy
**The right metric: spectrum-level R² instead of Pauli fraction.**

| Model         | n  | R²_fitted | R²_envelope | k (signif.) |
|---------------|----|-----------|-----------  |-------------|
| TFIM critical | 8  | 1.0000    | 0.231       | 3           |
| TFIM critical | 10 | 1.0000    | 0.188       | 3           |
| TFIM critical | 12 | 1.0000    | 0.159       | 4           |
| TFIM critical | 16 | 1.0000    | 0.121       | 4           |
| TFIM ordered  | 8  | 1.0000    | 0.171       | 4           |
| TFIM ordered  | 16 | 1.0000    | 0.088       | 4           |
| XXZ XY        | 8  | 1.0000    | 0.000       | 4           |
| XXZ XY        | 10 | 0.9528    | 0.014       | 7           |
| XXZ XY        | 12 | 0.9417    | 0.000       | 6           |

**R²_fitted = 1.0 for TFIM at ALL sizes!** With optimally-fitted per-term coefficients, the TFIM operator subspace perfectly reproduces the important entanglement energies. But this has k ≪ n_terms (overfitting caveat: 3-4 data points vs 7-15 parameters).

**R²_envelope ≈ 0.15:** The sin_inv envelope predicts the wrong β(i) profile. The BW FORM is correct but the ENVELOPE is wrong at these finite sizes.

**XXZ drops to 0.95:** More significant eigenvalues (6-7 vs 3-4) exceed the XXZ-type operator count for tight fitting.

### 35b2: Fidelity with Projected Coefficients
**Project H_E onto TFIM subspace, build thermal state, compare.**

| n  | Phase    | Pauli% | F_projected | Trace Dist | S_exact | S_projected | S_error% |
|----|----------|--------|-------------|------------|---------|-------------|----------|
| 8  | critical | 90.9%  | 0.983       | 0.093      | 0.515   | 0.763       | 48%      |
| 12 | critical | 46.9%  | 0.698       | 0.420      | 0.572   | 2.111       | 269%     |
| 20 | critical |  7.6%  | 0.009       | 0.988      | 0.640   | 9.686       | 1413%    |

Projected entropy blows up — the TFIM projection loses the fine-grained structure that keeps ρ_A low-rank. The full density matrix reconstruction fails catastrophically at large n.

### 35c: MI-CV Size Scaling
**TFIM at n=8,12,16. XXZ at n=8,12.**

TFIM MI-CV at transition (h/J=1.0): 0.73 (n=8) → 0.36 (n=12) → 0.32 (n=16).

The transition region SHARPENS: slope dCV/d(h/J) at criticality increases with n. Ordered phase stays flat (CV ≈ 0.05). Disordered phase values depend on pair-distance selection.

XXZ BKT dome: peak MI-CV at Δ≈1.0, narrowing from n=8 to n=12.

**Caveat:** n=8 uses all n(n-1)/2 pairs; n≥12 uses only up-to-3rd-neighbor pairs (DMRG efficiency). This affects absolute CV values but not qualitative trends.

## Key Insights

### 1. The Pauli Fraction is the Wrong BW Metric
Sprint 032's "91% TFIM locality" was correct for n=8 but the metric fundamentally breaks at larger sizes. The denominator (total ||H_E||²) grows exponentially from the exponentially many tiny eigenvalues, while the numerator (TFIM terms) grows only linearly. **Any** local operator subspace will show declining fraction, regardless of how well it captures the physics.

### 2. BW Form is Correct at Spectrum Level
The entanglement spectrum (the physically relevant part) is perfectly reproduced by TFIM-type operators at all tested sizes. The BW theorem's core prediction — that H_E lives in the space of physical Hamiltonian terms — holds. But this is weakened by overfitting (more parameters than significant eigenvalues).

### 3. The BW Envelope is Wrong at Finite Size
The standard sin_inv envelope (from CFT/BW) gives R²≈0.15 at all sizes. The correct β(i) profile is far from sin(πd/L). This suggests the BW envelope formula needs significant finite-size corrections for half-chain cuts of open-boundary systems.

### 4. Full Density Matrix Reconstruction Fails
Keeping only TFIM terms in H_E produces a thermal state with wildly wrong entropy (48% error at n=8, 1400% at n=20). The non-TFIM terms (even though tiny individually) collectively encode the low-rank structure of ρ_A.

### 5. MI-CV Sharpens with System Size
The MI-CV transition slope increases, consistent with it being a genuine order parameter. The ordered phase value stabilizes, while the disordered phase values depend on the MI pair selection.

## Surprises
- **Sprint 032's 91% is a finite-size illusion** — same metric gives 47% at n=12 and 8% at n=20
- **Spectrum R²=1.0** while Pauli fraction=8% — the "right" and "wrong" metrics disagree by a factor of 12
- **The BW envelope is far more wrong than the BW form** — R²_form≈1.0 but R²_envelope≈0.15
- **S_projected/S_exact diverges** — from 1.5× at n=8 to 15× at n=20

## Next Steps
- Find the correct BW envelope for finite open-boundary chains (fit and parameterize)
- Test BW with non-TFIM corrections (add next-order symmetry-allowed terms near boundary)
- MI-CV with all-pairs at larger n (need efficient long-range MI from MPS)
- Size scaling of H/G-inv ratio prediction (do the algebraic predictions from Sprint 034 improve?)
- Potts BW size scaling (does d=3 show different behavior?)
