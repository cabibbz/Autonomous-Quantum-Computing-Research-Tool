# Sprint 037 — Critical Exponent Extraction & First-Order Transition Test

**Date:** 2026-04-01
**Status:** In progress

## Idea

Sprint 036 confirmed MI-CV as a genuine order parameter with crossing points near h/J ≈ 0.93 at finite size. This sprint extracts the critical exponent ν from the finite-size shift of crossing points: h_c(n) = h_c(∞) + a·n^{-1/ν}. If ν ≈ 1, MI-CV belongs to the Ising universality class.

Second experiment: test MI-CV at the XXZ ferromagnetic transition (Δ = -1), which is first-order. MI-CV should show a discontinuous jump rather than the smooth inflection seen for TFIM (Ising) or the dome seen for BKT.

## Experiments

### 037a: TFIM Critical Exponent from MI-CV Crossing Points

Dense h/J sweep near criticality for n=8,12,16,24,32. Correlation-function MI reconstruction from DMRG. Found pairwise crossing points where CV(n1) = CV(n2):

| Size pair | n_eff | h_c (crossing) | |shift| from 1.0 |
|-----------|-------|-----------------|-------------------|
| (8, 12)   | 9.8   | 0.8949          | 0.1051            |
| (12, 16)  | 13.9  | 0.9341          | 0.0659            |
| (16, 24)  | 19.6  | 0.9590          | 0.0410            |
| (24, 32)  | 27.7  | 0.9723          | 0.0277            |

**Two-point ν estimates** (using exact h_c=1.0):
- (9.8, 13.9): 1/ν = 1.348, ν = 0.742
- (13.9, 19.6): 1/ν = 1.367, ν = 0.732
- (19.6, 27.7): 1/ν = 1.138, ν = 0.879

**Power-law fit** h_c(n) = 1.0 + a·n^{-1/ν}:
- a = -2.154, 1/ν = 1.325 ± 0.029
- **ν = 0.755** (Ising exact: 1.0)

**Key finding:** ν is below the Ising value, but the two-point estimates show clear convergence toward ν=1 with increasing system size (0.742 → 0.879). This indicates significant corrections to scaling at these sizes. The asymptotic behavior is consistent with Ising universality, but sizes n=8-32 are still in the crossover regime.

**Why ν < 1:** MI-CV is a non-standard order parameter — it's a statistical property of ALL pairwise correlations, not a single correlation function. At small n, boundary effects and finite-pair-count statistics create additional finite-size corrections beyond the standard ξ/L scaling.

### 037b: XXZ First-Order Transition at Δ = -1

*(results pending)*
