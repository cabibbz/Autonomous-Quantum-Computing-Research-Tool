# Sprint 095 — DMRG BW Fidelity: Does the Threshold Persist When n≫nA?

**Date:** 2026-04-02
**Status:** In progress

## Motivation

Sprint 094 found BW R² collapses at a threshold nA* (≈5 for q=2-4, ≈4 for q=5), with exponential rate B(q) = 0.48q + 1.09. But those experiments always used n = 2·nA (equal bipartition). The BW theorem applies in the limit nA ≪ n. Literature (Dalmonte et al., Giudici et al.) shows lattice BW trace distance scales as ℓ⁻² for large ℓ — accuracy *improves* with subsystem size when n→∞.

**Key question:** Is the Sprint 094 threshold a genuine lattice effect, or an artifact of the equal-bipartition setup (nA/n = 1/2)?

**Approach:** Use DMRG at n=20-24 to extract full ρ_A for small subsystems (nA=3-6) embedded in a much larger chain. Compare BW R² at fixed nA but varying n/nA ratio. If R² improves dramatically at large n, the threshold is a finite-n artifact.

**Literature context:**
- Giudici et al. (SciPost 2018): lattice BW works well for 1D critical chains, corrections ~ 1/ℓ²
- Dalmonte et al. (Ann. Phys. 2022): review of BW on lattice — best for nA/n small
- Our Sprint 094: threshold at nA=5 for q=2 with n=2·nA — contradicts large-ℓ improvement

## Experiments

### 095a — q=2 DMRG BW R² at n=12,16,20,24 with nA=3-6
[Results pending]

### 095b — q=5 DMRG BW R² at n=12,16,20,24 with nA=3-5
[Results pending]

### 095c — Compilation: R²(nA, n/nA) dependence
[Results pending]
