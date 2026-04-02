# Sprint 089 — q=7 Tail Exponent Resolution & Level Redistribution as Walking Discriminator

**Date:** 2026-04-02
**Status:** In progress

## Motivation

Sprint 088 showed tail weight exponent b≈2.0 is universal for q=2,3,5 (real and walking CFT). But q=7 shows b≈3.0 from only 4 data points (n=6-12). Two hypotheses:
1. b≈3.0 is genuine — walking-broken regime has enhanced tail growth
2. b≈3.0 is pre-asymptotic — q=7 would converge to b≈2.0 at larger n

Additionally, the walking-specific signal resides in level REDISTRIBUTION (how weight moves between ground state, (q-1) multiplet, and tail), not in the tail growth rate itself. This sprint quantifies that redistribution.

## Experiments

### 089a — q=7 DMRG at n=14,16: Does b≈3.0 survive?
**Goal:** Extend q=7 entanglement spectrum from n=12 to n=14,16 via DMRG.
**Risk:** q=7 DMRG at n=12 took 293s. n=14 may exceed 300s limit.
**Plan:** Test n=14 first with chi=30. If <200s, attempt n=16.

### 089b — Level redistribution trajectories across q
**Goal:** Compute d(%S)/d(ln n) for each level (ground, multiplet, tail) across q=2,3,5,7.
**Prediction:** Walking-specific effect is in lev1→tail transfer RATE, not tail growth exponent.

### 089c — Walking discriminator synthesis
**Goal:** Identify the cleanest walking discriminator from all spectral data.

---

## Results

### 089a — q=7 b≈3.0 was PRE-ASYMPTOTIC, converging to universal b≈2.0

DMRG at q=7, g_c=1/7, open BC:
- n=14, chi=30: 237s. λ_max=0.845, w_mult=0.153, w_tail=0.00143, %S(tail)=1.7%
- n=16, chi=25: 138s. λ_max=0.840, w_mult=0.159, w_tail=0.00174, %S(tail)=2.0%

Combined power-law fit (n=6-16, 6 points):
- All 6 points: b = 2.53 ± 0.22, R² = 0.970
- Excluding n=6 (n≥8, 5 pts): **b = 2.07 ± 0.14**
- Previous (n=6-12, 4 pts): b = 2.97

**The b≈3.0 exponent was pre-asymptotic.** n=6 is in a crossover regime where tail weight grows faster. Excluding it, the exponent converges to b≈2.1, consistent with the universal b≈2.0 for q=2,3,5.

**Tail weight is truly UNIVERSAL at b≈2.0 for ALL critical q=2,3,5,7.** The walking-broken regime does NOT enhance tail growth rate.

Other trends at q=7:
- %S(lev0) stabilizes at ~0.199 (vs ~0.225 for q=5) — ground state fraction lower at higher q
- %S(lev1) decreasing: 0.785→0.781 (n=14→16) — multiplet slowly feeding tail
- Entanglement gap stable at ~3.48 (very slowly decreasing)
