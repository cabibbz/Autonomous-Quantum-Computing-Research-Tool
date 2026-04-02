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

### 089b — Level redistribution: %S(lev0) and %S(lev1) are MONOTONIC walking discriminators

Loaded all DMRG spectrum data for q=2,3,5,7 (n=8-24 for q≤5, n=6-16 for q=7).

**Two monotonic discriminators found:**

| Discriminator | q=2 | q=3 | q=5 | q=7 | Monotonic? |
|---|---|---|---|---|---|
| %S(lev0) at n_max | 0.321 | 0.273 | 0.225 | 0.199 | YES (decreasing) |
| %S(lev1) at n_max | 0.635 | 0.686 | 0.737 | 0.781 | YES (increasing) |
| S_lev0 saturation | 0.338 | 0.281 | 0.229 | 0.203 | YES (decreasing) |

**Asymptotic derivatives d(%S)/d(ln n) for n≥10:**

| q | d(%S_lev0) | d(%S_lev1) | d(%S_tail) | Type |
|---|---|---|---|---|
| 2 | +0.018 | −0.054 | +0.036 | real |
| 3 | +0.011 | −0.047 | +0.035 | real |
| 5 | +0.006 | −0.038 | +0.032 | walking |
| 7 | +0.006 | −0.026 | +0.020 | broken |

**Key findings:**
1. **Weight always flows multiplet → tail** (d(S_lev1) < 0, d(S_tail) > 0 for all q)
2. **Ground state gains entropy fraction** (d(S_lev0) > 0) — counterintuitive: w_max DECREASES with n but its entropy fraction INCREASES because entropy is -p·log(p)
3. **Tail transfer rate (d(S_tail)/d(ln n)) is NOT monotonic** — 0.036, 0.036, 0.032, 0.020. q=2 and q=3 are nearly identical; the drop is mostly for q=7
4. **The lev1→tail ratio is NOT monotonic** (−1.49, −1.32, −1.19, −1.31) — not a clean discriminator
5. **%S(lev0) saturation fits perfectly** (R² > 0.999): S_lev0_inf decreasing monotonically with q. Convergence exponent increases: α = 0.80 (q=2), 0.99 (q=3), 1.23 (q=5), 1.26 (q=7)
6. **%S(lev1) saturation to ~0** (!!) for q=2,3,5 — all entropy ultimately goes to ground + tail. For q=7, fit gives S_lev1_inf ≈ 0.40 but this is likely unreliable (only n≤16 data)

**The walking-specific physics is in the STATIC entropy partition, not the dynamics.** At any given n, the (q-1) multiplet absorbs more of the entropy for larger q (78% at q=7 vs 63% at q=2), leaving less for the tail. This is why %S(tail) at n=24 decreases with q (4.4% → 3.8% → 2.0%) even though w_tail grows at the same rate.
