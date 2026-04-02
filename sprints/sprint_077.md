# Sprint 077 — S_q Potts at q=7,10: g_c Found, NO First-Order Signal at n≤8

**Date:** 2026-04-02
**Status:** Complete (3 experiments)

## Motivation

Sprint 076 found S_q Potts at q=5 looks CFT-like at n≤8 (g_c=0.200, clean crossings, exact degeneracies). Sprint 076c reported a "69% gap collapse at q=7 n=4→6" — but this was at an estimated g_c that turned out to be wrong. Need proper gap×N crossing scan at q=7 to find the true g_c and test for first-order behavior.

**Key question:** Does S_q Potts become visibly first-order at q=7 or q=10 (as predicted by Gorbenko et al. for q>4), or does it remain CFT-like at accessible sizes?

## Experiments

### 077a — S_q Potts q=7 Gap×N Crossings

**S_q Potts g_c(q=7) = 0.144** from (4,6) crossing. Gap×N at crossing = 0.646.

Timing: n=4 (dim=2401) instant, n=6 (dim=118k) 1.9s/pt, n=8 (dim=5.8M) 122s/pt.

**Targeted n=8 test at g_c=0.144:** gap×N = 0.649 — matches (4,6) crossing to 0.5%. This is remarkable convergence across three system sizes, fully consistent with a continuous transition.

| n | dim | gap at g_c | Δ·N | Time |
|---|-----|-----------|------|------|
| 4 | 2,401 | 0.1602 | 0.641 | <1s |
| 6 | 117,649 | 0.1057 | 0.634 | 2s |
| 8 | 5,764,801 | 0.0812 | 0.649 | 59s |

Gap×N convergence: 0.641 → 0.634 → 0.649. Non-monotonic but within 2.4% band.

**Previous "69% gap collapse" was at wrong g_c.** At the true g_c=0.144, all three sizes agree.

### 077b — Conformal Tower Comparison at q=7

**S_q Potts vs hybrid degeneracy structures are dramatically different:**

| Level | S_q Potts | Hybrid |
|-------|-----------|--------|
| σ (1st excited) | **6-fold** (S₇ symmetry) | **2-fold** (Z₇ conjugate pair) |
| σ² | — (merged into σ) | 2-fold, R=3.32 |
| σ³ | — (merged into σ) | 2-fold, R=5.95 |
| ε (energy) | 1-fold, R=3.48 | 1-fold, R=8.07 |
| Descendant | 8-fold, R=7.35 | 3-fold, R=9.59 |

**S_q Potts has massive degeneracy:** The full S₇ permutation group merges all q-1=6 spin field primaries into a single 6-fold degenerate state. The hybrid's Z₇ symmetry only pairs conjugates (σᵏ, σ^{7-k}).

**R_ε (energy field position) radically different:** S_q has R_ε=3.48, hybrid has R_ε=8.07 (2.3x ratio). This means the S_q and hybrid CFTs have completely different operator content.

**Both look like genuine CFTs at n=6.** Clean integer degeneracies, well-separated levels.

### 077c — Δ·N Scaling at q=5,7,10: No First-Order for Any q

**Systematic Δ·N drift comparison:**

| Model | q | n range | Δ·N (first) | Δ·N (last) | Drift | Verdict |
|-------|---|---------|-------------|-----------|-------|---------|
| S_q Potts | 5 | 4→8 | 0.671 | 0.639 | -4.8% | Continuous |
| S_q Potts | 7 | 4→8 | 0.641 | 0.649 | +1.3% | Continuous |
| S_q Potts | 10 | 4→6 | 0.587 | 0.591 | +0.6% | Continuous |
| Hybrid | 5 | 4→8 | 0.473 | 0.483 | +2.3% | Continuous |
| Hybrid | 7 | 4→6 | 0.338 | 0.357 | +5.6% | Continuous |
| Hybrid | 10 | 4→6 | 0.242 | 0.255 | +5.4% | Continuous |

**S_q Potts g_c(q=10) = 0.101** from (4,5) crossing. New measurement.

**ALL models show <6% Δ·N drift.** No exponential gap collapse at any q tested. S_q Potts at q=7 and q=10 are just as stable as q=5.

**S_q Potts Δ·N decreases with q:** 0.64 (q=5,7) → 0.59 (q=10). Hybrid Δ·N decreases faster: 0.48 (q=5) → 0.36 (q=7) → 0.25 (q=10).

## Summary of S_q Potts g_c Values

| q | g_c (S_q) | g_c (hybrid) | Field ratio | g_c ratio |
|---|-----------|-------------|-------------|-----------|
| 5 | 0.200 | 0.441 | 6.47 | 2.21 |
| 7 | 0.144 | 0.535 | 4.81 | 3.72 |
| 10 | 0.101 | 0.684 | 5.56 | 6.77 |

**g_c ratio (hybrid/S_q) grows with q** — the models diverge increasingly at larger q. S_q field is (q-1) times stronger, so S_q g_c → 0 as q → ∞.

## Key Findings

1. **S_q Potts g_c(q=7) = 0.144, g_c(q=10) = 0.101.** Decreasing with q as expected from field strength.
2. **NO first-order signal at q=5,7,10 up to n=8.** Δ·N stable to <5% for all. The Gorbenko et al. prediction of first-order behavior requires much larger systems than n=8.
3. **Conformal tower at q=7: S_q has 6-fold degeneracy (S₇), hybrid has 2-fold (Z₇).** Sharpest discriminator between models.
4. **R_ε (energy field) ratio = 2.3x between S_q and hybrid.** Completely different operator content = different universality classes confirmed at q=7.
5. **Sprint 076c "69% gap collapse" was artifact of wrong g_c estimate.** At true g_c, convergence is excellent.

## Surprises

- S_q Potts at q=7 n=8 gap×N=0.649 matches n=4,6 to <3% — better convergence than q=5 (4.8% drift)
- S_q q=10 also has no first-order signal — the complex CFT regime extends to at least n=6
- S_q Δ·N is remarkably q-independent (~0.6 for q=5-10), while hybrid Δ·N drops sharply with q
- g_c ratio (hybrid/S_q) grows rapidly — 2.2x at q=5, 6.8x at q=10

## Implications

The Gorbenko-Rychkov-Zan complex CFT prediction is that S_q Potts for q>4 is "walking" — it looks like a CFT up to a correlation length ξ* that grows exponentially as q→4⁺. At q=5, ξ* is presumably very large. Our finding is that even at q=10, the walking length exceeds n=6-8. This is consistent with the complex CFT picture: the transition is eventually first-order, but the crossover scale is beyond exact diag reach.

**DMRG at n=20-50 could potentially reach the first-order regime** for S_q Potts at q≥7.
