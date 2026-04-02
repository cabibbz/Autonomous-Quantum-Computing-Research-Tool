# Sprint 081 — ξ*(q=6) via DMRG: Is q=6 a Second Walking Case?

**Date:** 2026-04-02
**Status:** Complete (3 experiments)

## Motivation

Sprint 080 mapped the walking boundary: c_eff/Re(c) = exp(-0.105(q-5)) with q=5 as the exact threshold. q=6 is the marginal case — ratio 0.92 at n=8, and crucially c_eff was STABLE with n (+0.12% drift n=6→8). This stability is the hallmark of walking: the system hasn't yet "seen" its first-order nature at accessible sizes.

**Key question:** Does q=6 c_eff stay at ~1.15 out to n=12-24, or does it eventually drop?

## Literature Context

- Ma & He (PRB 99, 195130, 2019): Measured c_eff for q=5,6,7 S_q Potts, found approximate conformality but didn't track size-dependence systematically.
- Gorbenko, Rychkov & Zan (SciPost 2018): Complex CFT walking framework predicts walking length decreases with q.
- arXiv:2502.02001 (2025): Complex entanglement entropy for complex CFT — complementary non-Hermitian approach.

## Results

### 081a — DMRG c_eff(q=6) at n=10,12

DMRG with d=6 is very slow (n=8 chi=40 takes 172s). Redesigned to run n=10 (chi=40, 48s) and n=12 (chi=30, 205s).

**c_eff is DROPPING with system size:**

| n | method | chi | c_eff | R² | c/Re(c) | time |
|---|--------|-----|-------|----|---------|------|
| 6 | exact  | — | 1.146 | — | 0.915 | — |
| 7 | exact  | — | 1.148 | — | 0.917 | — |
| 8 | exact  | — | 1.147 | — | 0.916 | — |
| 10 | DMRG  | 40 | 1.130 | 0.99999 | 0.902 | 48s |
| 12 | DMRG  | 30 | 1.115 | 0.99998 | 0.890 | 205s |

c_eff drift n=8→12: -2.9%. DMRG at n=8 verified against exact diag (1.1473 vs 1.147).

### 081b — DMRG gap×N(q=6) at n=10

Excited state DMRG (orthogonal_to method) for energy gap.

| n | chi | gap×N | time |
|---|-----|-------|------|
| 10 | 40 | 2.101 | 313s |

gap×N = 2.10 at n=10, matching/exceeding Sprint 080's n=6 value (~2.01). Gap is healthy — breakdown is in entropy, not gap (confirms Sprint 080 finding).

n=12 infeasible: excited state DMRG takes ~2× ground state time.

### 081c — Three-way comparison q=5,6,7

**c_eff drift rate dc/d(ln n):**
- q=5: **+0.014** (rising — real CFT-like FSS overshoot)
- q=6: **-0.048** (dropping — walking degrading)
- q=7: **-0.091** (dropping — walking broken)

**c_eff/Re(c) at n=12:**
- q=5: 1.013 (above Re(c) — walking perfect)
- q=6: 0.890 (below but still >0.85)
- q=7: 0.784 (well below — walking broken)

**q=6 drift n=8→12 is CLOSER to q=7 (-4.6%) than q=5 (+1.0%):**
- Distance to q=5: 3.9 percentage points
- Distance to q=7: 1.8 percentage points

**Extrapolation:** c_eff crosses 0.85×Re(c) threshold at n ≈ 38.

## Key Findings

1. **q=6 walking is MARGINAL BREAKING, not stable.** c_eff drops 2.9% from n=8→12. The apparent stability at n=6-8 (Sprint 080) was because 6→8 is too small a range — the drop only becomes visible at n≥10.

2. **q=5 is truly unique as the walking sweet spot.** It's the only q where c_eff stays above Re(c) to n=20. q=6 is already 11% below Re(c) at n=12.

3. **Walking breakdown is hierarchical:**
   - q=5: ξ* > 20 (stable to largest sizes tested)
   - q=6: ξ* ≈ 38 (extrapolated) — walking extends further than q=7 but eventually breaks
   - q=7: ξ* ≈ 6 (ratio < 0.85 already at n=6)

4. **Gap remains healthy even as entropy degrades.** gap×N = 2.10 for q=6 at n=10. Consistent with Sprint 080: breakdown is in entanglement entropy, not energy spectrum.

5. **The walking boundary is not sharp.** q=5 (extends), q=6 (marginal), q=7 (broken) — there's a crossover region, not a phase boundary.

## Surprises

- q=6 c_eff stability at n=6-8 was MISLEADING — too small a range to detect -0.048 drift rate
- q=6 drift rate is 2× closer to q=7 than to q=5
- Extrapolated breakdown at n≈38 — just barely outside exact diag+DMRG reach
- DMRG at d=6 is dramatically slower than d=5 (n=8: 172s vs ~10s for q=5)

[Full data: results/sprint_081a_sq_dmrg_c_q6.json, results/sprint_081b_sq_dmrg_gap_q6.json, results/sprint_081c_walking_comparison.json]
