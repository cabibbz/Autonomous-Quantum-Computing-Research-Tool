# Sprint 047 — ν(q=7) Extraction: Does ν Decrease Past q=5?

**Status:** In progress
**Date:** 2026-04-01

## Motivation

The ν(q) picture so far: q=2 (1.0) → q=3 (5/6) → q=4 (≥2.2, BKT-like) → q=5 (2.0). The key question: does ν decrease for q>5, confirming q=4 as a peak? Or does it stay large?

If ν(q=7) < ν(q=5) ≈ 2.0, it confirms q=4 is a BKT-like peak separating standard 2nd-order (q≤3) from large-ν 2nd-order (q≥5), with ν eventually decreasing back toward mean-field values.

If ν(q=7) ≥ 2.0, the large-ν regime persists and q=4 may not be special — the entire q≥4 regime might have anomalously large ν.

## Method

- Direct MPS contraction for all-pairs MI (the only reliable method for d≥4)
- χ=20 (required for d≥7, Sprint 044)
- g_c(7) = 0.259 from Sprint 044 blind prediction (1.6% accurate)
- Dense sweep around g_c for n=8, 12, 16
- ν extraction: slope ratio dCV/dg at g_c across sizes, plus data collapse

## Experiments

### Exp 047a — q=7 n=8 timing test + MI-CV sweep

**Timing:** DMRG ~60-74s per g-point at d=7, n=8, χ=20. MI contraction ~2s. Total 743s for 11 points.

**Results (n=8):**
| g | CV |
|---|-----|
| 0.150 | 0.5973 |
| 0.200 | 0.5115 |
| 0.220 | 0.5189 |
| 0.240 | 0.5258 |
| 0.250 | 0.5292 |
| 0.260 | 0.4252 |
| 0.270 | 0.4288 |
| 0.280 | 0.4319 |
| 0.300 | 0.4378 |
| 0.350 | 0.4467 |
| 0.400 | 0.3070 |

**Key features:**
- CV dip at g=0.200 (local minimum), then rises to peak at g=0.250 (0.529)
- **Sharp drop at g=0.260 (0.425)** — transition clearly between g=0.25 and g=0.26, consistent with g_c=0.259
- Plateau at 0.425-0.447 for g=0.26-0.35 (disordered phase)
- Second drop to 0.307 at g=0.40 (deep disordered phase, approaching product state)
- NOT the same profile as q=4 (Sprint 046): q=7 shows a clear CV peak/drop, while q=4 had monotonically ordered curves

### Exp 047b — q=7 n=12 MI-CV sweep
**Status:** Running...
