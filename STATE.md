# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 059 — Conformal Tower Analysis: Genuine CFT Confirmed for q>4 Potts

## Active Research Thread
**Complete characterization of 1D quantum Potts CFT.** Six quantities now measured for q=2-10:

| q | g_c | c | ν | x₁ | c/x₁ | # spin fields | Tower? |
|---|-----|---|---|-----|------|---------------|--------|
| 2 | 0.250 | 0.500 | 1.00 | 0.125 | 4.0 | 1 | ✓ |
| 3 | 0.333 | 0.800 | 0.86 | 0.133 | 6.0 | 2 | ✓ |
| 4 | 0.392 | 1.000 | 0.82 | ~0.117 | 8.5 | 3 | ✓ |
| 5 | 0.441 | ~1.10 | 0.85 | ~0.101 | 10.8 | 4 | ✓ |
| 7 | 0.535 | ~1.30 | 0.97 | ~0.086 | 15.1 | 6 | ✓ |
| 10 | 0.684 | ~1.40 | 1.12 | ~0.083 | 16.9 | 9 | ✓ |

**Sprint 059 key finding:** Momentum-resolved conformal towers confirm GENUINE CFT for all q. Descendants at spin ±1 with exact degeneracy 4 (q≥3) or 2 (q=2). First measurement for Hermitian q>4 Potts.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **QPU validation of g_c(q=3)=1/3** — 580s budget. Prepare VQE-style ansatz near g_c, measure order parameter. First hardware test of Potts physics.
2. **OPE coefficients from 3-point functions** — Use operator overlaps to extract OPE coefficients. Tang et al. (2024) measured 9 for non-Hermitian q=5; we could measure them for Hermitian q=5,7,10.
3. **Large-q limit x₁(q→∞)** — Extend to q=15,20 at n=4,6. Does x₁ saturate? Does the tower gap approach a universal value? Connect to free boson prediction.

## What's Been Ruled Out
- c/x₁ = 2q for q≥4 (Sprint 058)
- Approximate conformality (Sprint 059: genuine towers confirmed)
- Free boson k² for harmonic ratios at any finite q (Sprint 057)
- iDMRG for c extraction at criticality (Sprint 055)
- Analytic continuation of Potts CFT for q>4 (Sprint 056)
- Small-scale QEC active correction (Sprints 026-028)
- Data collapse for ν at n≤10 (Sprint 053)
- All Sprints 038-048 Potts g_c values: WRONG (Sprint 049)

## Key Tools Available
- Exact diag with periodic BC: n≤12 (q=2), n≤10 (q=3,4), n≤8 (q=5), n≤6 (q=7,10)
- DMRG (TeNPy): finite chain, open BC only
- Energy gap Δ·N crossing: best method for g_c
- CFT spectrum + momentum: tower structure, operator content
- IBM QPU: 580s remaining
