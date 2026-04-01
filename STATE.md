# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 058 — x₁(q) Scaling Dimensions: Peak at q=3, Sub-linear c/x₁ Growth

## Active Research Thread
**Complete characterization of 1D quantum Potts CFT.** Four quantities now measured for q=2-10:

| q | g_c | c | ν | x₁ | c/x₁ | # spin fields |
|---|-----|---|---|-----|------|---------------|
| 2 | 0.250 | 0.500 | 1.00 | 0.125 | 4.0 | 1 |
| 3 | 0.333 | 0.800 | 0.86 | 0.133 | 6.0 | 2 |
| 4 | 0.392 | 1.000 | 0.82 | ~0.117 | 8.5 | 3 |
| 5 | 0.441 | ~1.10 | 0.85 | ~0.101 | 10.8 | 4 |
| 7 | 0.535 | ~1.30 | 0.97 | ~0.086 | 15.1 | 6 |
| 10 | 0.684 | ~1.40 | 1.12 | ~0.083 | 16.9 | 9 |

**Key finding:** x₁ peaks at q=3, c/x₁ grows sub-linearly (NOT 2q for q≥4). Sign flip in finite-size corrections at q≥5.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Conformal tower analysis** — Identify descendant structure and verify L₀ levels. Compute L₋₁σ and check if descendant gaps match x₁+1 prediction. This would confirm genuine CFT vs approximate conformality.
2. **QPU validation of g_c(q=3)=1/3** — 580s budget. Prepare VQE-style ansatz near g_c, measure order parameter. Would be first hardware test of Potts physics.
3. **x₁(q) at q=4 with DMRG trick** — Use open-BC DMRG at large n (n=40+) with entanglement entropy to independently extract x₁ via S(l) = (c/3)ln(l) + boundary terms. Compare to periodic-chain result.

## What's Been Ruled Out
- c/x₁ = 2q for q≥4 (Sprint 058: q=10 data 16.9 ≠ 20)
- Simple formula for c/x₁(q) — no closed form fits (Sprint 058)
- Free boson k² for harmonic ratios at any finite q (Sprint 057)
- iDMRG S-vs-ln(xi) for c extraction at criticality (Sprint 055)
- Analytic continuation of Potts CFT for q>4 (Sprint 056)
- Small-scale QEC active correction (Sprints 026-028)
- Data collapse for ν at n≤10 (Sprint 053)
- All Sprints 038-048 Potts g_c values: WRONG (Sprint 049)

## Key Tools Available
- Exact diag with periodic BC: n≤12 (q=2), n≤10 (q=3,4), n≤8 (q=5), n≤6 (q=7,10)
- DMRG (TeNPy): finite chain, open BC only
- Energy gap Δ·N crossing: best method for g_c
- CFT spectrum on periodic chain: scaling dimension ratios from gap ratios
- IBM QPU: 580s remaining
