# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 060 — OPE Coefficients: Structure Constants of the Novel q>4 Potts CFT

## Active Research Thread
**Complete characterization of 1D quantum Potts CFT.** Seven quantities now measured for q=2-10:

| q | g_c | c | ν | x₁ | c/x₁ | # fields | Tower | C_sse |
|---|-----|---|---|-----|------|----------|--------|-------|
| 2 | 0.250 | 0.500 | 1.00 | 0.125 | 4.0 | 1 | ✓ | 0.500 |
| 3 | 0.333 | 0.800 | 0.86 | 0.133 | 6.0 | 2 | ✓ | 0.540 |
| 4 | 0.392 | 1.000 | 0.82 | ~0.117 | 8.5 | 3 | ✓ | 0.464 |
| 5 | 0.441 | ~1.10 | 0.85 | ~0.101 | 10.8 | 4 | ✓ | 0.376 |
| 7 | 0.535 | ~1.30 | 0.97 | ~0.086 | 15.1 | 6 | ✓ | 0.272 |
| 10 | 0.684 | ~1.40 | 1.12 | ~0.083 | 16.9 | 9 | ✓ | ~0.21 |

**Sprint 060 key finding:** C_{sigma*,sigma,epsilon} peaks at q=3 (~0.54), then decreases monotonically as ~q^{-0.8}. First OPE measurement for Hermitian q>4 Potts.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Large-q limit x₁(q→∞) and C_sse(q→∞)** — Extend to q=15,20 at n=4. Does x₁ saturate? Does C_sse follow a simple power law? Connect to free boson predictions.
2. **QPU validation of Potts physics** — 580s budget. Measure order parameter or energy gap ratio for q=3 on hardware. First hardware test of Potts CFT.
3. **Full OPE coefficient matrix** — Extract C_{sigma^k,sigma^j,epsilon} for all k,j at q=5,7. Map out the complete structure constant matrix. Compare harmonic channel couplings.

## What's Been Ruled Out
- c/x₁ = 2q for q≥4 (Sprint 058)
- Approximate conformality (Sprint 059: genuine towers confirmed)
- Free boson k² for harmonic ratios at any finite q (Sprint 057)
- iDMRG for c extraction at criticality (Sprint 055)
- Analytic continuation of Potts CFT for q>4 (Sprint 056)
- Small-scale QEC active correction (Sprints 026-028)
- Data collapse for ν at n≤10 (Sprint 053)
- Naive Z_q charge detection for conjugate pairs (Sprint 060)
- Simple analytic formula for C_sse(q) (Sprint 060)

## Key Tools Available
- Exact diag with periodic BC: n≤12 (q=2), n≤10 (q=3,4), n≤8 (q=5), n≤6 (q=7,10)
- DMRG (TeNPy): finite chain, open BC only
- Energy gap Δ·N crossing: best method for g_c
- CFT spectrum + momentum: tower structure, operator content
- OPE from matrix elements: ratio method |<eps|Z|sigma>|/|<0|Z|sigma>|
- IBM QPU: 580s remaining
