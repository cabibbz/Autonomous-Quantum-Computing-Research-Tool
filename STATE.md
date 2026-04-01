# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 055 — Entropy Profile Method & c(q=5) > 1 Confirmed

## Active Research Thread
**Central charge c(q) for 1D quantum Potts.** Two methods now validated: FSS (Sprint 054) and entropy profile (Sprint 055). c(q=5) ≈ 1.10 ± 0.10, confirmed > 1 by both methods. c(q=4) has anomalous flat overshoot (log corrections). iDMRG failed at criticality.

| q | g_c | ν | c (profile, best) | c (exact/CFT) |
|---|-----|---|-------------------|---------------|
| 2 | 0.250/1.0 | 1.00 | 0.512 (n=64) | 0.500 |
| 3 | 0.333 | 0.86 | 0.827 (n=48) | 0.800 |
| 4 | 0.392 | 0.82 | 1.148 (n=24) | 1.000 |
| 5 | 0.441 | 0.85 | 1.261 (n=16) | ~1.10±0.10 |

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **c(q) formula** — Fit c(q) = f(q) using q=2-5 data. Predict c(q=7,10). Test with larger n when feasible.
2. **QPU validation of g_c(q=3)=1/3** — Prepare ground state on hardware near g_c. Budget allows. Strongest prediction available.
3. **Operator content at q=5 criticality** — If c>1, what CFT describes this? Measure scaling dimensions from energy spectrum. Gap ratios E_n/E_1 identify the CFT.

## What's Been Ruled Out
- iDMRG S-vs-ln(xi) for c extraction at criticality (Sprint 055)
- Small-scale QEC active correction (Sprints 026-028)
- Pauli fraction as BW metric (Sprint 035)
- Clock ≡ Potts for q≥4 (Sprint 041)
- Gell-Mann MI for q≥4 (Sprints 045-046)
- Raw MI-CV for d≥10 (Sprint 048)
- Entropy FSS for ν extraction (Sprint 049)
- All Sprints 038-048 Potts g_c values: WRONG (Sprint 049)
- Old g_c scaling law g_c ∝ (q-3)^{-0.85}: INVALIDATED
- Data collapse for ν at n≤10: UNRELIABLE (Sprint 053)
- q=7+ DMRG profile method: too slow at d≥7 (Sprint 055)

## Key Tools Available
- Exact diag: n≤10 (q=2), n≤8 (q=3,4), n≤6 (q=5,7)
- DMRG (TeNPy): finite only (iDMRG unreliable at criticality)
- Energy gap Δ·N crossing: best method for g_c
- Corrected power-law slope: best method for ν
- Entropy profile S(l) vs chord distance: best method for c (new Sprint 055)
- IBM QPU: 580s remaining
