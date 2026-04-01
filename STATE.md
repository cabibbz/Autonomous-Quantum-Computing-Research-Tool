# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 053 — ν(q) at True Critical Points: Method Validated, Old ν Values WRONG

## Active Research Thread
**ν(q) extraction via corrected energy gap slopes.** Validated to 3% accuracy for q=3 (exact 5/6) and <1% for q=2 (exact 1).

| q | g_c | ν (corrected) | Confidence |
|---|-----|---------------|------------|
| 2 | 0.250 | 1.00 | High (3 size pairs, <1% from exact) |
| 3 | 0.333 | 0.86 | High (4 sizes, 3% from exact 5/6) |
| 4 | 0.392 | 0.82 | High (3 size pairs, all agree) |
| 5 | 0.441 | 0.85 | Medium (1 size pair only) |
| 7 | 0.535 | 0.97 | Low (1 pair, b-calibration may not transfer) |
| 10 | 0.684 | 1.12 | Low (1 pair (4,5), suspect) |

**Key insight:** ν(q=3,4,5) ≈ 0.82-0.86, nearly constant. Old non-monotonic picture was artifact.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **ν(q=7) with more sizes** — Need n=8 (7^8=5.7M, too large for exact diag). Try DMRG with penalty method for excited state, or use iDMRG for direct ν from correlation length.
2. **Central charge c(q) at true g_c** — Extract c from entropy scaling S=(c/6)ln(n) at correct g_c for q=3,4,5. Known: c=4/5 (q=3), c=1 (q=4). Tests whether our Potts model reproduces CFT predictions.
3. **QPU validation of g_c(q=3)=1/3** — Prepare q=3 Potts ground state on hardware at g≈0.33 and g≈0.50. Measure order parameter. Budget allows it.

## What's Been Ruled Out
- Small-scale QEC active correction (Sprints 026-028)
- Pauli fraction as BW metric (Sprint 035)
- Clock ≡ Potts for q≥4 (Sprint 041)
- Gell-Mann MI for q≥4 (Sprints 045-046)
- Raw MI-CV for d≥10 (Sprint 048)
- Entropy FSS for ν extraction (Sprint 049)
- All Sprints 038-048 Potts g_c values: WRONG (Sprint 049)
- Old g_c scaling law g_c ∝ (q-3)^{-0.85}: INVALIDATED
- Self-duality for q≥4: BROKEN (Sprint 050)
- **Data collapse for ν extraction at n≤10: UNRELIABLE (Sprint 053, 43% error)**
- **MI-CV data collapse for ν: UNRELIABLE (Sprint 053)**
- **All Sprint 045-047 ν values for q≥3: WRONG (Sprint 053)**

## Key Tools Available
- Exact diag: n≤10 (q=2), n≤8 (q=3,4), n≤6 (q=5,7), n≤5 (q=10)
- DMRG (TeNPy): ground state only (excited states failed for Potts)
- Energy gap Δ·N crossing: best method for g_c (exact diag, FSS-corrected)
- Corrected power-law slope: best method for ν (need ≥3 sizes, b=0.86)
- IBM QPU: 580s remaining
