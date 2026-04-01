# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 056 — c(q) Formula: Logarithmic Growth Confirmed, Analytic Continuation Ruled Out

## Active Research Thread
**Central charge c(q) for 1D quantum Potts.** c grows monotonically and approximately as ln(q). Analytic continuation of Potts CFT (Coulomb gas/minimal model) WRONG for q>4 — gives c<1 while we measure c>1. Quadratic interpolation also ruled out by q=7,10 measurements.

| q | g_c | ν | c (raw, best n) | c (corrected) | c (exact) |
|---|-----|---|-----------------|---------------|-----------|
| 2 | 0.250/1.0 | 1.00 | 0.512 (n=64) | 0.500 | 0.500 |
| 3 | 0.333 | 0.86 | 0.827 (n=48) | 0.803 | 0.800 |
| 4 | 0.392 | 0.82 | 1.148 (n=24) | ~1.00 | 1.000 |
| 5 | 0.441 | 0.85 | 1.261 (n=16) | ~1.10 | — |
| 7 | 0.535 | 0.97 | 1.462 (n=8) | ~1.3 | — |
| 10 | 0.684 | 1.12 | 1.596 (n=6) | ~1.4 | — |

Best fit: c ≈ 0.40·ln(q-1) + 0.55

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Operator content / scaling dimensions at q=5** — Energy spectrum gap ratios E_n/E_1 from exact diag identify the CFT. If c>1, is it a free boson + something? Or a non-unitary CFT? This would explain WHY c grows logarithmically.
2. **c(q=4) with larger n** — q=4 has anomalous flat overshoot from log corrections. Need n=32+ to see convergence. DMRG feasible (d=4 is fast). Would calibrate overshoot model for q>4.
3. **QPU validation of g_c(q=3)=1/3** — Prepare ground state near g_c on hardware. 580s budget available. Strongest prediction.

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
- Analytic continuation of Potts CFT for q>4: WRONG (Sprint 056)
- Quadratic c(q) interpolation: WRONG — fails at q=7,10 (Sprint 056)

## Key Tools Available
- Exact diag: n≤10 (q=2), n≤8 (q=3,4), n≤6 (q=5,7), n≤6 (q=10, 42s)
- DMRG (TeNPy): finite only (iDMRG unreliable at criticality)
- Energy gap Δ·N crossing: best method for g_c
- Corrected power-law slope: best method for ν
- Entropy profile S(l) vs chord distance: best method for c
- IBM QPU: 580s remaining
