# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 054 — Central Charge c(q) at True Critical Points

## Active Research Thread
**c(q) from entropy scaling at corrected g_c.** Validated for q=2,3 (matches CFT). q=4 has anomalous flat overshoot (logarithmic corrections). q=5 shows c>1 — potentially novel.

| q | g_c | ν | c (raw) | c (CFT) |
|---|-----|---|---------|---------|
| 2 | 0.250/1.0 | 1.00 | 0.516 | 0.500 |
| 3 | 0.333 | 0.86 | 0.884 | 0.800 |
| 4 | 0.392 | 0.82 | 1.229 | 1.000 |
| 5 | 0.441 | 0.85 | 1.335 | ? |

**Key insight:** c(q) increases with q. For q≤4 this matches the Potts CFT minimal model series. For q=5, c>1 places the transition outside minimal models.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **c(q=5) at larger n** — Need n=32,48 to confirm c>1 and measure convergence rate. Very expensive (d=5, ~8 min per point at chi=40). Try iDMRG instead.
2. **c(q=7,10) from iDMRG** — iDMRG gives c directly from transfer matrix eigenvalues, avoiding FSS entirely. Would extend c(q) trend and bypass slow convergence.
3. **QPU validation of g_c(q=3)=1/3** — Prepare q=3 Potts ground state on hardware. Budget allows it.

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
- Data collapse for ν at n≤10: UNRELIABLE (Sprint 053)
- All Sprint 045-047 ν values for q≥3: WRONG (Sprint 053)

## Key Tools Available
- Exact diag: n≤10 (q=2), n≤8 (q=3,4), n≤6 (q=5,7), n≤5 (q=10)
- DMRG (TeNPy): ground state only. Chi=20 sufficient for entropy at q=3.
- Energy gap Δ·N crossing: best method for g_c
- Corrected power-law slope: best method for ν
- Entropy scaling S=(c/6)ln(n): works for c, needs n≥32 for q≥4
- IBM QPU: 580s remaining
