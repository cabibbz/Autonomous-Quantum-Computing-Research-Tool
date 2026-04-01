# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 052 — g_c(q) Scaling Law: √(q-1) Formula + q=10 Verification

## Active Research Thread
**g_c(q) ≈ (1/5)√(q-1) + 1/20 for 1D quantum Potts.** Fits 6 data points (q=2-10) to <5% error. Potentially novel — no prior measurement found.

| q | g_c (corrected) | Formula prediction | Error |
|---|-----------------|-------------------|-------|
| 2 | 0.250 (exact) | 0.250 | 0.0% |
| 3 | 0.333 (exact) | 0.333 | 0.05% |
| 4 | 0.392 | 0.396 | 1.2% |
| 5 | 0.441 | 0.450 | 2.1% |
| 7 | 0.535 | 0.540 | 0.9% |
| 10 | 0.684 | 0.650 | 4.9% |

**FSS corrections are size-pair dependent:** n=4,6 → 4.8%, n=6,8 → 2.5%.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **ν(q=3) at true g_c=1/3** — Data collapse near g=0.33. Known exact ν=5/6. Validates ν extraction at correct g_c.
2. **g_c(q) formula verification at q=15 or q=20** — Energy gap at q=15 (dim 15^4=50k at n=4, 15^6=11M at n=6). Tests formula prediction g_c(15)≈0.80.
3. **ν(q=5) recheck** — Old ν≈2.0 measured near g≈0.45, close to true g_c≈0.44. Quick validation with energy gap + MI-CV at correct g_c.

## What's Been Ruled Out
- Small-scale QEC active correction (Sprints 026-028)
- Pauli fraction as BW metric (Sprint 035)
- Clock ≡ Potts for q≥4 (Sprint 041)
- Gell-Mann MI for q≥4 (Sprints 045-046)
- Raw MI-CV for d≥10 (Sprint 048)
- Entropy FSS for ν extraction (Sprint 049)
- All Sprints 038-048 Potts g_c values: WRONG (Sprint 049)
- Old g_c scaling law g_c ∝ (q-3)^{-0.85}: INVALIDATED — g_c increases
- Self-duality for q≥4: BROKEN (Sprint 050)
- Logarithmic g_c(q) growth: too slow, √(q-1) is correct (Sprint 052)

## Key Tools Available
- Exact diag: n≤10 (q=2), n≤8 (q=3,4), n≤6 (q=5,7,10), n≤4 (q=15,20)
- DMRG (TeNPy): n>10, 1D only. PottsChain has NO conservation (slow).
- Energy gap Δ·N crossing: best method for g_c (exact diag, FSS-corrected)
- Direct MPS contraction MI: all-pairs MI for ANY d
- IBM QPU: 580s remaining
