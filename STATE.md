# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 051 — Energy Gap Method: g_c INCREASES with q, Reversing Old Picture

## Active Research Thread
**Energy gap Δ·N crossing gives true g_c for 1D quantum Potts.** Validated on q=2,3 (known exact). Applied to q=4,5,7:

| q | g_c | Method | Accuracy |
|---|-----|--------|----------|
| 2 | 0.250 | Self-duality (exact) | exact |
| 3 | 0.333 | Self-duality (exact) | exact |
| 4 | ~0.39 | Gap crossing n=6,8 + 2.5% correction | ~3% |
| 5 | ~0.44 | Gap crossing n=6,8 + 2.5% correction | ~3% |
| 7 | ~0.52 | Gap crossing n=4,6 + 2.5% correction | ~5% |

**g_c INCREASES with q.** Physical: q-fold degeneracy needs stronger X+X† field. This reverses 10 sprints of wrong results. q=5 old estimate (0.41) was closest to truth (0.44).

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **ν(q=3) at true g_c=1/3** — MI-CV data collapse with DMRG n=8,12,16 near g=0.33. Known exact: ν=5/6. First validation of ν extraction at correct critical point.
2. **ν(q=5) check** — Old ν≈2.0 was measured near g≈0.45, close to true g_c≈0.44. Quick validation: does data collapse still work at g_c=0.44?
3. **g_c(q) scaling law — NEW** — Fit g_c(q) for q=2-7. Sublinear growth. Can we predict g_c(10)?

## What's Been Ruled Out
- Small-scale QEC active correction (Sprints 026-028)
- Pauli fraction as BW metric (Sprint 035)
- Clock ≡ Potts for q≥4 (Sprint 041)
- Gell-Mann MI for q≥4 (Sprints 045-046)
- Raw MI-CV for d≥10 (Sprint 048)
- Entropy FSS for ν extraction (Sprint 049)
- **All Sprints 038-048 Potts g_c values: WRONG (Sprint 049)**
- **Old g_c scaling law g_c(q) = 0.87(q-3)^{-0.85}: INVALIDATED — g_c increases, not decreases**
- **Self-duality for q≥4: BROKEN (Sprint 050)**
- **Central charge from n≤12 for q=3 Potts: massive overshoot (Sprint 050)**

## Key Tools Available
- Exact diag: n≤10 (q=2), n≤8 (q=3,4), n≤6 (q=5,7), n≤4 (q=10)
- DMRG (TeNPy): n>10, 1D only. PottsChain has NO conservation (slow).
- **Energy gap Δ·N crossing: best method for g_c** (exact diag, ~2-3% error, no MI issues)
- Direct MPS contraction MI: all-pairs MI for ANY d
- IBM QPU: 580s remaining
