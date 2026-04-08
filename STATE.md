# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 120 — Hybrid q=4 chi_F spectral decomposition. **g_c(hybrid,q=4)=0.3933. z_m=1.004 — exactly marginal.** Hybrid alpha peaks at q=4 (1.55) then declines: alpha(q)=[1.40, **1.55**, 1.41, 0.95, 0.55] for q=[3,4,5,7,10]. S_q alpha monotonically increases. Models diverge starting at q=4 (12% gap). Pairwise alpha drifts downward (1.64→1.49), opposite to S_q q=4 expectation (upward toward 2.0). Different universality classes even at q=4.

## ⚠ CRITICAL: Model Identity Correction (April 2026 Audit)
**All experiments from Sprint 076 onward use the STANDARD S_q Potts model**, not the "Potts-clock hybrid." Sprints 119-120 return to the hybrid model for comparison. See KNOWLEDGE.md for full details.

## Active Research Thread
**chi_F spectral decomposition as universal diagnostic.** z_m crosses 1 at q=4 for the hybrid, marking the walking→continuous boundary. The S_q model maintains z_m>1 (walking) at all q tested. The z_m=1 crossing is the microscopic mechanism for the universality class split.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 95 SPRINTS

## Top 3 Next Experiments
1. **z_m(q) continuous fit** — with 5 hybrid data points (q=3,4,5,7,10), fit z_m(q) to find precise q_cross where z_m=1. Should be q_cross ≈ 3.5-4.0.
2. **S_q q=4 extended sizes** — extend Sprint 118 to n=4-11 (same as hybrid) to compare pairwise alpha drift. Literature predicts upward drift toward 2.0.
3. **Log correction fit at q=4** — fit chi_F = A·N^α·(ln N)^{-p} for both models. Extract p and compare.

## What's Been Ruled Out
- Hybrid chi_F super-scaling — hybrid alpha < 2 for all q≥5 (Sprint 119)
- alpha(q=4) = 2.0 for S_q — measured 1.77 at accessible sizes (Sprint 118), log corrections (literature)
- alpha(q) quadratic, linear, sqrt, power-law for S_q — dAIC>=7.7 (Sprint 117)
- DMRG open-BC chi_F — boundary effects (Sprint 113)
- "Our model is novel" for sprints 076+ — it's the standard S_q Potts (April 2026 audit)

## Key Tools Available
- chi_F spectral (periodic BC): q=2 n<=18, q=3 n<=14, q=4 n<=11, q=5 n<=10, q=6 n<=9, q=7 n<=8, q=8-9 n<=7, q=10 n<=7, q=12 n<=6, q=15 n<=5, q=20 n<=5, q=25 n<=5, q=30 n<=4
- Hybrid chi_F spectral: q=3 n<=12, q=4 n<=11, q=5 n<=10, q=7 n<=8, q=10 n<=7
- GPU eigsh practical limit: ~10M dim
- IBM QPU: 580s remaining
