# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 121 — z_m(q) fit & S_q q=4 extended sizes. **q_cross = 3.58 ± 0.04** (all 5 fit forms converge). S_q q=4: alpha=1.777, z_m=1.092 — firmly walking. Hybrid q=4: alpha=1.55, z_m=1.004 — barely walking. Alpha drift: S_q converges to 1.77, hybrid still decreasing (1.49 at n=10-11). 19% divergence at largest sizes.

## ⚠ CRITICAL: Model Identity Correction (April 2026 Audit)
**All experiments from Sprint 076 onward use the STANDARD S_q Potts model**, not the "Potts-clock hybrid." Sprints 119-121 return to the hybrid model for comparison. See KNOWLEDGE.md for full details.

## Active Research Thread
**chi_F spectral decomposition as universal diagnostic.** z_m crosses 1 at q_cross=3.58 for the hybrid, marking the walking→continuous boundary. The S_q model maintains z_m>1 (walking) at all q tested (1.09 at q=4, slowly drifting toward 1). The two models are in different universality classes for all q≥4.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 96 SPRINTS

## Top 3 Next Experiments
1. **Log correction fit at q=4** — fit S_q chi_F = A·N²·(ln N)^{-p} to extract p. Compare with Salas-Sokal predictions. Also fit hybrid to see if different correction form.
2. **QPU hardware test** — strongest prediction: q=2 Ising chi_F ~ N^1.0 at g_c=0.5. Test on 5-8 qubits. Budget: 580s available.
3. **Hybrid q=3.5 interpolation** — non-integer q test to directly verify q_cross prediction.

## What's Been Ruled Out
- Hybrid chi_F super-scaling — hybrid alpha < 2 for all q≥5 (Sprint 119)
- alpha(q=4) = 2.0 for S_q — measured 1.77 at accessible sizes (Sprint 118, 121), log corrections
- alpha(q) quadratic, linear, sqrt, power-law for S_q — dAIC>=7.7 (Sprint 117)
- DMRG open-BC chi_F — boundary effects (Sprint 113)
- "Our model is novel" for sprints 076+ — it's the standard S_q Potts (April 2026 audit)
- q_cross at integer q=4 — actually q_cross=3.58±0.04 (Sprint 121)

## Key Tools Available
- chi_F spectral (periodic BC): q=2 n<=18, q=3 n<=14, q=4 n<=11, q=5 n<=10, q=7 n<=8, q=10 n<=7
- Hybrid chi_F spectral: q=3 n<=12, q=4 n<=11, q=5 n<=10, q=7 n<=8, q=10 n<=7
- GPU eigsh practical limit: ~10M dim
- IBM QPU: 580s remaining
