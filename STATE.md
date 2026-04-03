# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 115 — q=20 χ_F: α=4.59, **logarithmic α(q) ≈ 1.78 ln(q) − 0.79 now AIC-BEST** (ΔAIC=9.4 over quadratic). Quadratic definitively ruled out (predicted 3.82, measured 4.59). z_m=2.05 crosses 2.0 for first time. Single-multiplet dominance confirmed through q=20.

## Active Research Thread
**Walking regime: six confirmed novel findings, α(q) extended to q=20.** Casimir-Re(c) (098), χ_F scaling α(q) (103/110-115), χ_F spectral mechanism (107), entropy concentration, Rényi mapping, multiplet dominance. Sprint 115: logarithmic α(q) confirmed as both AIC-best and physically preferred.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 90 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unspent for 90 sprints. q=2 Ising χ_F peak or Casimir on ~10 qubits. Strongest prediction: χ_F ~ N^1.0 at g_c with known ν=1.
2. **α(q) at q=25 or q=30** — Log predicts 4.93/5.26. Dim n=4: 390k/810k (GPU), n=5: 9.8M/24.3M (may exceed GPU). Would test continued logarithmic growth.
3. **q=4 BKT convergence via iDMRG** — Open-BC DMRG ruled out (Sprint 113). Try iDMRG for bulk χ_F density, bypassing boundary effects.

## What's Been Ruled Out
- α(q) quadratic — predicted 3.82 at q=20, measured 4.59 (Sprint 115). ΔAIC=+9.4
- α(q) linear (0.260q+0.815 from Sprint 110) — ruled out by q=12 data (Sprint 112)
- α(q) √q and power-law forms — ΔAIC≥11.3 with 9 pts (Sprint 115)
- DMRG open-BC χ_F for α(q) extension — boundary effects too large (Sprint 113)
- Sprint 108's 1/ln(N) correction form for q=2,3 (retracted Sprint 109)
- Energy-entropy hierarchy universality (Sprint 104) — walking-specific
- χ_F super-scaling universality (Sprint 105) — walking-specific
- SREE as walking discriminator (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)

## Key Tools Available
- χ_F spectral (periodic BC): q=2 n≤18, q=3 n≤14, q=4 n≤11, q=5 n≤10, q=6 n≤9, q=7 n≤8, q=8 n≤7, q=9 n≤7, q=10 n≤7, q=12 n≤6, q=15 n≤5, q=20 n≤5
- DMRG overlap χ_F (open BC): validated but useless for α(q) due to BC effects
- Vectorized S_q Potts builder (fast: q=6 n=9 build in 30s vs 184s with loop)
- GPU eigsh limits: q=6 n=9 (10M, 305s), q=10 n=7 (10M, 700s), q=12 n=6 (3M, 160s), q=15 n=5 (759k, 56s), q=20 n=5 (3.2M, 333s)
- IBM QPU: 580s remaining
