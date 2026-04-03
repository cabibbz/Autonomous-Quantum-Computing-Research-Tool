# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 114 — q=15 χ_F: α=3.92, ALL prior models overpredict. Quadratic AIC-best but peaks at q≈17 (unphysical). **Logarithmic α(q) ≈ 1.73 ln(q) − 0.68 emerging as true functional form.** Components z_m and β_me both logarithmic in q. Single-multiplet dominance confirmed through q=15.

## Active Research Thread
**Walking regime: six confirmed novel findings, α(q) extended to q=15.** Casimir-Re(c) (098), χ_F scaling α(q) (103/110-114), χ_F spectral mechanism (107), entropy concentration, Rényi mapping, multiplet dominance. Sprint 114: logarithmic α(q) growth preferred over polynomial.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 89 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unspent for 89 sprints. q=2 Ising χ_F peak or Casimir on ~10 qubits. Strongest prediction: χ_F ~ N^1.0 at g_c with known ν=1.
2. **α(q) at q=20** — Log predicts 4.49, quadratic predicts 3.83. Dim n=4 = 160k (CPU), n=5 = 3.2M (GPU). Would decisively test log vs quadratic.
3. **q=4 BKT convergence via iDMRG** — Open-BC DMRG ruled out (Sprint 113). Try iDMRG for bulk χ_F density, bypassing boundary effects.

## What's Been Ruled Out
- α(q) linear (0.260q+0.815 from Sprint 110) — ruled out by q=12 data (ΔAIC=33.5 with 8 pts)
- α(q) √q and power-law forms — ruled out by q=15 (ΔAIC≥26)
- α(q) quadratic as physical form — peaks at q≈17, unphysical (but best local AIC)
- DMRG open-BC χ_F for α(q) extension — boundary effects too large (Sprint 113)
- Sprint 108's 1/ln(N) correction form for q=2,3 (retracted Sprint 109)
- Energy-entropy hierarchy universality (Sprint 104) — walking-specific
- χ_F super-scaling universality (Sprint 105) — walking-specific
- SREE as walking discriminator (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)

## Key Tools Available
- χ_F spectral (periodic BC): q=2 n≤18, q=3 n≤14, q=4 n≤11, q=5 n≤10, q=6 n≤9, q=7 n≤8, q=8 n≤7, q=9 n≤7, q=10 n≤7, q=12 n≤6, q=15 n≤5
- DMRG overlap χ_F (open BC): validated but useless for α(q) due to BC effects
- Vectorized S_q Potts builder (fast: q=6 n=9 build in 30s vs 184s with loop)
- GPU eigsh limits: q=6 n=9 (10M, 305s), q=10 n=7 (10M, 700s), q=12 n=6 (3M, 160s), q=15 n=5 (759k, 56s)
- IBM QPU: 580s remaining
