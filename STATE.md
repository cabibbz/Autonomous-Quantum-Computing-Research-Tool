# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 113 — DMRG overlap χ_F: boundary conditions matter. **Open BC α is fundamentally different from periodic BC α for q≥4 (Δα up to 0.49).** DMRG overlap method validated (exact match at all sizes), but open-BC DMRG cannot extend α(q) to larger sizes. Periodic/open χ_F ratio diverges with N — boundary fidelity susceptibility enhanced at walking transitions.

## Active Research Thread
**Walking regime: six confirmed novel findings, α(q) extended to q=12.** Casimir-Re(c) (098), χ_F scaling α(q) (103/110/111/112), χ_F spectral mechanism (107), entropy concentration, Rényi mapping, multiplet dominance. Sprint 113: DMRG open-BC path to larger sizes ruled out.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 88 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unspent for 88 sprints. q=2 Ising χ_F peak or Casimir on ~10 qubits. Strongest prediction: χ_F ~ N^1.0 at g_c with known ν=1.
2. **q=4 BKT convergence via periodic DMRG (iDMRG)** — Open-BC DMRG ruled out (Sprint 113). Try iDMRG to get bulk χ_F density directly, bypassing boundary effects.
3. **α(q) at q=15** — Quadratic predicts 4.15 vs √q predicts 4.33. Exact diag at n=4,5 (dims 51k, 759k) feasible on CPU.

## What's Been Ruled Out
- DMRG open-BC χ_F for α(q) extension — boundary effects too large (Sprint 113)
- α(q) linear (0.260q+0.815 from Sprint 110) — ruled out by q=12 data (ΔAIC=12.2)
- α(q) power-law (0.69q^0.69 from Sprint 111) — also overpredicts q=12 (ΔAIC=5.5)
- α(q)=0.315q+0.469 (Sprint 103) — biased by q=4 BKT inclusion
- Sprint 108's 1/ln(N) correction form for q=2,3 (retracted Sprint 109)
- Energy-entropy hierarchy universality (Sprint 104) — walking-specific
- χ_F super-scaling universality (Sprint 105) — walking-specific
- SREE as walking discriminator (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)

## Key Tools Available
- χ_F spectral (periodic BC): q=2 n≤18, q=3 n≤14, q=4 n≤11, q=5 n≤10, q=6 n≤9, q=7 n≤8, q=8 n≤7, q=9 n≤7, q=10 n≤7, q=12 n≤6
- DMRG overlap χ_F (open BC): validated but useless for α(q) due to BC effects
- Vectorized S_q Potts builder (fast: q=6 n=9 build in 30s vs 184s with loop)
- GPU eigsh limits: q=6 n=9 (10M, 305s), q=10 n=7 (10M, 700s), q=12 n=6 (3M, 160s)
- IBM QPU: 580s remaining
