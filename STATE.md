# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 112 — q=12 χ_F measurement. Result: **α(q=12)≈3.73 (last-pair). Quadratic α(q)=-0.010q²+0.41q+0.30 strongly AIC-preferred (ΔAIC=12.2 over linear).** Single-multiplet dominance universal through q=12.

## Active Research Thread
**Walking regime: six confirmed novel findings, α(q) extended to q=12.** Casimir-Re(c) (098), χ_F scaling α(q) (103/110/111/112), χ_F spectral mechanism (107), entropy concentration, Rényi mapping, multiplet dominance. Sprint 112 rules out linear α(q), establishes sublinear growth.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 87 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unspent for 87 sprints. q=2 Ising χ_F peak or Casimir on ~10 qubits. Strongest prediction: χ_F ~ N^1.0 at g_c with known ν=1.
2. **q=4 BKT convergence** — DMRG at n=12-20 to test whether α continues toward 2.0. q=4 is stuck at 1.77.
3. **α(q) √q vs quadratic discrimination** — q=15 (dim=15^6=11.4M, borderline GPU) or q=20 (too large for exact diag). Quadratic predicts 4.15 vs √q predicts 4.33 at q=15.

## What's Been Ruled Out
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
- χ_F spectral: q=2 n≤18, q=3 n≤14, q=4 n≤11, q=5 n≤10, q=6 n≤9, q=7 n≤8, q=8 n≤7, q=9 n≤7, q=10 n≤7, q=12 n≤6
- Vectorized S_q Potts builder (fast: q=6 n=9 build in 30s vs 184s with loop)
- Power-law correction fit (Sprint 109): recovers exact ν for q=2,3
- GPU eigsh limits: q=6 n=9 (10M, 305s), q=10 n=7 (10M, 700s), q=12 n=6 (3M, 160s)
- IBM QPU: 580s remaining
