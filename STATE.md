# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 111 — q=10 χ_F measurement. Result: **α(q=10)≈3.42 (last-pair), single-multiplet dominance universal through q=10.** Power-law α(q)=0.69·q^0.69 AIC-preferred over linear (ΔAIC=4.1), but discrimination weak at accessible q range.

## Active Research Thread
**Walking regime: six confirmed novel findings, α(q) extended to q=10.** Casimir-Re(c) (098), χ_F scaling α(q) (103/110/111), χ_F spectral mechanism (107), entropy concentration, Rényi mapping, multiplet dominance. Sprint 111 adds q=10 data point and power-law functional form.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 86 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unspent for 86 sprints. q=2 Ising χ_F peak or Casimir on ~10 qubits. Strongest prediction: χ_F ~ N^1.0 at g_c with known ν=1.
2. **q=4 BKT convergence** — DMRG at n=12-20 to test whether α continues toward 2.0. q=4 is stuck at 1.77.
3. **α(q) sublinear test at q=12-15** — Power-law (0.69q^0.69) vs linear (0.26q+0.83) diverge strongly at q≥12. q=12 n=6 (dim=2.99M, GPU feasible) would give α_linear=3.95 vs α_power=4.50 — but only 2 sizes available, so discrimination is limited.

## What's Been Ruled Out
- α(q)=0.315q+0.469 (Sprint 103) — biased by q=4 BKT inclusion
- α(q)=0.262q+0.815 (Sprint 110) — still valid but power-law 0.69q^0.69 slightly preferred (ΔAIC=4.1)
- Sprint 108's 1/ln(N) correction form for q=2,3 (retracted Sprint 109)
- Energy-entropy hierarchy universality (Sprint 104) — walking-specific
- χ_F super-scaling universality (Sprint 105) — walking-specific
- SREE as walking discriminator (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)

## Key Tools Available
- χ_F spectral: q=2 n≤18, q=3 n≤14, q=4 n≤11, q=5 n≤10, q=6 n≤9, q=7 n≤8, q=8 n≤7, q=9 n≤7, q=10 n≤7
- Vectorized S_q Potts builder (fast: q=6 n=9 build in 30s vs 184s with loop)
- Power-law correction fit (Sprint 109): recovers exact ν for q=2,3
- GPU eigsh limits: q=6 n=9 (10M, 305s), q=10 n=7 (10M, 700s)
- IBM QPU: 580s remaining
