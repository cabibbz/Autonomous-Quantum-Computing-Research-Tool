# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 110 — Extend α(q) to q=5-9. Result: **Sprint 103 formula α=0.315q+0.469 REVISED to α≈0.26q+0.81** (pure walking data, q=5-9). q=6 now has 4 sizes (was 2). Walking mechanism universal through q=9.

## Active Research Thread
**Walking regime: six confirmed novel findings, formula refined.** Casimir-Re(c) (098), χ_F scaling α(q) (103/110), χ_F spectral mechanism (107), entropy concentration, Rényi mapping, multiplet dominance. Sprint 110 refines α(q) coefficients using 5 walking-only q values.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 85 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unspent for 85 sprints. q=2 Ising χ_F peak or Casimir on ~10 qubits.
2. **q=4 BKT convergence** — DMRG at n=12-20 to test whether α continues toward 2.0. q=4 is stuck at 1.77.
3. **q=10 measurement** — Test α(q) at q=10 (n=6: dim=1M, n=7: dim=10M GPU). Predicted α≈3.41 (linear) or α≈3.28 (quadratic). Would distinguish linear from sublinear.

## What's Been Ruled Out
- α(q)=0.315q+0.469 (Sprint 103) — biased by q=4 BKT inclusion. Correct: α≈0.26q+0.81
- Sprint 107 component fits (z_m, β_me) — based on q=2-7 mixed. Walking-only: z_m=0.082q+0.741, β_me=0.098q+0.333
- 1/ln(N) correction form for q=2,3 (Sprint 109)
- Energy-entropy hierarchy universality (Sprint 104) — walking-specific
- χ_F super-scaling universality (Sprint 105) — walking-specific
- SREE as walking discriminator (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)

## Key Tools Available
- χ_F spectral: q=2 n≤18, q=3 n≤14, q=4 n≤11, q=5 n≤10, q=6 n≤9, q=7 n≤8, q=8 n≤7, q=9 n≤7
- Vectorized S_q Potts builder (fast: q=6 n=9 build in 30s vs 184s with loop)
- Power-law correction fit (Sprint 109): recovers exact ν for q=2,3
- GPU eigsh limits: q=6 n=9 (10M, 305s), q=9 n=7 (4.8M, 217s)
- IBM QPU: 580s remaining
