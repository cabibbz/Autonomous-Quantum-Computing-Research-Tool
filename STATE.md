# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 109 — α(q) correction form: power-law vs logarithmic. Result: **Sprint 108's 1/ln(N) extrapolation RETRACTED for q=2,3.** Power-law 1/N² corrections recover exact ν to 0.0-0.2%. Walking (q≥5) confirmed zero corrections. q=4 BKT slow convergence confirmed.

## Active Research Thread
**Walking regime: six confirmed novel findings.** Casimir-Re(c) (098), χ_F scaling α(q) (103), χ_F spectral mechanism (107), entropy concentration, Rényi mapping, multiplet dominance. Sprint 109 refines: α(q)=0.315q+0.469 confirmed as TRUE asymptotic formula for q≥5 (zero FSS corrections). q=2,3 exact ν recovered via power-law fits.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 84 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unspent. q=2 Ising χ_F peak or Casimir on ~10 qubits.
2. **q=4 BKT convergence** — DMRG at n=12-20 to test whether α continues drifting toward 2.0. Or: dedicated BKT log correction fit with known p exponent.
3. **Walking χ_F at q=6,7 GPU extension** — Extend q=6 to n=9-10, q=7 to n=8-9 to test if α(q)=0.315q+0.469 holds at more sizes.

## What's Been Ruled Out
- 1/ln(N) correction form for q=2,3 (Sprint 109) — power-law 1/N² is correct
- Sprint 108 extrapolated α_∞: q=3→1.22 WRONG (correct: 1.40), q=2 WRONG (correct: 1.00)
- Energy-entropy hierarchy universality (Sprint 104) — walking-specific
- χ_F super-scaling universality (Sprint 105) — walking-specific, BKT invisible
- SREE as walking discriminator (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)
- Entanglement gap as proxy for energy multiplet gap (Sprint 107)

## Key Tools Available
- χ_F spectral decomposition: q=2 n≤18, q=3 n≤14, q=4 n≤11, q=5 n≤10 (all GPU)
- Power-law correction fit infrastructure (Sprint 109): recovers exact ν
- Combined scaling: z_m(q) = 0.065q + 0.845, β_me(q) = 0.188q − 0.238
- J1-J2 chain Sz-sector exact diag (Sprint 104-105): N≤20 in <100s
- Vectorized S_q Potts builder (Sprint 098/103/107/109)
- GPU eigsh: q=3 n=14 (466s), q=4 n=11 (503s), q=5 n=10 (279s)
- IBM QPU: 580s remaining
