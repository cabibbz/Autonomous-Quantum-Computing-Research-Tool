# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 108 — q=4 BKT crossover spectral decomposition. Result: **Three-regime log correction structure found.** Walking (q≥5) has zero log correction to α, continuous (q=3) has strong +0.46/ln(N), BKT (q=4) moderate +0.24/ln(N). q=4 α_∞ ≈ 1.66 (not 1.77). Gap ratio (mult/spec) decreases monotonically with q: 6.2→5.2→4.6.

## Active Research Thread
**Walking regime: six confirmed novel findings.** Casimir-Re(c) (098), χ_F scaling α(q) (103), χ_F spectral mechanism (107), entropy concentration, Rényi mapping, multiplet dominance. All scoped as walking-specific (104-105). Sprint 108 adds: linear formula α=0.315q+0.469 valid ONLY for walking (q≥5); q≤4 inflated by log corrections.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 83 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unspent. q=2 Ising χ_F peak at g_c measurable on ~10 qubits. Or test Casimir scaling.
2. **H_E operator content at walking** — Sprint 097 showed H_E becomes non-compact at threshold. Does the energy multiplet (which dominates χ_F) have a signature in H_E?
3. **DMRG extension of log correction test** — q=3 at n=12-24 (DMRG) could nail down α_∞ independently. q=5 DMRG at n=12-16 could confirm zero log correction.

## What's Been Ruled Out
- Energy-entropy hierarchy universality (Sprint 104) — walking-specific
- χ_F super-scaling universality (Sprint 105) — walking-specific, BKT invisible
- SREE as walking discriminator (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)
- Entanglement gap as proxy for energy multiplet gap (Sprint 107) — different z exponents
- Linear α(q) formula for q≤4 (Sprint 108) — log corrections invalidate finite-size values

## Key Tools Available
- χ_F spectral decomposition infrastructure (Sprint 106-108): q≤7, n≤10 (GPU)
- Log correction fits (Sprint 108): α_∞ extrapolation via 1/ln(N)
- Combined scaling: z_m(q) = 0.065q + 0.845, β_me(q) = 0.188q − 0.238
- J1-J2 chain Sz-sector exact diag (Sprint 104-105): N≤20 in <100s
- Vectorized S_q Potts builder (Sprint 098/103/107)
- GPU eigsh: q=5 n=10 (279s), q=4 n=10 (96s), q=7 n=8 (185s)
- S_q Potts DMRG: q=2 (fast, n≤24+), q=5 (n≤12, chi≤50)
- IBM QPU: 580s remaining
