# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 106 — χ_F spectral decomposition. Result: **Walking super-scaling mechanism identified.** A single (q-1)-fold multiplet captures 100% of χ_F. α = β_me + 2z_m - 1, where z_m (multiplet gap exponent) and β_me (matrix element growth) are both linear in q. Spectral gap is symmetry-forbidden from contributing.

## Active Research Thread
**Walking regime: five confirmed novel findings + mechanism.** Casimir (098), χ_F scaling (103), entropy concentration, Rényi mapping, multiplet dominance. Sprint 104-105: scope constraints. Sprint 106: mechanism = accelerated multiplet gap closing + matrix element growth.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 81 SPRINTS

## Top 3 Next Experiments
1. **Harden χ_F mechanism** — The z_m and β_me decomposition is potentially novel. Need more sizes at q=5 (currently 3 sizes). Use GPU for n=10 (dim=9.8M, ~160s). Also cross-check with finite-difference χ_F from Sprint 103.
2. **Hardware validation** — 580s QPU unspent. Could test q=2 Ising χ_F peak at g_c (measurable on ~10 qubits).
3. **Connect mechanism to entanglement spectrum** — The (q-1)-fold multiplet in χ_F is the SAME as the (q-1)-fold entanglement spectrum multiplet from Sprint 084. Are z_m and the entanglement gap Δξ related?

## What's Been Ruled Out
- Energy-entropy hierarchy universality (Sprint 104) — walking-specific
- χ_F super-scaling universality (Sprint 105) — walking-specific, BKT invisible
- SREE as walking discriminator (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)

## Key Tools Available
- χ_F spectral decomposition infrastructure (Sprint 106): q≤7, n≤9 (GPU)
- J1-J2 chain Sz-sector exact diag (Sprint 104-105): N≤20 in <100s
- Vectorized S_q Potts builder (Sprint 098/103)
- GPU eigsh: q=5 n=10 (156s), q=7 n=8 (90s)
- S_q Potts DMRG: q=2 (fast, n≤24+), q=5 (n≤12, chi≤50)
- IBM QPU: 580s remaining
