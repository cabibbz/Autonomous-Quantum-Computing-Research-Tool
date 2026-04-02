# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 105 — χ_F scaling at J1-J2 BKT transition. Result: **Walking χ_F super-scaling is unique.** BKT gives α→0 (invisible at N≤20), MG first-order gives saturating peak, only Potts walking gives persistent α>2 with upward-converging pairwise exponents.

## Active Research Thread
**Walking regime: five confirmed novel findings + two scope constraints.** Casimir (098), χ_F scaling (103), entropy concentration, Rényi mapping, multiplet dominance. Sprint 104: energy-entropy hierarchy walking-specific. Sprint 105: χ_F super-scaling walking-specific. Both confirmed novel results now properly scoped to walking mechanism only.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 80 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unspent. Strongest prediction: q=2 Ising χ_F at g_c with exact α=0.98 (measurable on ~10 qubits via VQE).
2. **Walking χ_F mechanism** — WHY does walking give α>2? Is it the entanglement spectrum reorganization (multiplet dominance) that drives super-susceptibility? Test: decompose χ_F into spectral contributions.
3. **2D walking test** — Can we detect walking signatures in 2D S_q Potts? We have g_c for Ly=2,3 cylinders. Even a single 2D χ_F measurement would be novel.

## What's Been Ruled Out
- Energy-entropy hierarchy universality (Sprint 104) — walking-specific
- χ_F super-scaling universality (Sprint 105) — walking-specific, BKT invisible
- SREE as walking discriminator (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)

## Key Tools Available
- J1-J2 chain Sz-sector exact diag (Sprint 104-105): N≤20 in <100s
- Vectorized S_q Potts builder (Sprint 098/103)
- χ_F infrastructure: coupling/field split for fast g-scan
- GPU eigsh: q=5 n=10 (156s), q=7 n=8 (90s)
- S_q Potts DMRG: q=2 (fast, n≤24+), q=5 (n≤12, chi≤50)
- IBM QPU: 580s remaining
