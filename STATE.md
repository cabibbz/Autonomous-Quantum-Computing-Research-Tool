# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 103 — Harden χ_F Scaling + α(q) Mapping. q=5 confirmed with 4 sizes: α=2.091±0.002 (pairwise stable to 0.3%). New: q=6 α=2.37, q=7 α=2.65. α(q) = 0.315q + 0.469 linear for q≥4. Super-first-order (α>2) for all q≥5. Upgraded from POTENTIALLY NOVEL to CONFIRMED NOVEL.

## Active Research Thread
**Walking regime: five confirmed novel findings now accumulated.** Casimir (Sprint 098), χ_F scaling (Sprint 103), plus entropy concentration, Rényi mapping, and multiplet dominance. The energy-entropy hierarchy is the strongest single result — should test in a different model for universality.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 78 SPRINTS

## Top 3 Next Experiments
1. **Test energy-entropy hierarchy in J1-J2 chain** — Does the Casimir/entropy decoupling appear at the J1-J2 deconfined critical point? If yes, this is universal (not Potts-specific), upgrading finding to PRL level.
2. **χ_F DMRG at q=5** — Can DMRG fidelity extend to n=12-20? Would confirm α stays at 2.09 at larger sizes.
3. **Hardware validation** — 580s QPU. Strongest prediction: q=2 Ising fidelity at g_c on 5-10 qubits, compare with exact α=0.98.

## What's Been Ruled Out
- Entanglement asymmetry as walking probe (ΔS_A=0 for symmetric ground states) — Sprint 102
- SREE as walking discriminator at accessible sizes (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)

## Key Tools Available
- Vectorized S_q Potts builder (Sprint 098/103): fast for n≤10 q=5, n≤8 q=7
- χ_F infrastructure: coupling/field split for fast g-scan
- GPU eigsh: q=5 n=10 (156s), q=7 n=8 (90s)
- S_q Potts DMRG: q=2 (fast, n≤24+), q=5 (n≤12, chi≤50)
- IBM QPU: 580s remaining
