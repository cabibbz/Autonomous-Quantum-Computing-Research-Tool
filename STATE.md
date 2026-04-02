# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 104 — Energy-Entropy Hierarchy Universality Test. Tested in J1-J2 spin-1/2 chain (XX, Heisenberg, BKT). Result: **NOT universal.** Hierarchy direction is model-dependent. Potts walking gives O(1) entropy deviations (26% at q=7); J1-J2 gives only O(1%) differences. The Casimir-Re(c) finding is walking-specific.

## Active Research Thread
**Walking regime: five confirmed novel findings + one scope constraint.** Casimir (098), χ_F scaling (103), entropy concentration, Rényi mapping, multiplet dominance. Sprint 104 established that energy-entropy decoupling is NOT universal — constrains Casimir finding to walking-specific mechanism.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 79 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unspent. Strongest predictions: (a) Heisenberg chain c_eff on 5-10 qubits, (b) q=2 Ising χ_F at g_c with exact α=0.98.
2. **χ_F in J1-J2 chain** — Does χ_F show BKT-specific scaling? α should be different from Potts walking. Quick test at J2c with same infrastructure.
3. **Publish preparation** — Five confirmed novel results. Write a unified summary of the walking research program with scope constraints from Sprint 104.

## What's Been Ruled Out
- Energy-entropy hierarchy universality (Sprint 104) — NOT universal, walking-specific
- Entanglement asymmetry as walking probe (Sprint 102)
- SREE as walking discriminator (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)

## Key Tools Available
- J1-J2 chain Sz-sector exact diag (Sprint 104): N≤20 in <20s
- Vectorized S_q Potts builder (Sprint 098/103)
- χ_F infrastructure: coupling/field split for fast g-scan
- GPU eigsh: q=5 n=10 (156s), q=7 n=8 (90s)
- S_q Potts DMRG: q=2 (fast, n≤24+), q=5 (n≤12, chi≤50)
- IBM QPU: 580s remaining
