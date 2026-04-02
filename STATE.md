# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 102 — Fidelity Susceptibility Across Walking Boundary. χ_F validates ν for q=2 (1.009) and q=3 (0.841). Walking q=5: scaling exponent α=2.09 (first-order-like, expected 1.41 from gap ν). χ_F grows ~2.88× per unit q at fixed n. Fourth observable showing walking-specific behavior.

## Active Research Thread
**Fidelity susceptibility: strong result, needs hardening.** α≈2 at q=5 from only 2 sizes — need more data points. Also need q=6,7 with multiple sizes to map the crossover.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 77 SPRINTS

## Top 3 Next Experiments
1. **Harden χ_F scaling at q=5** — Get n=10 on GPU (~5min single calc) or use DMRG fidelity for larger n. Also measure q=6,7 at n=6,8 to map α(q) across walking boundary.
2. **Universality test in different model** — Does Casimir-entropy hierarchy appear in J1-J2 or SU(N) chains? Could upgrade finding from PRB to PRL.
3. **Hardware validation** — 580s QPU expiring. Strongest prediction: Ising ground state fidelity at g_c on 5-8 qubits.

## What's Been Ruled Out
- Entanglement asymmetry as walking probe (ΔS_A=0 for symmetric ground states) — Sprint 102
- SREE as walking discriminator at accessible sizes (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=7
- χ_F scan infrastructure: build_sq_potts_parts (coupling/field split for fast g-scan)
- S_q Potts DMRG: q=2 (fast, n≤24+), q=5 (n≤12, chi≤50)
- IBM QPU: 580s remaining
