# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 107 — Harden χ_F spectral decomposition mechanism. Result: **CONFIRMED NOVEL.** α = β_me + 2z_m - 1 exact to machine precision. q=5 now has 5 sizes (n=6-10). Cross-check via finite-diff confirms multiplet dominance. Energy and entanglement multiplet gaps are independent (different z exponents).

## Active Research Thread
**Walking regime: six confirmed novel findings.** Casimir-Re(c) (098), χ_F scaling α(q) (103), χ_F spectral mechanism (107), entropy concentration, Rényi mapping, multiplet dominance. All scoped as walking-specific (104-105).

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 82 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unspent. q=2 Ising χ_F peak at g_c is measurable on ~10 qubits. Could also test fidelity susceptibility scaling directly.
2. **Walking at q=4 BKT point** — Sprint 106c showed q=4 as intermediate (α=1.69). Detailed spectral decomp at q=4 with more sizes could reveal the BKT→walking crossover mechanism.
3. **Entanglement Hamiltonian operator content at walking** — Sprint 097 showed H_E becomes non-compact at threshold. Does the (q-1)-fold energy multiplet show up in H_E operator content?

## What's Been Ruled Out
- Energy-entropy hierarchy universality (Sprint 104) — walking-specific
- χ_F super-scaling universality (Sprint 105) — walking-specific, BKT invisible
- SREE as walking discriminator (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)
- Entanglement gap as proxy for energy multiplet gap (Sprint 107) — different z exponents

## Key Tools Available
- χ_F spectral decomposition infrastructure (Sprint 106-107): q≤7, n≤10 (GPU)
- Combined scaling: z_m(q) = 0.065q + 0.845, β_me(q) = 0.188q − 0.238, α(q) = 0.318q + 0.452
- J1-J2 chain Sz-sector exact diag (Sprint 104-105): N≤20 in <100s
- Vectorized S_q Potts builder (Sprint 098/103/107)
- GPU eigsh: q=5 n=10 (279s), q=7 n=8 (185s)
- S_q Potts DMRG: q=2 (fast, n≤24+), q=5 (n≤12, chi≤50)
- IBM QPU: 580s remaining
