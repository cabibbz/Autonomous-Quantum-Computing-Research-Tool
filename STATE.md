# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 101 — Symmetry-Resolved Entanglement Entropy (SREE). First SREE measurement for q=2-10 S_q Potts. S_number/S_total ≈ 0.908 universal across q at n=6 — q-independent at fixed geometry. Charge-0 enrichment increases with q (p(0)*q from 1.5 to 5.2). No walking-specific signature in SREE at accessible sizes.

## Active Research Thread
**SREE thread: one sprint, findings documented.** Key result: S_n/S_t universality across q. Not a strong walking discriminator. Ready for new direction.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 76 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU expiring. Strongest prediction: q=2 Ising ground state at g_c on 5-8 qubits. Test: entanglement entropy scaling, or BW R²>0.999 at nA=3 via state tomography.
2. **Non-Hermitian S_q Potts** — Add non-Hermitian perturbation to access complex CFT fixed points directly. Extract both Re(c) and Im(c). Would extend Tang et al. (2024) to our model.
3. **Entanglement asymmetry** — Ares et al. (2023): measures how much symmetry is broken in subsystem. Different from SREE — probes spontaneous symmetry breaking at finite size. Could distinguish walking (approximate SSB) from real CFT (no SSB at criticality).

## What's Been Ruled Out
- SREE as walking discriminator at accessible sizes (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q≥5 (Sprint 100)
- All BW correction approaches (Sprint 097)
- All previously ruled-out items still apply

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=7
- Dense periodic Casimir: q=2 N=4-14, q=5 N=4-10, q=7 N=4-8
- S_q Potts DMRG: q=2 (fast, n≤24+), q=5 (n≤12, chi≤50), q=7 (n≤8 exact only)
- SREE projectors: works up to dimA ≈ 5000 (nA*q combinations)
- IBM QPU: 580s remaining
