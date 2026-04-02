# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 100 — DMRG Casimir: Open vs Periodic BC. Open-BC extraction fails for q≥5 (boundary corrections 36-64%). Periodic-BC pairwise c/Re(c) ≈ 1.00 for all q=2-7 at exact diag sizes (N≤14). Casimir-Re(c) result confirmed but size-limited. Im(c) oscillation detection impossible at any accessible size (literature confirms).

## Active Research Thread
**Casimir/walking thread FULLY CLOSED.** Three confirmed findings:
1. Casimir energy tracks Re(c) across walking boundary (periodic BC, Sprint 098)
2. Walking breakdown is exclusively an entropy phenomenon (Sprints 082-084)
3. Open-BC extension to DMRG sizes fails due to boundary corrections (Sprint 100)

**Need a NEW research direction.**

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 75 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU expiring. Strongest prediction: q=2 Ising ground state at g_c on 5-8 qubits. Test: BW R²>0.999 at nA=3 via state tomography. Or: entanglement entropy scaling.
2. **Non-Hermitian S_q Potts** — Following Tang et al. (2024): add non-Hermitian perturbation to access complex CFT fixed points directly. Extract both Re(c) and Im(c). Would be first measurement on our hybrid model.
3. **Entanglement asymmetry / symmetry resolution** — New observable: decompose entanglement entropy by symmetry sector. Recent literature (Ares et al. 2023) on symmetry-resolved entanglement in critical chains. Could reveal new structure at walking boundary.

## What's Been Ruled Out
- Im(c) oscillation detection (any method at accessible sizes): RULED OUT (Sprints 099-100)
- Open-BC Casimir extraction for q≥5: RULED OUT (Sprint 100)
- All BW correction approaches: RULED OUT (Sprint 097)
- All previously ruled-out items still apply

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=7
- Dense periodic Casimir: q=2 N=4-14, q=5 N=4-10, q=7 N=4-8
- S_q Potts DMRG: q=2 (fast, n≤24+), q=5 (n≤12, chi≤50), q=7 (n≤8 exact only)
- IBM QPU: 580s remaining
