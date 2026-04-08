# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 124 -- DMRG chi_F for q=4 at larger sizes. iDMRG overlap FAILED (non-abelian symmetry breaks convergence). Finite DMRG extended to n=20 (open BC): pairwise alpha drifts UP from 1.505 to 1.523, corrected alpha=1.524+/-0.002. Log-corrected alpha=2 still disfavored but drift direction is consistent with eventual convergence. q=4 log correction regime remains beyond our tools (need n>100).

## CRITICAL: Model Identity Correction (April 2026 Audit)
**All experiments from Sprint 076 onward use the STANDARD S_q Potts model**, not the "Potts-clock hybrid." Sprints 119-122 compare both models. See KNOWLEDGE.md.

## Active Research Thread
**chi_F spectral decomposition -- model comparison nearly complete.** Key results across Sprints 119-124:
- Hybrid z_m crosses 1 at q_cross=3.58 (walking -> continuous boundary)
- S_q q=4: alpha_periodic=1.777, alpha_open=1.524 (drifting upward)
- Log correction fits: alpha=2 form fails at n<=20 (open BC)
- No prior chi_F log correction measurement in literature
- iDMRG overlap chi_F doesn't work for non-abelian models
QPU thread: circuits ready, need credentials.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s -- BLOCKED on credentials (qiskit-ibm.json empty)

## Top 3 Next Experiments
1. **Compile hybrid model results** -- systematic chi_F data across q=[3,4,5,7,10], z_m(q) crossing, model comparison. Paper-worthy dataset.
2. **q=5 DMRG chi_F extension** -- walking regime open-BC alpha drift (expect different behavior than q=4 BKT).
3. **Configure QPU credentials** and submit prepared TFIM circuits.

## What's Been Ruled Out
- iDMRG overlap chi_F for S_q Potts -- non-abelian symmetry causes random convergence (Sprint 124)
- Hybrid chi_F super-scaling -- alpha < 2 for all q>=5
- alpha(q=4)=2.0 for S_q at n<=20 -- measured 1.52 open, 1.76 periodic
- Log-corrected alpha=2 fit at n<=20 (open BC) -- worst R^2 (Sprint 124)
- VQE for critical GS -- fidelity ~50% at n=4 (Sprint 123b)
- Adiabatic Trotter prep -- too few steps for critical gap (Sprint 123c)
- q_cross at integer q=4 -- actually 3.58+/-0.04 (Sprint 121)

## Key Tools Available
- chi_F spectral (periodic BC): q=2 n<=18, q=3 n<=14, q=4 n<=11, q=5 n<=10, q=7 n<=8, q=10 n<=7
- chi_F DMRG (open BC): q=2 n<=24, q=4 n<=20, q=5 n<=10
- Hybrid chi_F spectral: q=3 n<=12, q=4 n<=11, q=5 n<=10, q=7 n<=8, q=10 n<=7
- GPU eigsh practical limit: ~10M dim
- IBM QPU: 580s remaining (credentials needed)
- QPU circuits: 6 TFIM circuits ready (11 CX each)
