# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 123 — QPU feasibility for TFIM critical state. VQE and adiabatic prep both fail; Qiskit `initialize` decomposition gives only 11 CX gates for n=4 (88% signal at 1% noise). Noisy sim resolves phase transition. QPU submission BLOCKED — IBM credentials empty.

## CRITICAL: Model Identity Correction (April 2026 Audit)
**All experiments from Sprint 076 onward use the STANDARD S_q Potts model**, not the "Potts-clock hybrid." Sprints 119-122 compare both models. See KNOWLEDGE.md.

## Active Research Thread
**chi_F spectral decomposition — model comparison nearly complete.** Key results across Sprints 119-122:
- Hybrid z_m crosses 1 at q_cross=3.58 (walking -> continuous boundary)
- S_q q=4: alpha=1.777, z_m=1.092 (firmly walking)
- Log correction fits: alpha=2 form fails at n<=11
- No prior chi_F log correction measurement in literature
QPU thread: circuits ready, need credentials.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — BLOCKED on credentials (qiskit-ibm.json empty)

## Top 3 Next Experiments
1. **Compile hybrid model results** — systematic chi_F data across q=[3,4,5,7,10], z_m(q) crossing, model comparison. This is paper-worthy: unstudied model with quantitative predictions.
2. **DMRG chi_F for S_q q=4** — reach n=50+ to probe log correction regime. Needs ground state overlap method.
3. **Configure QPU credentials** and submit prepared TFIM circuits.

## What's Been Ruled Out
- Hybrid chi_F super-scaling — alpha < 2 for all q>=5
- alpha(q=4)=2.0 for S_q at n<=11 — measured 1.76
- Log-corrected alpha=2 fit at n<=11 — worst AIC model (Sprint 122)
- VQE for critical GS — fidelity ~50% at n=4 (Sprint 123b)
- Adiabatic Trotter prep — too few steps for critical gap (Sprint 123c)
- q_cross at integer q=4 — actually 3.58+/-0.04 (Sprint 121)

## Key Tools Available
- chi_F spectral (periodic BC): q=2 n<=18, q=3 n<=14, q=4 n<=11, q=5 n<=10, q=7 n<=8, q=10 n<=7
- Hybrid chi_F spectral: q=3 n<=12, q=4 n<=11, q=5 n<=10, q=7 n<=8, q=10 n<=7
- GPU eigsh practical limit: ~10M dim
- IBM QPU: 580s remaining (credentials needed)
- QPU circuits: 6 TFIM circuits ready (11 CX each)
