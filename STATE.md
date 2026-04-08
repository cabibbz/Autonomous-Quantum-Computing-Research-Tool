# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 122 — Log correction fits at q=4. Power+1/N^2 correction best fit for both models at n=4-11. S_q: alpha=1.757, hybrid: alpha=1.460. Log-corrected alpha=2 form is the WORST fit (dAIC=42-70). Asymptotic log correction regime beyond exact diag reach. No prior chi_F log correction measurement in literature.

## CRITICAL: Model Identity Correction (April 2026 Audit)
**All experiments from Sprint 076 onward use the STANDARD S_q Potts model**, not the "Potts-clock hybrid." Sprints 119-122 compare both models. See KNOWLEDGE.md for full details.

## Active Research Thread
**chi_F spectral decomposition as universal diagnostic.** Model comparison nearly complete. Key results: q_cross=3.58 for hybrid z_m crossing, S_q alpha=1.76 at q=4 with no sign of approaching 2.0 at accessible sizes. Literature search: no prior chi_F log correction for q=4, no prior study of hybrid model z_m crossing.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 97 SPRINTS

## Top 3 Next Experiments
1. **QPU hardware test** — q=2 Ising chi_F ~ N^1.0 at g_c=0.5 (exact alpha=1.0). 5-8 qubits. Budget: 580s. Test strongest simulator prediction on real hardware.
2. **DMRG chi_F for S_q q=4** — reach n=50+ to probe log correction regime. Need method to compute chi_F from DMRG (ground state overlap or spectral decomposition via excited states).
3. **Compile hybrid model paper dataset** — systematic chi_F data across q=[3,4,5,7,10], z_m(q) crossing at q_cross=3.58, model comparison. Unstudied model in the literature.

## What's Been Ruled Out
- Hybrid chi_F super-scaling — hybrid alpha < 2 for all q>=5 (Sprint 119)
- alpha(q=4) = 2.0 for S_q at n<=11 — measured 1.76, no sign of approaching 2.0 (Sprints 118, 121, 122)
- Log-corrected alpha=2 fit at n<=11 — worst AIC model for both S_q and hybrid (Sprint 122)
- alpha(q) polynomial forms for S_q — log+loglog marginally best (Sprint 116-117)
- DMRG open-BC chi_F — boundary effects (Sprint 113)
- "Our model is novel" for sprints 076+ — standard S_q Potts (April 2026 audit)
- q_cross at integer q=4 — actually q_cross=3.58+/-0.04 (Sprint 121)

## Key Tools Available
- chi_F spectral (periodic BC): q=2 n<=18, q=3 n<=14, q=4 n<=11, q=5 n<=10, q=7 n<=8, q=10 n<=7
- Hybrid chi_F spectral: q=3 n<=12, q=4 n<=11, q=5 n<=10, q=7 n<=8, q=10 n<=7
- GPU eigsh practical limit: ~10M dim
- IBM QPU: 580s remaining
