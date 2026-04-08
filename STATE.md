# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 125 — Methodology stress-test: spectral chi_F captures exactly 50% of true chi_F due to missing factor-2 prefactor in the Lehmann sum. After 2x correction, spectral matches exact to 0-2.2%. Increasing k from 12 to 200 adds zero new states — truncation is NOT the issue. Dominant fraction is 92-99%, not the claimed 1.000. Scaling exponents (alpha, z_m, beta_me) are UNAFFECTED by the prefactor.

## CRITICAL: Factor-2 Error in Spectral chi_F (Sprint 125)
All spectral chi_F VALUES in Sprints 106-124 are off by factor 2. Scaling EXPONENTS (alpha, z_m, beta_me) are unaffected (constant prefactor cancels in log-log slopes). "frac=1.000" was wrong — true dominant fraction is 0.92-0.99.

## CRITICAL: Model Identity Correction (April 2026 Audit)
**All experiments from Sprint 076 onward use the STANDARD S_q Potts model**, not the "Potts-clock hybrid." Sprints 119-122 compare both models. See KNOWLEDGE.md.

## Active Research Thread
**chi_F methodology now validated.** The spectral exponents are reliable despite the prefactor error. Key remaining question: for hybrid q=5 (dominance=92%), does the 8% non-dominant contribution affect exponent extraction at the margin?

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — BLOCKED on credentials (qiskit-ibm.json empty)

## Top 3 Next Experiments
1. **Fix spectral formula + revalidate hybrid q=5** — multiply by 2, verify exponents unchanged where dominance is lowest (92%)
2. **Compile hybrid model results for paper** — systematic data, corrected methodology, honest error bars via fss_utils
3. **q=5 DMRG chi_F extension** — walking regime open-BC alpha drift

## What's Been Ruled Out
- "frac=1.000" single-multiplet dominance — actual frac is 0.92-0.99 (Sprint 125)
- k-truncation as source of spectral error — k=200 same as k=12, only 2-6 states contribute (Sprint 125)
- iDMRG overlap chi_F for S_q Potts — non-abelian symmetry breaks convergence (Sprint 124)
- Hybrid chi_F super-scaling — alpha < 2 for all q>=5
- alpha(q=4)=2.0 for S_q at n<=20 — measured 1.52 open, 1.76 periodic
- VQE for critical GS — fidelity ~50% at n=4 (Sprint 123b)
- q_cross at integer q=4 — actually 3.58+/-0.04 (Sprint 121)

## Key Tools Available
- chi_F spectral (periodic BC): q=2 n<=18, q=3 n<=14, q=4 n<=11, q=5 n<=10, q=7 n<=8, q=10 n<=7
- chi_F DMRG (open BC): q=2 n<=24, q=4 n<=20, q=5 n<=10
- Hybrid chi_F spectral: q=3 n<=12, q=4 n<=11, q=5 n<=10, q=7 n<=8, q=10 n<=7
- **NEW**: hamiltonian_utils.py (shared builders), fss_utils.py (error bars via lmfit)
- GPU eigsh practical limit: ~10M dim
- IBM QPU: 580s remaining (credentials needed)
