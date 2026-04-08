# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 116 -- q=25 chi_F: alpha=5.17, **log+loglog now marginally AIC-best** (dAIC=1.4 over pure log). alpha(q) = 2.62 ln(q) - 1.77 ln(ln(q)) - 1.26. Pure log underpredicts by 4.7% at q=25. z_m=2.29. Single-multiplet dominance confirmed through q=25. n=5 at 9.8M dim took 1616s (GPU limit).

## Active Research Thread
**Walking regime: six confirmed novel findings, alpha(q) extended to q=25.** Casimir-Re(c) (098), chi_F scaling alpha(q) (103/110-116), chi_F spectral mechanism (107), entropy concentration, Renyi mapping, multiplet dominance. Sprint 116: log+loglog marginally preferred over pure log.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s -- UNUSED FOR 91 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** -- 580s QPU unspent for 91 sprints. q=2 Ising chi_F peak or Casimir on ~10 qubits. Strongest prediction: chi_F ~ N^1.0 at g_c with known nu=1.
2. **alpha(q) at q=30** -- log+loglog predicts 5.47, pure log predicts 5.37. Dim n=4: 810k (GPU easy), n=5: 24.3M (exceeds GPU). Would test subleading correction.
3. **q=4 BKT convergence via iDMRG** -- Open-BC DMRG ruled out (Sprint 113). Try iDMRG for bulk chi_F density, bypassing boundary effects.

## What's Been Ruled Out
- alpha(q) quadratic -- unphysical peak at q=32.6 (Sprint 116). DAIC=+7.7
- alpha(q) linear (0.260q+0.815 from Sprint 110) -- ruled out by q=12 data (Sprint 112)
- alpha(q) sqrt and power-law forms -- dAIC>=4.2 with 10 pts (Sprint 116)
- DMRG open-BC chi_F for alpha(q) extension -- boundary effects too large (Sprint 113)
- Sprint 108's 1/ln(N) correction form for q=2,3 (retracted Sprint 109)
- Energy-entropy hierarchy universality (Sprint 104) -- walking-specific
- chi_F super-scaling universality (Sprint 105) -- walking-specific
- SREE as walking discriminator (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q>=5 (Sprint 100)
- All BW correction approaches (Sprint 097)

## Key Tools Available
- chi_F spectral (periodic BC): q=2 n<=18, q=3 n<=14, q=4 n<=11, q=5 n<=10, q=6 n<=9, q=7 n<=8, q=8 n<=7, q=9 n<=7, q=10 n<=7, q=12 n<=6, q=15 n<=5, q=20 n<=5, q=25 n<=5
- DMRG overlap chi_F (open BC): validated but useless for alpha(q) due to BC effects
- Vectorized S_q Potts builder (fast: q=6 n=9 build in 30s vs 184s with loop)
- GPU eigsh limits: q=6 n=9 (10M, 305s), q=10 n=7 (10M, 700s), q=12 n=6 (3M, 160s), q=15 n=5 (759k, 56s), q=20 n=5 (3.2M, 333s), q=25 n=5 (9.8M, 1616s)
- IBM QPU: 580s remaining
