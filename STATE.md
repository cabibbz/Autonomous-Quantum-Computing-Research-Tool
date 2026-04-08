# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 118 -- q=4 chi_F extended to 8 sizes (n=4-11). **alpha(q=4) = 1.771 +/- 0.001 (converged).** NOT consistent with exact q=4 Potts alpha=2.0 (11.5% off). Confirms hybrid model q=4 is NOT in standard Potts universality class. Logarithmic formula from q>=5 underpredicts (1.62 vs 1.77).

## Active Research Thread
**Walking regime findings stabilized.** alpha(q) logarithmic for q>=5 with 11 pts q=5-30. q=4 anomalous exponent 1.77 (distinct from both exact Potts 2.0 and log formula 1.62). Six confirmed novel findings.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s -- UNUSED FOR 93 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** -- 580s QPU unspent. q=2 Ising chi_F or Casimir on ~10 qubits.
2. **alpha(q) unified formula q=2-30** -- Can we find a single formula covering q=2 (exact 1.0), q=3 (exact 1.4), q=4 (1.77), q=5-30 (logarithmic)?
3. **Cross-model validation** -- Test chi_F in Z_q clock model at its BKT critical point. BKT should give different alpha behavior.

## What's Been Ruled Out
- alpha(q=4) = 2.0 (exact Potts) -- measured 1.77, stable 8 sizes (Sprint 118)
- alpha(q) log+loglog subleading -- dAIC narrowed to 0.8 (Sprint 117)
- alpha(q) quadratic, linear, sqrt, power-law -- dAIC>=7.7 (Sprint 117)
- DMRG open-BC chi_F -- boundary effects (Sprint 113)
- Plus all prior ruled-out items

## Key Tools Available
- chi_F spectral (periodic BC): q=2 n<=18, q=3 n<=14, q=4 n<=11, q=5 n<=10, q=6 n<=9, q=7 n<=8, q=8-9 n<=7, q=10 n<=7, q=12 n<=6, q=15 n<=5, q=20 n<=5, q=25 n<=5, q=30 n<=4
- GPU eigsh practical limit: ~10M dim
- IBM QPU: 580s remaining
