# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 117 -- q=30 chi_F: pair(3,4) alpha=5.38. dAIC gap between log+loglog and pure log narrowed to 0.8 (from 1.4). LOO cross-validation nearly identical. **Pure logarithmic alpha(q) = 1.87 ln(q) - 0.97 is the preferred model** (simpler, stable). 11 data points q=5-30. Diminishing returns on extending further.

## Active Research Thread
**Walking regime: six confirmed novel findings, alpha(q) extended to q=30.** Casimir-Re(c) (098), chi_F scaling alpha(q) (103/110-117), chi_F spectral mechanism (107), entropy concentration, Renyi mapping, multiplet dominance. Sprint 117: pure logarithmic confirmed as preferred model.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s -- UNUSED FOR 92 SPRINTS

## Top 3 Next Experiments
1. **Hardware validation** -- 580s QPU unspent for 92 sprints. q=2 Ising chi_F peak on ~10 qubits. Strongest prediction: chi_F ~ N^1.0 at g_c with known nu=1. Design: prepare ground state via VQE or adiabatic, measure fidelity at g_c +/- delta.
2. **q=4 BKT convergence via iDMRG** -- Open-BC DMRG ruled out (Sprint 113). Try iDMRG for bulk chi_F density, bypassing boundary effects.
3. **Cross-model chi_F test** -- Test alpha(q) in pure Z_q clock model. If clock also shows logarithmic alpha(q), it's a universal feature of BKT-adjacent transitions; if not, it's specific to our hybrid.

## What's Been Ruled Out
- alpha(q) log+loglog subleading correction -- dAIC narrowed to 0.8 with 11 pts (Sprint 117). Not significant.
- alpha(q) quadratic -- unphysical (Sprint 116)
- alpha(q) linear, sqrt, power-law -- dAIC>=11.3 (Sprint 117)
- DMRG open-BC chi_F for alpha(q) extension -- boundary effects too large (Sprint 113)
- Sprint 108's 1/ln(N) correction form for q=2,3 (retracted Sprint 109)
- Energy-entropy hierarchy universality (Sprint 104) -- walking-specific
- chi_F super-scaling universality (Sprint 105) -- walking-specific
- SREE as walking discriminator (Sprint 101)
- Im(c) oscillation detection (Sprints 099-100)
- Open-BC Casimir extraction for q>=5 (Sprint 100)
- All BW correction approaches (Sprint 097)

## Key Tools Available
- chi_F spectral (periodic BC): q=2 n<=18, q=3 n<=14, q=4 n<=11, q=5 n<=10, q=6 n<=9, q=7 n<=8, q=8-9 n<=7, q=10 n<=7, q=12 n<=6, q=15 n<=5, q=20 n<=5, q=25 n<=5, q=30 n<=4
- GPU eigsh practical limit: ~10M dim (q=25 n=5 at 1616s)
- IBM QPU: 580s remaining
