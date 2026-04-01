# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 037 — Critical exponent extraction & first-order transition test

## Active Research Thread
MI-CV as universal phase transition classifier. Three transition types have distinct MI-CV signatures: crossing (Ising), step function (first-order), dome (BKT). Critical exponent ν from crossing points converges toward Ising (ν=1) but slowly at n=8-32.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **Data collapse** — Plot CV vs (h-h_c)·n^{1/ν} with ν=1 (Ising). If curves collapse, MI-CV is in Ising universality class despite effective ν<1 at small n.
2. **Potts MI-CV** — q=3 Potts in 1D (effectively second-order). Does MI-CV behave like Ising or show new signature? Needs TeNPy Potts model.
3. **Entanglement Hamiltonian at first-order** — Is BW valid at the FM→XY transition? First-order may break BW locality.

## What's Been Ruled Out
- Small-scale QEC active correction: exhaustively proven to fail for [[5,1,3]] (Sprints 026-028)
- Pauli fraction as BW metric: fails at large n (Sprint 035)
- NN-only MI-CV: dominated by boundary effects (Sprint 036)
- BW sin_inv envelope: wrong at finite size (Sprint 035)
- Simple power-law fit for ν at n≤32: corrections to scaling too large (Sprint 037)

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- Correlation-function MI reconstruction: all-pairs MI from MPS at any size
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
