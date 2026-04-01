# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 036 — MI-CV finite-size scaling (n=8→50 for TFIM, n=8→32 for XXZ)

## Active Research Thread
MI-CV as a genuine phase transition order parameter. Confirmed divergence at criticality, convergence in ordered phases. Correlation-function technique enables all-pairs MI at any DMRG-accessible size.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **Critical exponent extraction** — fit h_c(n) = h_c(∞) + a·n^{-1/ν} from MI-CV crossing points at n=8,16,32,50. Does ν match Ising (ν=1)?
2. **First-order transition test** — XXZ at Δ=-1 (FM transition). MI-CV should show discontinuous jump, not smooth crossover. Would complete the classification.
3. **Potts MI-CV** — q=3 Potts has first-order transition. Test if MI-CV jump is sharper than BKT dome. Need q>2 TeNPy model.

## What's Been Ruled Out
- Small-scale QEC active correction: exhaustively proven to fail for [[5,1,3]] (Sprints 026-028)
- Pauli fraction as BW metric: fails at large n (Sprint 035)
- NN-only MI-CV: dominated by boundary effects, doesn't capture phase transition (Sprint 036)
- BW sin_inv envelope: wrong at finite size (Sprint 035)

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- Correlation-function MI reconstruction: all-pairs MI from MPS at any size (Sprint 036 technique)
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
