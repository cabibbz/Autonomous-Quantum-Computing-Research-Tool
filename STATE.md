# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 035 — BW locality size scaling (TFIM, using DMRG for n>10)

## Active Research Thread
Bisognano-Wichmann entanglement Hamiltonian locality — how symmetry, local dimension, and system size control whether H_E = physical Hamiltonian × entanglement temperature.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **Finite-size scaling of MI-CV** — use DMRG at n=8,16,32,50 for TFIM and XXZ. Does the dome/inflection/jump distinction sharpen? This upgrades Sprint 030's observation to a result.
2. **BW locality vs system size** — does TFIM's 9% non-BW gap shrink with n? Does the H/G-inv ratio prediction hold at n=50?
3. **Potts q-sweep** — BW locality at q=2,3,4,5 to add data points to the H/G-inv predictor curve.

## What's Been Ruled Out
- Small-scale QEC active correction: exhaustively proven to fail for [[5,1,3]] (Sprints 026-028)
- Code comparisons under symmetric noise at n≤10 (Sprints 015-024)
- Non-FT, flag-FT, and repeated syndrome extraction (Sprints 026-028)
- Bias-tailored codes at small scale (Sprints 022-023)

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
- Qiskit + Aer for circuit simulation with gate-level noise
