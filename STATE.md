# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 039 — Potts Data Collapse: ν=5/6 vs ν=1

## Active Research Thread
MI-CV as universal phase transition classifier. Data collapse now distinguishes universality classes: Potts ν=5/6 gives 14% better collapse than Ising ν=1 (joint-optimized). Slope exponent (n^1.36 Potts vs n^1.1 Ising) provides second discriminator.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **q=4 Potts MI-CV** — q=4 Potts is marginal (BKT-like). Does MI-CV show dome signature instead of crossing? Would complete the MI-CV classification: crossing=2nd order, step=1st order, dome=BKT.
2. **Potts n=24-32 with Z₃ conserve** — Need Z₃ symmetry in DMRG to reach larger sizes. Would sharpen ν convergence from 0.87→5/6.
3. **Universal scaling functions** — After collapse, plot TFIM and Potts scaling functions F(x) on same axes. If visually distinct, MI-CV captures not just ν but the full universal function.

## What's Been Ruled Out
- Small-scale QEC active correction: exhaustively proven to fail for [[5,1,3]] (Sprints 026-028)
- Pauli fraction as BW metric: fails at large n (Sprint 035)
- NN-only MI-CV: dominated by boundary effects (Sprint 036)
- BW sin_inv envelope: wrong at finite size (Sprint 035)
- Simple power-law fit for ν at n≤32: corrections to scaling too large (Sprint 037)
- ν=0.755 from crossing points: artifact (Sprint 038)
- Fixed g_c collapse: misleading due to finite-size g_c shift (Sprint 039)
- n=24 qutrit DMRG without symmetry: too slow (>60s per point) (Sprint 039)

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- Correlation-function MI reconstruction: all-pairs MI from MPS at any size (d=2 Pauli, d=3 Gell-Mann)
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
