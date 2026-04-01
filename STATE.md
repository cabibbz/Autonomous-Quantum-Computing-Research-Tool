# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 040 — q=4 Potts MI-CV: Marginal transition shows crossing, not dome

## Active Research Thread
MI-CV as universal phase transition classifier. q=4 Potts (marginal point where 2D Potts goes first-order) still shows crossing curves at n=8,12, same as q=3 and TFIM. Marginal character not visible at these sizes — possibly hidden by logarithmic corrections.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **q=5 Potts MI-CV** — q=5 is definitively first-order (above marginal q=4). Should show step function (no crossings). This would be the definitive test: crossings at q≤4, step at q≥5.
2. **q=4 at n=16** — Test if logarithmic corrections at q=4 start broadening the crossing region at larger sizes. Would reveal marginal character.
3. **Universal scaling functions** — Plot collapsed F(x) for q=2,3,4 on same axes. Visually distinct curves = MI-CV captures full universal function.

## What's Been Ruled Out
- Small-scale QEC active correction: exhaustively proven to fail for [[5,1,3]] (Sprints 026-028)
- Pauli fraction as BW metric: fails at large n (Sprint 035)
- NN-only MI-CV: dominated by boundary effects (Sprint 036)
- BW sin_inv envelope: wrong at finite size (Sprint 035)
- Simple power-law fit for ν at n≤32: corrections to scaling too large (Sprint 037)
- ν=0.755 from crossing points: artifact (Sprint 038)
- Fixed g_c collapse: misleading due to finite-size g_c shift (Sprint 039)
- n=24 qutrit DMRG without symmetry: too slow (>60s per point) (Sprint 039)
- q=4 marginal showing dome: shows crossing at n=8,12 (Sprint 040)

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- Correlation-function MI reconstruction: all-pairs MI from MPS at any size (d=2 Pauli, d=3 Gell-Mann, d=4 SU(4) generators)
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
- SU(d) generalized Gell-Mann basis: implemented for d=3,4. ~6.5s/point at d=4 n=8, ~18s at d=4 n=12
