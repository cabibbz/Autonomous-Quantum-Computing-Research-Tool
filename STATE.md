# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 038 — Data Collapse & Potts MI-CV

## Active Research Thread
MI-CV as universal phase transition classifier. Data collapse confirms Ising universality (optimal ν→1 with increasing n). Potts q=3 shows same crossing signature, extending MI-CV to d>2 via Gell-Mann basis.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **Potts n=16 + data collapse** — Get n=16 Potts MI-CV (needs ~40s), extract ν from crossing points. Compare Potts collapse with ν=5/6 vs ν=1 to distinguish universality classes via MI-CV.
2. **2D Potts q=5 MI-CV** — q=5 Potts is first-order in 2D. Does MI-CV show step function (no crossing) for d>2? Would validate that crossing/no-crossing diagnostic is universal.
3. **Entanglement Hamiltonian at first-order** — Is BW valid at the FM→XY transition (Δ=-1 in XXZ)? First-order may break BW locality. Connect BW to MI-CV.

## What's Been Ruled Out
- Small-scale QEC active correction: exhaustively proven to fail for [[5,1,3]] (Sprints 026-028)
- Pauli fraction as BW metric: fails at large n (Sprint 035)
- NN-only MI-CV: dominated by boundary effects (Sprint 036)
- BW sin_inv envelope: wrong at finite size (Sprint 035)
- Simple power-law fit for ν at n≤32: corrections to scaling too large (Sprint 037)
- ν=0.755 from crossing points: artifact of corrections to scaling; data collapse gives ν≈1 (Sprint 038)

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- Correlation-function MI reconstruction: all-pairs MI from MPS at any size (d=2 Pauli, d=3 Gell-Mann)
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
