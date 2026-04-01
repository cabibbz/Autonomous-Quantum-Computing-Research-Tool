# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 042 — True q=5 Potts MI-CV: Second-order crossings at g_c≈0.41, NOT first-order

## Active Research Thread
MI-CV as universal phase transition classifier. True Potts q=5 (Kronecker-delta coupling) shows crossing curves at g_c≈0.41, disproving the prediction of first-order step function. Both Clock and Potts q=5 are second-order in 1D quantum, despite 2D classical Potts being first-order for q>4. Potts slope (4.5-6.1) is 5.7x steeper than Clock (0.8-1.0) but same qualitative signature.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **Potts q=5 data collapse** — With n=8,12 data, add n=16 point (if feasible ~50s) and attempt collapse. Does optimal ν match any known universality? If ν differs from both Ising (1) and 3-state Potts (5/6), it's a new universality class.
2. **Potts q=3 vs q=5 slope comparison** — Potts q=3 has slope ~n^1.36 (Sprint 039). Potts q=5 has slope ratio 1.35 at n=8→12. Similar? Would test whether 1D quantum Potts is always Ising-like or has q-dependent exponents.
3. **First-order search: large q Potts** — q=5 is second-order. At what q does 1D quantum Potts actually become first-order? Try q=10 or q=20 at n=8. If still crossing, 1D quantum Potts may NEVER be first-order.

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
- q=5 clock showing no crossings: still shows crossings at g≈0.67 (Sprint 041)
- ClockChain as Potts proxy for q≥4: clock ≠ Potts at q≥4 (Sprint 041)
- **q=5 Potts first-order: shows second-order crossings at g≈0.41 (Sprint 042)**

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- Correlation-function MI reconstruction: all-pairs MI from MPS at any size (d=2-5 tested)
- Custom PottsChain model: Kronecker-delta coupling for any q (Sprint 042)
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
- SU(d) generalized Gell-Mann basis: implemented for d=2-5. ~17s/point at d=5 n=8
