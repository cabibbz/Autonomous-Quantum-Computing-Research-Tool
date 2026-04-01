# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 041 — q=5 Clock MI-CV: Crossings persist, shifted dramatically to g≈0.67

## Active Research Thread
MI-CV as universal phase transition classifier. q=5 clock model still shows crossing curves at n=8,12, disproving the prediction that crossings would vanish above q=4. The crossing point shifts dramatically: 0.93 (q=2) → 0.923 (q=3) → 0.893 (q=4) → 0.673 (q=5). Slope halves from q=4→5. Key insight: ClockChain ≠ Potts for q≥4.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **True q=5 Potts model** — Build custom Potts model with Kronecker-delta coupling (H = -J Σ δ(s_i,s_j) - g Σ s_i). This IS first-order and should show step function. Would definitively separate clock vs Potts universality.
2. **q=5 clock n=16** — Check if crossing sharpens or broadens. If slope stays flat (doesn't grow with n), confirms BKT character. If slope grows, it's second-order.
3. **q=6,8 clock** — Map out how crossing point and slope evolve for large q. Does g_c → 0? Does slope → 0? Would establish the clock model's approach to the XY limit.

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

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- Correlation-function MI reconstruction: all-pairs MI from MPS at any size (d=2-5 tested)
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
- SU(d) generalized Gell-Mann basis: implemented for d=2-5. ~16.5s/point at d=5 n=8
