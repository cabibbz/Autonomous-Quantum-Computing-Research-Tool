# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 098 — Harden Casimir Finding: GPU-Extended Sizes and Convergence Test. Two experiments: (1) GPU-extended Casimir energy at q=3 n=12, q=5 n=10 (9.8M), q=7 n=8 (5.7M), q=8 n=7. Combined with Sprint 083 for 4-5 point fits per q. (2) Pairwise convergence analysis: c_implied/Re(c) from consecutive pairs converges to 1.00 for ALL q. Casimir 16× more consistent with Re(c) than entropy. Upgraded from POTENTIALLY NOVEL to CONFIRMED NOVEL.

## Active Research Thread
**Casimir/walking thread COMPLETE.** Two confirmed novel findings:
1. Casimir energy obeys Re(c) across walking boundary (16× more consistent than entropy)
2. Walking breakdown is exclusively an entropy phenomenon at finite size

**BW thread (Sprints 091-097) also COMPLETE.** BW is optimal compact H_E approximation; breakdown is fundamental.

Need a NEW research direction.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU budget unused for 73 sprints. Strongest predictions: (a) S(nA) at g_c for q=2 on QPU, (b) Casimir E₀ on QPU, (c) BW R²>0.999 at nA=3 on hardware.
2. **New direction: entanglement Hamiltonian in 2D** — BW story complete in 1D. Test on Ly=2 cylinder (DMRG accessible). Does BW threshold shift? Is non-Potts content different?
3. **Complex CFT: detect Im(c) oscillations** — Im(c) should produce oscillations in ln(N) with period ~80. High-precision multi-size measurement.

## What's Been Ruled Out
- All BW correction approaches: RULED OUT (Sprint 097)
- Entanglement spectrum as BW predictor: smooth across threshold (Sprint 096)
- Casimir hardening: DONE, confirmed (Sprint 098)
- All previously ruled-out items still apply

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=7
- S_q Potts DMRG: q=2,3 (fast, n≤24+), q=5 (n≤24, chi≤60)
- Periodic exact diag: q=2 n≤18, q=3 n≤12, q=5 n≤10
- Vectorized Hamiltonian builder: ~30s build + 130s eigsh at 10M dim
- IBM QPU: 580s remaining
