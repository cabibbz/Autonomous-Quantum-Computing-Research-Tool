# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 081 — DMRG c_eff(q=6) at n=10,12. q=6 walking is MARGINAL BREAKING: c_eff drops 2.9% (n=8→12), drift rate closer to q=7 than q=5. Extrapolated breakdown at n≈38. q=5 confirmed truly unique walking sweet spot.

## Active Research Thread
**S_q Potts walking regime: q=6 resolved as marginal breaking.**

Complete c_eff convergence picture at n=12 (best common size):

| q | Re(c) | c_eff(n=12) | c/Re(c) | dc/d(ln n) | Status |
|---|-------|-------------|---------|------------|--------|
| 5 | 1.138 | 1.152 | 1.013 | +0.014 | WALKING — stable |
| 6 | 1.253 | 1.115 | 0.890 | -0.048 | MARGINAL — slowly breaking |
| 7 | 1.351 | 1.059 | 0.784 | -0.091 | BROKEN |

Key insight: Walking boundary is a CROSSOVER, not a sharp transition. q=5 uniquely stable. q=6 drifts 3.4x faster than q=5. Breakdown extrapolates to n≈38 for q=6.

gap×N remains healthy for all q (2.0-2.1) even as entropy degrades.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Walking length ξ*(q=7) from DMRG** — At n=12, c_eff=1.06 and falling. Push to n=16 to quantify breakdown rate. Compare dc/d(ln n) = -0.091.
2. **Hardware validation** — 580s QPU, 56 sprints since last use. Best candidate: gap×N at g_c for q=5 (walking signature) vs q=7 (non-walking) on n=6-8 qubits. Encode q-state system in qubits.
3. **Non-Hermitian walking extension** — arXiv:2502.02001 shows complex entanglement entropy recovers Re(c) even at q>5. Could extend our walking curve with non-Hermitian perturbation.

## What's Been Ruled Out
- q=6 as second walking case (c_eff drops 2.9% n=8→12, drift closer to q=7)
- Sharp walking-to-first-order transition (smooth crossover q=5→6→7)
- c_eff stability at n=6-8 as evidence of walking (too short a range)
- DMRG for q≥6 at n>12 with chi>30 (timing infeasible)

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=5 (fast, n≤24), q=6 (slow, n≤12 chi=30), q=7 (n≤12 chi=56)
- Exact g_c = 1/q for S_q Potts
- Complex CFT formula verified q=2-5
- IBM QPU: 580s remaining
