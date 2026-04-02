# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 080 — c_eff at q=6,8,9 from exact diag. Walking boundary mapped: c_eff/Re(c) = exp(-0.105(q-5)), ratio=1 at q≈5. q=6 marginal (0.92, stable with n). q≥7 breaking (ratio ≤0.78, c_eff decreasing with n).

## Active Research Thread
**S_q Potts: walking regime boundary fully characterized.** q=5 is the exact threshold where complex CFT Re(c) matches c_eff.

Complete c_eff/Re(c) at best available sizes:

| q | Re(c) | c_eff | n_best | c/Re(c) | Size trend |
|---|-------|-------|--------|---------|------------|
| 3 | 0.800 | 0.893 | 24 | 1.12 | stable (real CFT) |
| 5 | 1.138 | 1.147 | 20 | 1.01 | stable |
| 6 | 1.253 | 1.147 | 8 | 0.92 | STABLE (+0.1%) |
| 7 | 1.351 | 1.059 | 12 | 0.78 | DEGRADING (−5%/size) |
| 8 | 1.438 | 1.062 | 7 | 0.74 | degrading (−1%) |
| 9 | 1.515 | 1.012 | 6 | 0.67 | — |
| 10 | 1.584 | 0.946 | 6 | 0.60 | — |

Fit: c_eff/Re(c) = 1.004 × exp(−0.105(q−5))

Key insight: gap×N INCREASES with q even as walking breaks down — breakdown is in entropy, not gap.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **ξ*(q=6) via DMRG** — q=6 is marginal walking (0.92 ratio, stable with n at n=8). Push to n=12-24 with DMRG to see if c_eff stays at 1.15 or drops toward Re(c)=1.25. This would confirm/deny q=6 as a second walking case.
2. **Walking length ξ*(q=7) from DMRG** — At n=12, c_eff=1.06 and falling. Push to n=16-24 to find where it stabilizes.
3. **Hardware validation** — 580s QPU remaining (55 sprints since last use). Best candidate: conformal tower degeneracy (q-1 fold S_q vs 2-fold Z_q) on n=6-8 qubits. Or: gap×N at g_c for q=5,6 to test walking on real hardware.

## What's Been Ruled Out
- Complex CFT c prediction at q≥7 (c_eff 22-40% below Re(c))
- Sharp walking-to-first-order transition (degradation is smooth exponential)
- q-independent FSS overshoot (overshoot varies from +12% to −40%)
- c_eff as discriminator between q values at finite n (all q≥5 give c≈1.0-1.15)

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6, n≤7 for q=8
- S_q Potts DMRG: working for q=3,5,7 (needs testing for q=6,8)
- Exact g_c = 1/q for S_q Potts
- Complex CFT formula verified q=2-5
- IBM QPU: 580s remaining
