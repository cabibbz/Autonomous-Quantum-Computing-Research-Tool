# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 079 — c(q=7,10) from exact diag + DMRG. Complex CFT Re(c) matches c_eff at q=5 (1%) but FAILS at q=7 (18% deficit) and q=10 (40% deficit). Walking regime ξ*(q) shrinks with q. DMRG = exact at q=7 n=8 (zero chi truncation). q=5 is the unique "sweet spot" for complex CFT.

## Active Research Thread
**S_q Potts: walking regime breakdown characterized.** Complex CFT c_eff ≈ Re(c) works only for q near 5. At q≥7, walking correlation length ξ* ≤ ~10 sites; system "sees" first-order nature.

c_eff at n=8 (exact diag, matched sizes):

| q | Re(c)_CFT | c_eff(n=8) | c/Re(c) | Walking? |
|---|-----------|------------|---------|----------|
| 3 | 0.800 | 0.894 | 1.118 | (real CFT, +12% FSS) |
| 5 | 1.138 | 1.141 | 1.003 | YES — sweet spot |
| 7 | 1.351 | 1.111 | 0.822 | BREAKING DOWN |
| 10 | 1.584 | 0.946* | 0.597 | NO — first-order visible |

*q=10 at n=6 only

DMRG data at larger n:
- q=5 n=24: c=1.15, ratio 1.01
- q=7 n=12: c=1.06, ratio 0.78 (worse with increasing n)

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Walking length ξ*(q=7) from DMRG** — Push q=7 to n=16-24 with chi=120-200. Find where c_eff stabilizes or gap×N starts decreasing. Need to test chi convergence first at n=12.
2. **c(q=6,8) exact diag** — Fill in the c_eff/Re(c) curve at q=6 (just above walking threshold?) and q=8 (between 7 and 10). Map the walking breakdown more precisely.
3. **Hardware validation** — 580s QPU remaining (54 sprints since last use). Best candidate: conformal tower degeneracy (q-1 fold S_q vs 2-fold Z_q) on n=6-8 qubits.

## What's Been Ruled Out
- Complex CFT c prediction at q≥7 (c_eff 18-40% below Re(c))
- DMRG chi truncation as explanation for low c at q=7 (exact diag confirms)
- q-independent FSS overshoot (overshoot is +12% for q=3, +0% for q=5, -18% for q=7)
- All items from Sprint 078 STATE.md

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤5, GPU: n≤10 for q≤5, n≤8 for q=7, n≤6 for q=10
- S_q Potts DMRG: working for q=3,5,7 (n=12 at q=7 takes 350s)
- Exact g_c = 1/q for S_q Potts
- Complex CFT formula: √Q=2cos(π/p), c=1-6/[p(p-1)], verified q=2-5
- IBM QPU: 580s remaining
