# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 085 — Rényi entropies across walking boundary. Computed c_α from size pairs for α=0.5,1,2,3,5,10,∞ across q=2-8. α=1 closest to Re(c) for walking q≤6 (q=5: 0.06% accuracy). α=3 uniquely recovers Re(c) in walking-broken regime (q=8: 0.9%). α=2 always overshoots by 5-13% but most stable. Rényi spread (c₂-c_∞) is new monotonic walking discriminator. Original hypothesis (c_∞ recovers Re(c)) WRONG.

## Active Research Thread
**S_q Potts walking regime: Rényi entropy decomposition complete.**

Walking breakdown is visible in ALL Rényi entropies, but with different signatures:

| q | c_1/Rec | c_2/Rec | c_3/Rec | c_∞/Rec | best α |
|---|---------|---------|---------|---------|--------|
| 2 | 0.996 | 1.055 | 1.114 | 1.101 | 0.5 |
| 5 | 0.999 | 1.126 | 1.132 | 1.041 | 1 |
| 6 | 0.995 | 1.116 | 1.106 | 1.010 | 1 |
| 7 | 0.938 | 1.050 | 1.027 | 0.932 | 3 |
| 8 | 0.922 | 1.024 | 0.991 | 0.896 | 3 |

The "optimal α" shifts from 0.5-1 (walking) to 3 (broken). Rényi spread (c₂-c_∞) monotonically increases with q.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU, 60 sprints since last use. Encode q=2 S_q Potts at n=6-8, measure entanglement spectrum or gap×N.
2. **Entanglement spectrum scaling with n** — Does level 1 concentration GROW with n at q=7? DMRG for n=10-16.
3. **Rényi entropy at larger n via DMRG** — Does c_3 remain closest to Re(c) at q=7-8 for n>8? Test at n=10,12 with DMRG.

## What's Been Ruled Out
- c_∞ (min-entropy) as Re(c) recovery probe at large q — c_∞ degrades like c_1
- Oscillatory correlator corrections at r ≤ 7 (Im(x_σ) too small)
- Open-BC raw power law for x extraction (inflates η by 5×)
- x_σ as walking discriminator (nearly constant across q)
- Casimir energy as walking discriminator (follows Re(c) for all q)
- Entanglement gap as walking discriminator (INCREASES with q)
- q=6 as second walking case (c_eff drops 2.9% n=8→12)
- Sharp walking-to-first-order transition (smooth crossover)

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=5 (fast, n≤24), q=6 (slow, n≤12), q=7 (n≤12)
- Periodic exact diag correlator: q=2 n≤14, q=3 n≤10, q=5 n≤8, q=7 n≤7
- Entanglement spectrum + Rényi entropies: same size limits as exact diag
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
