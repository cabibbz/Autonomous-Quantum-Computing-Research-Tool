# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 089 — q=7 tail exponent resolution + level redistribution as walking discriminator. Three key results: (1) q=7 tail exponent b≈3.0 was pre-asymptotic → converges to b≈2.07 ± 0.14 with n≥8 data, confirming universal b≈2.0 for ALL q. (2) Static entropy partition (%S(lev0), %S(lev1)) is monotonic in q and the cleanest walking discriminator. (3) Multiplet dominance ratio M/[(q-1)/q] crosses 1.0 at q≈4 — the real-to-complex CFT boundary.

## Active Research Thread
**Entanglement spectrum structure as walking probe — complete picture.**

Universal tail weight scaling confirmed (w_tail ~ n^2.0 for all q=2,3,5,7):

| q | b (n≥8) | b_err | R² | M(n=16) | M/[(q-1)/q] | type |
|---|---------|-------|------|---------|-------------|------|
| 2 | 1.976 | 0.137 | 0.986 | 0.676 | 1.328 | real |
| 3 | 2.013 | 0.138 | 0.986 | 0.724 | 1.073 | real |
| 5 | 2.012 | 0.138 | 0.986 | 0.771 | 0.957 | walking |
| 7 | 2.065 | 0.141 | 0.986 | 0.797 | 0.930 | broken |

Walking discriminators (monotonic in q): %S(lev0), %S(lev1), S_lev0_inf, M.
M/[(q-1)/q] > 1 for real CFT, < 1 for complex CFT. Crosses at q≈4.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unused for 64 sprints. q=2 Ising at n=6-8: measure entanglement entropy on real hardware. Prediction: S(half-chain) at g_c matches simulator to ~10%.
2. **q=4 entanglement spectrum** — M/[(q-1)/q] crossover at q=4 is predicted but unmeasured. q=4 is the boundary case (S_q transition changes character). Test whether M/[(q-1)/q] ≈ 1.0 exactly.
3. **Entanglement spectrum at q=6,8,9** — Fill in the walking/broken boundary. Does M/[(q-1)/q] decrease smoothly? Exact diag at q=6,8 n=6-8 would be fast.

## What's Been Ruled Out
- Tail weight growth as walking-specific (universal n^2 for all q, Sprint 088-089)
- Sprint 087 exponents 1.72, 2.52 (real-space fitting artifact)
- q=7 b≈3.0 exponent (pre-asymptotic, converges to b≈2.1, Sprint 089a)
- Dynamic transfer rates as walking discriminator (not monotonic, Sprint 089b)
- α=3 as Re(c) recovery probe (single-size extraction artifact, Sprint 086)
- c_∞ (min-entropy) as Re(c) recovery probe (degrades like c_1)
- Open-BC CC profile fit for Rényi c_α extraction (R²<0.9 at n≤12)
- Oscillatory correlator corrections at r ≤ 7 (Im(x_σ) ≈ 0.02, period ~80)
- Open-BC raw power law for x extraction (inflates η by 5×)
- x_σ as walking discriminator (nearly constant across q)
- Casimir energy as walking discriminator (follows Re(c) for all q)
- Entanglement gap as walking discriminator (INCREASES with q at fixed n)
- q=6 as second walking case (c_eff drops 2.9% n=8→12)

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=2,3 (very fast, n≤24+), q=5 (fast, n≤24), q=6 (slow, n≤12), q=7 (n≤16 confirmed Sprint 089)
- Periodic exact diag correlator: q=2 n≤14, q=3 n≤10, q=5 n≤8, q=7 n≤7
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
