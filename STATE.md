# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 086 — Rényi entropy scaling via DMRG. Sprint 085's α=3 finding was a single-size extraction artifact. Size-pair extraction (cancels non-universal c'_α) shows α=1 is consistently best for both walking and broken regimes. Optimal α shifts gradually: 0.5 (real CFT) → 1 (walking) → 2 (broken). No Rényi index recovers Re(c) for q≥7. Open-BC profile fits unreliable (R²<0.9).

## Active Research Thread
**S_q Potts walking regime: Rényi analysis corrected, entropy-only breakdown fully characterized.**

Walking breakdown is confirmed in ALL Rényi entropies, with NO α that recovers Re(c) for q≥7:

| q | pair(periodic) best | c_best/Rec | pair(DMRG) best | c_best/Rec |
|---|---------------------|------------|-----------------|------------|
| 2 | α=0.5 | 1.000 | — | — |
| 5 | α=1 | 1.004 | α=0.5 | 0.97 |
| 7 | α=2 | 0.904 | α=1 | 0.999→0.83 (degrades) |

Key correction: Sprint 085's α=3 was contaminated by non-universal constant c'_α.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU, 61 sprints since last use. q=2 S_q Potts at n=6-8 on QPU, measure gap or entanglement.
2. **Entanglement spectrum scaling at q=7 n=16-24** — DMRG. Does tail weight saturate or grow indefinitely?
3. **Rényi spread from size pairs** — Recompute (c₂-c_∞) spread using size-pair extraction. Still a walking discriminator?

## What's Been Ruled Out
- α=3 as Re(c) recovery probe (single-size extraction artifact, Sprint 086)
- c_∞ (min-entropy) as Re(c) recovery probe (degrades like c_1)
- Open-BC CC profile fit for Rényi c_α extraction (R²<0.9 at n≤12)
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
