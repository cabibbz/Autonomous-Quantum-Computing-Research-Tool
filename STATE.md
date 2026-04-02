# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 087 — Entanglement spectrum scaling via DMRG. Tail weight grows as UNBOUNDED power law: q=5 ~ n^1.72, q=7 ~ n^2.52. No saturation. Entanglement gap closes as -0.43·ln(n). Energy-entropy decoupling is TEMPORARY — breaks at n≈150 (q=5) or n≈70 (q=7). Weight flows from (q-1) multiplet to tail, not from ground state.

## Active Research Thread
**S_q Potts walking regime: entanglement spectrum fully characterized across sizes.**

Tail weight scaling (DMRG midchain, open BC):

| q | power law | exponent | R² | w_tail=10% at n≈ |
|---|-----------|----------|----|-------------------|
| 5 | 1.9e-5 × n^b | 1.72 | 0.993 | 147 |
| 7 | 2.1e-6 × n^b | 2.52 | 0.989 | 72 |

Key insight: the "entropy-only" breakdown from Sprints 082-083 is a finite-size effect. At large n, ALL observables will see walking breakdown. Different observables couple to different entanglement spectrum levels, creating a hierarchy of crossover scales.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unused for 62 sprints. q=2 Ising at n=6-8, measure gap or entanglement entropy on QPU.
2. **Rényi spread from size pairs** — Recompute (c₂-c_∞) spread using size-pair extraction with DMRG data. Still a walking discriminator?
3. **Tail weight at q=3 (real CFT)** — Does tail weight grow as power law for q=3 too, or is the growth specific to walking/complex CFT?

## What's Been Ruled Out
- Tail weight saturation at accessible sizes (power law with R²>0.99, Sprint 087)
- α=3 as Re(c) recovery probe (single-size extraction artifact, Sprint 086)
- c_∞ (min-entropy) as Re(c) recovery probe (degrades like c_1)
- Open-BC CC profile fit for Rényi c_α extraction (R²<0.9 at n≤12)
- Oscillatory correlator corrections at r ≤ 7 (Im(x_σ) too small)
- Open-BC raw power law for x extraction (inflates η by 5×)
- x_σ as walking discriminator (nearly constant across q)
- Casimir energy as walking discriminator (follows Re(c) for all q)
- Entanglement gap as walking discriminator (INCREASES with q at fixed n)
- q=6 as second walking case (c_eff drops 2.9% n=8→12)
- Sharp walking-to-first-order transition (smooth crossover)

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=5 (fast, n≤24), q=6 (slow, n≤12), q=7 (n≤12)
- Periodic exact diag correlator: q=2 n≤14, q=3 n≤10, q=5 n≤8, q=7 n≤7
- Entanglement spectrum + Rényi entropies: same size limits as exact diag
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
