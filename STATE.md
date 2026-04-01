# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 064 — Entanglement Hamiltonian at q>4 Potts: BW locality monotonically degrades with q

## Active Research Thread
**1D quantum Potts CFT characterization complete through q=30.** Key data table:

| q | g_c | c | ν | x₁ | c·x₁ | C_sse | BW peak |
|---|-----|---|---|-----|------|-------|---------|
| 2 | 0.250 | 0.500 | 1.00 | 0.125 | 0.062 | 0.50 | 91% |
| 3 | 0.333 | 0.800 | 0.86 | 0.133 | 0.107 | 0.54 | 76.5% |
| 4 | 0.392 | ~1.00 | 0.82 | 0.117 | 0.112 | 0.46 | 56% |
| 5 | 0.441 | ~1.10 | 0.85 | ~0.10 | 0.112 | 0.38 | 42% |
| 10 | 0.684 | ~1.40 | 1.12 | ~0.083 | 0.116 | ~0.21 | 39%* |

*n_A=3, not comparable to n_A=4 values above

**Sprint 064 findings:** BW locality decreases monotonically with q at fixed n_A. H/G-inv ratio predicts ordering. G-inv fraction = exactly 1/q. Unruh-like gradient persists for all q but captures less of H_E.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **QPU validation** — 580s unspent. Strongest prediction: q=3 Potts energy gap ratio at g_c=1/3, or Ising x₁=1/8 on hardware. Worth burning budget.
2. **c·x₁ mechanism** — WHY is c·x₁ ≈ constant for Potts? Test degeneracy-weighted sum rule across (q-1) primaries.
3. **BW non-local corrections** — What ARE the 44-58% non-BW contributions at q=4,5? Decompose into k-body operators. Do they have CFT structure (descendant towers)?

## What's Been Ruled Out
- c·x₁ = 1/9 exact (Sprint 063c: q=3 exact is 8/75)
- c·x₁ universal to Z_q models (Sprint 063b: clock ≠ Potts)
- Single compact boson for q≥5 (Sprint 062)
- Luttinger liquid for q≥4 Potts (Sprint 062)
- BW locality improving with q (Sprint 064: monotonic decrease)
- sin_inv envelope for q≥4 (Sprint 064: linear wins)

## Key Tools Available
- Exact diag: n=4 for q≤30, n=6 for q≤10, n=8 for q≤5
- DMRG: q≤15 at n≤8 (chi=30), q≤5 at n≤24
- IBM QPU: 580s remaining
