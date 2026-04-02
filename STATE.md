# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 073 — Cylinder Ly convergence: q=3 Ly=3 (g_c=0.797, 49.7% to 2D) and q=2 Ly=4 (g_c=0.688, 84.0% to 2D). q=3 converges 1.9x slower than q=2. Convergence is q-dependent, NOT universal.

## Active Research Thread
**Cylinder dimensional crossover — characterization mature for q=2, still early for q=3.**

Full Ly convergence dataset:

| q | Ly | g_c(cyl) | Progress to 2D | Crossing pairs |
|---|-----|----------|---------------|----------------|
| 2 | 1 | 0.250 | 0% | — |
| 2 | 2 | 0.451 | 38.6% | 3 |
| 2 | 3 | 0.655 | 77.7% | 4 |
| 2 | 4 | 0.688 | 84.0% | 2 |
| 2 | 2D | 0.771 | 100% | — |
| 3 | 1 | 0.333 | 0% | — |
| 3 | 2 | 0.565 | 24.8% | 3 |
| 3 | 3 | 0.797 | 49.7% | 1 |
| 3 | 2D | 1.267 | 100% | — |

Exponential fit g_c(Ly) = g_c(2D) - A·exp(-Ly/B):
- q=2: B = 1.60, A = 0.99. Predicts Ly=5: 0.728, Ly=6: 0.748.
- q=3: B = 3.03, A = 1.31. Predicts Ly=4: 0.917, Ly=5: 1.015.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU remaining, last used Sprint 025 (48 sprints ago). Best candidate: 1D q=2 conformal tower (n=6-8 qubits) or q=3 gap scaling.
2. **q=5 Ly=3 cylinder** — dim=125^Lx. Lx=3: 1.95M (feasible), Lx=4: 244M (infeasible). Single-size only — need alternative method (e.g., compare gap at fixed Lx to interpolate g_c).
3. **q=3 Ly=4 cylinder** — dim=81^Lx. Lx=3: 531k, Lx=4: 43M (GPU borderline). Would give 2 crossing pairs and test if q=3 convergence follows the B=3.03 prediction (expect ~62% at Ly=4).

## What's Been Ruled Out
- Floating phase in 1D hybrid model (Sprint 067)
- Weakly first-order in 1D at q=10 (Sprint 066)
- Hybrid = clock universality (Sprint 065+067)
- c·x₁ = 1/9 exact (Sprint 063c)
- L=2 in scaling regime for 2D torus (Sprints 068-070)
- DMRG for q≥5 cylinders at chi<50 (Sprint 071b)
- First-order transition on Ly=2 cylinder for q=2,3,5 (Sprints 071c, 072b)
- Universal Ly convergence rate across q (Sprint 073) — convergence is q-dependent

## Key Tools Available
- 2D torus exact diag: L×L up to ~10 sites total
- Cylinder exact diag: Ly=2 Lx≤6 for q=3, Ly=3 Lx≤4 for q=3, Ly=4 Lx≤5 for q=2
- 1D exact diag CPU: n≤8 for q≤5. GPU extends to n≤10.
- DMRG: q≤5 1D chains at n≤24. q=2 cylinder at Lx≤20. q≥5 cylinder IMPRACTICAL.
- IBM QPU: 580s remaining
