# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 074 — Cylinder convergence limits: q=3 Ly=4 infeasible (838s/pt), q=5 Ly=3 g_c=0.974 (46.5% to 2D). Convergence rate saturates for q≥3 (q=3: 49.7%, q=5: 46.5% at Ly=3). Exact diag cylinder study near its limits.

## Active Research Thread
**Cylinder dimensional crossover — mature for q=2, characterization complete for q=3,5 at accessible sizes.**

Full Ly convergence dataset:

| q | Ly | g_c(cyl) | Progress to 2D | Crossing pairs | Notes |
|---|-----|----------|---------------|----------------|-------|
| 2 | 1 | 0.250 | 0% | — | |
| 2 | 2 | 0.451 | 38.6% | 3 | |
| 2 | 3 | 0.655 | 77.7% | 4 | |
| 2 | 4 | 0.688 | 84.0% | 2 | |
| 2 | 2D | 0.771 | 100% | — | |
| 3 | 1 | 0.333 | 0% | — | |
| 3 | 2 | 0.565 | 24.8% | 3 | |
| 3 | 3 | 0.797 | 49.7% | 1 | |
| 3 | 2D | 1.267 | 100% | — | Ly=4 infeasible |
| 5 | 1 | 0.441 | 0% | — | |
| 5 | 2 | 0.714 | 23.8% | 1 | |
| 5 | 3 | 0.974 | 46.5% | 1 | (Lx=2,3) crossing |
| 5 | 2D | 1.588 | 100% | — | |

Exponential fit g_c(Ly) = g_c(2D) - A·exp(-Ly/B):
- q=2: B = 1.19. Predicts Ly=5: 0.745, Ly=6: 0.760.
- q=3: B = 2.49. Predicts Ly=4: 0.917, Ly=5: 1.057.
- q=5: B = 2.83. Predicts Ly=4: 1.144, Ly=5: 1.285.

**Key: Convergence saturates for q≥3.** B(q=2)=1.19, B(q=3)=2.49, B(q=5)=2.83. The big jump is q=2→q=3.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU remaining, last used Sprint 025 (49 sprints ago). Best candidate: 1D conformal tower (q=2 or q=3) on n=6-8 qubits. Overdue.
2. **q=2 Ly=5 cylinder** — dim=32^Lx. Lx=3: 32768 (trivial), Lx=4: 1M, Lx=5: 33M. Would test whether q=2 saturation continues (predict g_c≈0.745, 95% to 2D).
3. **New research direction** — The cylinder convergence thread is near exhaustion with exact diag. Consider: (a) entanglement entropy on cylinders, (b) order parameter Binder cumulant for first-order test, (c) return to CFT operator content at larger q, (d) universality class comparison with literature predictions.

## What's Been Ruled Out
- Floating phase in 1D hybrid model (Sprint 067)
- Weakly first-order in 1D at q=10 (Sprint 066)
- Hybrid = clock universality (Sprint 065+067)
- c·x₁ = 1/9 exact (Sprint 063c)
- L=2 in scaling regime for 2D torus (Sprints 068-070)
- DMRG for q≥5 cylinders at chi<50 (Sprint 071b)
- First-order transition on Ly=2 cylinder for q=2,3,5 (Sprints 071c, 072b)
- Universal Ly convergence rate across q (Sprint 073) — convergence is q-dependent
- q=3 Ly=4 cylinder with exact diag (Sprint 074a) — 838s/pt at Lx=4
- q=5 Ly≥4 cylinder with exact diag (Sprint 074b) — Lx=4 dim=244M infeasible

## Key Tools Available
- 2D torus exact diag: L×L up to ~10 sites total
- Cylinder exact diag: Ly=2 Lx≤6 for q=3, Ly=3 Lx≤4 for q=3, Ly=4 Lx≤5 for q=2
- 1D exact diag CPU: n≤8 for q≤5. GPU extends to n≤10.
- DMRG: q≤5 1D chains at n≤24. q=2 cylinder at Lx≤20. q≥5 cylinder IMPRACTICAL.
- IBM QPU: 580s remaining
