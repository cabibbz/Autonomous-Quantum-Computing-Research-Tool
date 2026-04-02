# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 075 — q=2 Ly=5 cylinder g_c=0.701 (86.6% to 2D). Convergence is 1/Ly² (power-law), NOT exponential. Cylinder entropy c_eff grows with Ly (FSS overshoot), S per bond cut decreasing (area law emerging).

## Active Research Thread
**Cylinder dimensional crossover — q=2 complete (Ly=1-5), q=3,5 characterized.**

Full Ly convergence dataset (q=2):

| Ly | g_c(cyl) | Progress to 2D | Crossing pairs | Increment |
|----|----------|---------------|----------------|-----------|
| 1 | 0.250 | 0% | — | — |
| 2 | 0.451 | 38.6% | 3 | +0.201 |
| 3 | 0.655 | 77.7% | 4 | +0.204 |
| 4 | 0.688 | 84.0% | 2 | +0.033 |
| 5 | 0.701 | 86.6% | 1 | +0.013 |
| 2D | 0.771 | 100% | — | — |

Best fit: g_c(Ly) = 0.771 - 1.285/Ly^2.03 (power-law, RMS=0.016). Predicts Ly=8: 0.752, Ly=10: 0.759.

Other q values (from prior sprints):

| q | Ly=2 | Ly=3 | 2D |
|---|------|------|----|
| 3 | 0.565 | 0.797 | 1.267 |
| 5 | 0.714 | 0.974 | 1.588 |

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU remaining, last used Sprint 025 (50 sprints ago). Best candidate: 1D conformal tower or gap crossing at q=2-3 on n=6-8 qubits. OVERDUE.
2. **New research direction** — Cylinder thread is complete for exact diag. Options: (a) universality class comparison with S_q Potts and clock operator content, (b) entanglement Hamiltonian on cylinders, (c) return to large-q CFT operator content, (d) Binder cumulant for first-order test.
3. **q=3 or q=5 cylinder entropy** — extend 075b to other q values. Would test whether area-law emergence (S/Ly decreasing) is universal.

## What's Been Ruled Out
- Floating phase in 1D hybrid model (Sprint 067)
- Weakly first-order in 1D at q=10 (Sprint 066)
- Hybrid = clock universality (Sprint 065+067)
- c·x₁ = 1/9 exact (Sprint 063c)
- L=2 in scaling regime for 2D torus (Sprints 068-070)
- DMRG for q≥5 cylinders at chi<50 (Sprint 071b)
- First-order transition on Ly=2 cylinder for q=2,3,5 (Sprints 071c, 072b)
- Universal Ly convergence rate across q (Sprint 073)
- q=3 Ly=4 cylinder with exact diag (Sprint 074a) — 838s/pt
- q=5 Ly≥4 cylinder with exact diag (Sprint 074b)
- Exponential Ly convergence for q=2 (Sprint 075c) — power-law 1/Ly² fits better
- q=2 Ly=5 Lx≥5 with exact diag (Sprint 075a) — 541s/pt

## Key Tools Available
- 2D torus exact diag: L×L up to ~10 sites total
- Cylinder exact diag: Ly=5 Lx≤4 for q=2, Ly=3 Lx≤4 for q=3, Ly=3 Lx≤3 for q=5
- 1D exact diag CPU: n≤8 for q≤5. GPU extends to n≤10.
- DMRG: q≤5 1D chains at n≤24. q=2 cylinder at Lx≤20. q≥5 cylinder IMPRACTICAL.
- IBM QPU: 580s remaining
