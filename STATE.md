# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 072 — q=3 Ly=2 cylinder + Ly=3 extension for q=2. g_c(cyl, q=3)=0.565 (3 crossing pairs). Cyl/1D ratio monotonically decreasing: 1.80→1.69→1.62 for q=2,3,5. Ly=3 cylinder for q=2: g_c=0.655 (4 pairs), 63.7% of the way from Ly=2 to 2D. All transitions smooth — no first-order signal.

## Active Research Thread
**Cylinder geometry characterization of hybrid model — now well-established.**

Complete Ly=2 cylinder dataset:

| q | g_c(1D) | g_c(cyl,Ly=2) | g_c(2D) | cyl/1D | Crossing pairs |
|---|---------|---------------|---------|--------|----------------|
| 2 | 0.250 | 0.451 | 0.771 | 1.80 | (4,5),(5,6),(6,7) |
| 3 | 0.333 | 0.565 | 1.267 | 1.69 | (3,4),(4,5),(5,6) |
| 5 | 0.441 | 0.714 | 1.588 | 1.62 | (3,4) |

Ly convergence for q=2:

| Ly | g_c | Progress to 2D |
|----|-----|---------------|
| 1 | 0.250 | 0% |
| 2 | 0.451 | 38.6% |
| 3 | 0.655 | 77.7% |
| ∞ | 0.771 | 100% |

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU remaining. Best candidate: 1D q=2 conformal tower (n=6-8 qubits). Clearest numerical prediction to test. Every ~10 sprints should consider QPU use (last was Sprint 025, 47 sprints ago).
2. **Cylinder Ly=4 for q=2** — dim=2^{4*Lx}=16^Lx. Lx=3: 4096, Lx=4: 65536, Lx=5: 1M. Would test if Ly convergence accelerates (63.7% at Ly=3, expect ~85% at Ly=4).
3. **q=3 Ly=3 cylinder** — dim=3^{3*Lx}=27^Lx. Lx=3: 19683, Lx=4: 531k. Compare Ly convergence between q=2 and q=3.

## What's Been Ruled Out
- Floating phase in 1D hybrid model (Sprint 067)
- Weakly first-order in 1D at q=10 (Sprint 066)
- Hybrid = clock universality (Sprint 065+067)
- c·x₁ = 1/9 exact (Sprint 063c)
- L=2 in scaling regime for 2D torus (Sprints 068-070)
- DMRG for q≥5 cylinders at chi<50 (Sprint 071b)
- Entropy peak as g_c proxy on finite-width cylinders (Sprint 071a)
- First-order transition on Ly=2 cylinder for q=2,3,5 (Sprints 071c, 072b)

## Key Tools Available
- 2D torus exact diag: L×L up to ~10 sites total
- Cylinder exact diag: Ly=2 Lx≤6 for q=3 (531k), Ly=3 Lx≤7 for q=2 (2M)
- 1D exact diag CPU: n≤8 for q≤5. GPU extends to n≤10.
- DMRG: q≤5 1D chains at n≤24. q=2 cylinder at Lx≤20. q≥5 cylinder IMPRACTICAL.
- IBM QPU: 580s remaining
