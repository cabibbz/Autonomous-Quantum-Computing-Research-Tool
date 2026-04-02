# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 070 — 2D energy derivatives and fidelity susceptibility. q=2 confirmed continuous (d²E~L^0.16, χ_F~L^0.94). q=5 INCONCLUSIVE — only L=2,3 accessible, L=2 out of scaling. No latent heat detected. Need larger lattices.

## Active Research Thread
**2D extension of Potts-clock hybrid — transition nature unknown at q=5.**

2D critical points (square lattice, periodic BC):

| q | g_c (1D) | g_c (2D) | 2D/1D | Sizes | gap×L |
|---|----------|----------|-------|-------|-------|
| 2 | 0.250 | 0.771 | 3.08 | L=2,3,4 | 2.363 |
| 3 | 0.333 | 1.267 | 3.80 | L=2,3 | 6.10 |
| 5 | 0.441 | 1.588 | 3.60 | L=2,3 | 3.50 |

Sprint 070 diagnostics at q=2 (L=3→4, validation):
- d²E/dg² peak scaling: L^0.16 (consistent with α=0, continuous)
- χ_F/N peak scaling: L^0.94 (consistent with continuous ν=1)
- dE/dg smooth across g_c, no latent heat

Sprint 070 at q=5 (L=2,3 only):
- L=2 out of scaling (10-50x off) — no valid L=2→3 scaling
- dE/dg smooth at g_c (no latent heat signal)
- F_min = 0.985 at L=3 (lower than q=2,3, but not diagnostic)
- d²E/dg² peak smaller at q=5 than q=2,3

**Fundamental bottleneck:** q=5 L=4 has dim = 5^16 ≈ 10^11 — infeasible.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Ly=2 cylinder DMRG** — 2D hybrid on cylinder geometry. Ly=2 with varying Lx=4-20. TeNPy can handle d=q² per rung site (d=25 for q=5). Extract entropy scaling, energy convergence, correlation length. Bypasses L_max=3 limitation.
2. **Hardware validation** — 580s QPU remaining. Best candidate: 1D q=2 conformal tower (n=6-8 physical qubits). Clear numerical prediction to test.
3. **2D order parameter** — Compute ⟨δ(s_i,s_j)⟩ across the transition at q=5 L=3. First-order shows discontinuous jump.

## What's Been Ruled Out
- Floating phase in 1D hybrid model (Sprint 067)
- Weakly first-order in 1D at q=10 (Sprint 066)
- Hybrid = clock universality (Sprint 065+067)
- c·x₁ = 1/9 exact (Sprint 063c)
- L=2 in scaling regime for 2D (Sprints 068-070)
- Entropy peak at g_c in 2D at small L (Sprint 069)
- Eigenstate-sum χ_F at dim > 5000 (truncation too severe, Sprint 070b)

## Key Tools Available
- 2D exact diag: L×L up to ~10 sites total (q=5 L=3 in 17-70s/pt, q=2 L=4 in 1-4s/pt)
- 1D exact diag CPU: n≤8 for q≤5. GPU extends to n≤10.
- DMRG: q≤15 at n≤8 (chi=30), q≤5 at n≤24
- IBM QPU: 580s remaining
- results.db: all key measurements from sprints 50-70
