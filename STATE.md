# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 076 — S_q Potts vs Hybrid head-to-head at q=5. S_q g_c=0.200, hybrid g_c=0.438. Degeneracy structure is the sharpest discriminator: S_q has 4-fold (S₅), hybrid has 2-fold (Z₅) with split at x₃=2.41. Both look like CFTs at n≤8 — no first-order signal for S_q at accessible sizes.

## Active Research Thread
**Direct universality class comparison — S_q Potts vs hybrid vs clock (new thread).**

Key results from Sprint 076:

| Property | S_q Potts q=5 | Hybrid q=5 |
|----------|--------------|------------|
| g_c | 0.200 | 0.438 |
| Field symmetry | S₅ (full permutation) | Z₅ (cyclic) |
| 1st excited degeneracy | 4-fold | 2-fold |
| x₃ | 1.000 (degenerate) | 2.413 (split) |
| Δ·N drift (n=4→8) | 4.8% | 0.46% |
| FSS corrections | 10x larger | Nearly converged |
| z_eff | 1.071 | 1.007 |

Prior Ly convergence dataset (q=2, Sprints 071-075): g_c(Ly) = 0.771 - 1.285/Ly^2.03 (power-law).

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **S_q Potts at q=7 — proper g_c measurement.** The q=7 gap collapsed 69% at estimated g_c=0.12 — need proper gap crossing scan to find true g_c and test first-order at larger q.
2. **Hardware validation** — 580s QPU remaining, last used Sprint 025 (51 sprints ago). Best candidate: gap crossing or conformal tower at q=2-3 on n=6-8 qubits. OVERDUE.
3. **S_q Potts large-n via DMRG** — S_q Potts at q=5 looks CFT-like at n≤8. DMRG at n=20-50 could reveal whether first-order nature emerges at larger N (ξ >> 8?).

## What's Been Ruled Out
- Floating phase in 1D hybrid model (Sprint 067)
- Weakly first-order in 1D hybrid at q=10 (Sprint 066)
- Hybrid = clock universality (Sprint 065+067)
- c·x₁ = 1/9 exact (Sprint 063c)
- L=2 in scaling regime for 2D torus (Sprints 068-070)
- DMRG for q≥5 cylinders at chi<50 (Sprint 071b)
- Exponential Ly convergence for q=2 (Sprint 075c) — power-law 1/Ly² fits better
- First-order S_q Potts signature at q=5 n≤8 (Sprint 076c) — both models look CFT-like
- S_q = hybrid universality at q≥4 (Sprint 076b) — different degeneracy structures

## Key Tools Available
- 2D torus exact diag: L×L up to ~10 sites total
- Cylinder exact diag: Ly=5 Lx≤4 for q=2, Ly=3 Lx≤4 for q=3, Ly=3 Lx≤3 for q=5
- 1D exact diag CPU: n≤8 for q≤5. GPU extends to n≤10.
- DMRG: q≤5 1D chains at n≤24. q=2 cylinder at Lx≤20. q≥5 cylinder IMPRACTICAL.
- S_q Potts Hamiltonian: now implemented (exp_076a). Same Potts coupling, S_q-symmetric field.
- IBM QPU: 580s remaining
