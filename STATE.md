# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 077 — S_q Potts at q=7,10: g_c found, NO first-order signal at n≤8. S_q g_c(q=7)=0.144, g_c(q=10)=0.101. Conformal tower at q=7 has 6-fold degeneracy (S₇) vs hybrid 2-fold (Z₇). R_ε ratio 2.3x between models. Δ·N stable <5% for all q=5,7,10 — complex CFT regime extends beyond n=8.

## Active Research Thread
**S_q Potts vs hybrid universality comparison (extended from Sprint 076).**

S_q Potts critical points (all from gap×N crossings):

| q | g_c (S_q) | g_c (hybrid) | g_c ratio | 1st exc. degen | R_ε (S_q) | R_ε (hybrid) |
|---|-----------|-------------|-----------|----------------|-----------|-------------|
| 5 | 0.200 | 0.441 | 2.21 | 4 vs 2 | ~3.9 | ~7.1 |
| 7 | 0.144 | 0.535 | 3.72 | 6 vs 2 | 3.48 | 8.07 |
| 10 | 0.101 | 0.684 | 6.77 | — | — | — |

Both models look CFT-like at accessible sizes. First-order behavior (predicted for S_q q>4) not visible at n≤8.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **S_q Potts large-n via DMRG** — S_q Potts at q=7,10 looks CFT-like at n≤8. DMRG at n=20-50 could reveal first-order regime (ξ* crossover). Most impactful: test if walking length is finite.
2. **Hardware validation** — 580s QPU remaining, last used Sprint 025 (52 sprints ago). Best candidate: conformal tower degeneracy (4-fold vs 2-fold at q=5) on n=6-8 qubits. OVERDUE.
3. **S_q Potts g_c(q) formula** — Have g_c at q=3,5,7,10. Fit to find scaling law analogous to hybrid's g_c ≈ (1/5)√(q-1)+1/20.

## What's Been Ruled Out
- Floating phase in 1D hybrid model (Sprint 067)
- Weakly first-order in 1D hybrid at q=10 (Sprint 066)
- Hybrid = clock universality (Sprint 065+067)
- c·x₁ = 1/9 exact (Sprint 063c)
- L=2 in scaling regime for 2D torus (Sprints 068-070)
- DMRG for q≥5 cylinders at chi<50 (Sprint 071b)
- Exponential Ly convergence for q=2 (Sprint 075c)
- S_q = hybrid universality at q≥4 (Sprints 076-077)
- First-order S_q Potts at q=5,7,10 at n≤8 (Sprints 076-077) — walking behavior
- "69% gap collapse" at q=7 S_q (Sprint 076c) — artifact of wrong g_c

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤5. GPU: n≤10 for q≤5, n≤8 for q=7, n≤6 for q=10.
- S_q Potts Hamiltonian: implemented (exp_076a, exp_077a). g_c known for q=3,5,7,10.
- DMRG: q≤5 1D chains at n≤24. q=2 cylinder at Lx≤20.
- IBM QPU: 580s remaining
