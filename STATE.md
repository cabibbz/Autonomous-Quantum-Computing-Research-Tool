# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 068 — 2D hybrid model: first study. g_c found for q=2,3,5. No first-order signal at q=5 but inconclusive (only L=2,3).

## Active Research Thread
**2D extension of Potts-clock hybrid — promising but inconclusive.**

2D critical points (square lattice, periodic BC):

| q | g_c (1D) | g_c (2D) | 2D/1D | Sizes | gap×L | Verdict |
|---|----------|----------|-------|-------|-------|---------|
| 2 | 0.250 | 0.771 | 3.08 | L=2,3,4 | 2.363 | Continuous confirmed |
| 3 | 0.333 | 1.267 | 3.80 | L=2,3 | 6.10 | Consistent |
| 5 | 0.441 | 1.588 | 3.60 | L=2,3 | 3.50 | No 1st-order signal |

Key finding: Z_q conjugate pair degeneracy (gap₂=gap₁) is exact in 2D for all q.
Caveat: L=2,3 only — too small for conclusive first-order test. Need QMC or tensor networks.

1D hybrid universality (Sprints 053-067, fully characterized):

| q | g_c | c | ν | x₁ | c·x₁ | C_sse |
|---|-----|---|---|-----|------|-------|
| 2 | 0.250 | 0.500 | 1.00 | 0.125 | 0.062 | 0.50 |
| 3 | 0.333 | 0.800 | 0.86 | 0.133 | 0.107 | 0.54 |
| 5 | 0.441 | ~1.10 | 0.83 | ~0.10 | 0.112 | 0.38 |
| 10 | 0.684 | ~1.40 | 1.12? | ~0.083 | 0.116 | ~0.21 |

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **2D with larger lattices** — Use DMRG on 2D strips (cylinder geometry) at q=5 to get more sizes. Ly=3 cylinder with Lx=4,6,8 could test gap scaling. Also: quantum Monte Carlo (stochastic series expansion) is sign-problem-free for this Hamiltonian.
2. **Hardware validation** — QPU budget 580s unspent. Strongest prediction: q=2 2D Ising g_c=0.771 on a 2×2 plaquette (4 qubits). Or 1D q=2 conformal tower on real hardware.
3. **Entanglement entropy in 2D** — Compute entanglement entropy at 2D g_c to extract effective c. Does the area-law→log-law transition occur? Compare with 1D c values.

## What's Been Ruled Out
- Floating phase in 1D hybrid model (Sprint 067)
- Weakly first-order in 1D at q=10 (Sprint 066)
- Hybrid = clock universality (Sprint 065+067)
- c·x₁ = 1/9 exact (Sprint 063c)
- L=2 in scaling regime for 2D Ising (Sprint 068: gap×L=4.19 vs 2.36)

## Key Tools Available
- 2D exact diag: L×L up to ~10 sites total (q=5 L=3 in 25s/pt, q=3 L=3 in 0.2s/pt)
- 1D exact diag CPU: n≤8 for q≤5. GPU would extend to n≤10.
- DMRG: q≤15 at n≤8 (chi=30), q≤5 at n≤24
- IBM QPU: 580s remaining
- results.db: all key measurements from sprints 50-68
