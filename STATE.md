# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 069 — 2D entanglement entropy at criticality. Ordered phase has S=ln(q) exactly (GHZ-like). Critical entropy follows area law but is much smaller than ln(q) at accessible sizes. Entropy peak is NOT at g_c.

## Active Research Thread
**2D extension of Potts-clock hybrid — entropy characterized, need larger lattices.**

2D critical points (square lattice, periodic BC):

| q | g_c (1D) | g_c (2D) | 2D/1D | Sizes | gap×L | S_crit(w=1,L=3) |
|---|----------|----------|-------|-------|-------|------------------|
| 2 | 0.250 | 0.771 | 3.08 | L=2,3,4 | 2.363 | 0.392 |
| 3 | 0.333 | 1.267 | 3.80 | L=2,3 | 6.10 | 0.235 |
| 5 | 0.441 | 1.588 | 3.60 | L=2,3 | 3.50 | 0.321 |

Key findings Sprint 069:
- Ordered phase: S = ln(q) exactly for ALL q, L (GHZ-like symmetry-preserving ground state)
- Critical entropy: area law, S/boundary ≈ 0.04-0.07, non-monotonic in q
- Entropy peak is in ordered phase (NOT at g_c) for L≤4 — needs L > ln(q)/(2α) ≈ 10-16 to cross over
- Sharp entropy drop at g ≈ 0.6·g_c, not at g_c itself

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **2D cylinder DMRG** — Ly=2 cylinder (d=q per site), vary Lx. Can get ground state energy and entropy (not gap). Extract area-law coefficient α(q) and Casimir energy. Feasible for q=5 with Ly=2, Lx=4-12.
2. **Hardware validation** — 580s QPU unspent. Best prediction: 1D q=2 conformal tower on real hardware (n=6-8 qubits). Or 2D q=2 on 2×2 plaquette (4 qubits).
3. **2D energy derivatives** — dE₀/dg and d²E₀/dg² at g_c for q=5 L=3. The second derivative should diverge at a continuous critical point. Compare q=2 (known) with q=5.

## What's Been Ruled Out
- Floating phase in 1D hybrid model (Sprint 067)
- Weakly first-order in 1D at q=10 (Sprint 066)
- Hybrid = clock universality (Sprint 065+067)
- c·x₁ = 1/9 exact (Sprint 063c)
- L=2 in scaling regime for 2D (Sprints 068-069)
- Entropy peak at g_c in 2D at small L (Sprint 069: peak is in ordered phase)

## Key Tools Available
- 2D exact diag: L×L up to ~10 sites total (q=5 L=3 in 17s/pt, q=2 L=4 in 0.8s/pt)
- 1D exact diag CPU: n≤8 for q≤5. GPU extends to n≤10.
- DMRG: q≤15 at n≤8 (chi=30), q≤5 at n≤24
- IBM QPU: 580s remaining
- results.db: all key measurements from sprints 50-69
