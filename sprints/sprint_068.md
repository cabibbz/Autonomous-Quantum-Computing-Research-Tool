# Sprint 068 — 2D Hybrid Model: Does the Continuous Transition Persist?

## Motivation

The Potts-clock hybrid H = -Σ δ(s_i,s_j) - g Σ(X+X†) has continuous second-order transitions for ALL q=2-10 in 1D (Sprints 053-067). But in 2D:
- Standard S_q Potts: **first-order** for q>4 (rigorously proven)
- Z_q clock: **two BKT transitions** for q>4
- Our hybrid: **UNKNOWN** in 2D

This sprint extends the hybrid model to a 2D square lattice to determine whether the continuous transition persists or becomes first-order.

**Literature search:** No prior study of the Potts-δ coupling + clock-X field hybrid on 2D lattices found. The 2D classical Potts first-order transition (q>4) is proven (Duminil-Copin et al. 2017). Recent work (arXiv:2511.11919) shows extended-range interactions can change transition order. Our hybrid has nearest-neighbor coupling but different symmetry (Z_q vs S_q), which could change the transition order.

## Experiments

### 068a — 2D Hamiltonian Validation (q=2,3)
**Status:** Complete.

Built 2D periodic square lattice Hamiltonian with Potts δ-coupling + clock field (X+X†). Validated on q=2 (Ising) and q=3 with lattices 2×2, 2×3, 2×4, 3×3.

Gap×N crossings found. Strip geometries (2×N) are quasi-1D — crossings depend heavily on lattice shape. Need square lattices with gap×L scaling.

### 068b — 2D Phase Transition on Square Lattices (q=2,3,5)
**Status:** Complete.

Used gap×L scaling on L×L periodic lattices (gap ~ L^{-z} with z=1 at criticality).

**q=2 (2D Ising):** L=2,3,4 (dim up to 65536).
- g_c(2D) ≈ 0.771 (L=3 vs L=4 crossing)
- Gap*L = 2.363 at crossing — **both L=3 and L=4 agree to 4 decimal places**
- Confirms continuous transition, validates method
- 2D/1D ratio: 0.771/0.250 = 3.08

**q=3:** L=2,3 (dim up to 19683).
- g_c(2D) ≈ 1.267
- Gap*L = 6.10 at crossing (L=2 and L=3 agree)
- 2D/1D ratio: 1.267/0.333 = 3.80

**q=5 (KEY TEST):** L=2 (dim=625), L=3 (dim=1.95M, ~25s/pt).
- g_c(2D) ≈ 1.588
- gap2/gap1 = 1.000 everywhere — Z₅ conjugate pair degeneracy maintained in 2D
- 2D/1D ratio: 1.588/0.441 = 3.60
- Only 2 sizes — need first-order diagnostic

**2D/1D g_c ratio:** ~3-4 for all q (consistent with doubling coordination number z=2→4).

### 068c — First-Order vs Continuous Diagnostic
**Status:** [pending]
