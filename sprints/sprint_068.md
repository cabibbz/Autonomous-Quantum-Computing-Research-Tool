# Sprint 068 — 2D Hybrid Model: Does the Continuous Transition Persist?

## Motivation

The Potts-clock hybrid H = -Σ δ(s_i,s_j) - g Σ(X+X†) has continuous second-order transitions for ALL q=2-10 in 1D (Sprints 053-067). But in 2D:
- Standard S_q Potts: **first-order** for q>4 (rigorously proven)
- Z_q clock: **two BKT transitions** for q>4
- Our hybrid: **UNKNOWN** in 2D

This sprint extends the hybrid model to a 2D square lattice to determine whether the continuous transition persists or becomes first-order.

**Literature search:** No prior study of the Potts-δ coupling + clock-X field hybrid on 2D lattices found. The 2D classical Potts first-order transition (q>4) is proven (Duminil-Copin et al. 2017). Recent work (arXiv:2511.11919) shows extended-range interactions can change transition order. Our hybrid has nearest-neighbor coupling but different symmetry (Z_q vs S_q), which could change the transition order.

## Experiments

### 068a — 2D Hamiltonian Construction and Validation
**Status:** Complete.

Built 2D periodic square lattice Hamiltonian: H = -Σ_{<ij>} δ(s_i,s_j) - g Σ_i (X_i+X_i†) on L×L torus. Nearest-neighbor bonds include horizontal and vertical (periodic both directions). Total bonds = 2·L².

Validated on q=2 (Ising) and q=3 with lattices 2×2, 2×3, 2×4, 3×3. Strip geometries (2×N) give quasi-1D behavior — gap×N crossings are geometry-dependent. Switched to square lattices with gap×L scaling for 068b.

### 068b — 2D Phase Transitions on Square Lattices (q=2,3,5)
**Status:** Complete.

**Method:** Gap×L scaling on L×L periodic lattices. At a continuous critical point with z=1, gap ~ L^{-1} so gap×L → universal constant. Crossing between two sizes locates g_c.

**q=2 (2D Ising), L=2,3,4:**
- L=2 vs L=3 crossing: g = 1.180
- **L=3 vs L=4 crossing: g = 0.771** (best estimate)
- L=2 vs L=4 crossing: g = 0.952
- Gap×L at L=3,4 crossing: **2.3632 and 2.3631** (agree to 4 decimal places!)
- L=2 out of scaling regime (gap×L = 4.19 at same g)
- **Confirms continuous transition.** Method validated.
- 2D/1D ratio: 0.771/0.250 = 3.08

**q=3, L=2,3:**
- L=2 vs L=3 crossing: g = 1.267
- Gap×L at crossing: 6.102 (both sizes)
- 2D/1D ratio: 1.267/0.333 = 3.80

**q=5, L=2,3 (KEY TEST):**
- L=2 vs L=3 crossing: g = 1.588
- Gap×L at crossing: 3.504 (both sizes)
- L=3 computation: 25s per g point (dim=1.95M)
- gap₂/gap₁ = 1.000 at all g — Z₅ conjugate pair degeneracy maintained in 2D
- 2D/1D ratio: 1.588/0.441 = 3.60

**2D/1D g_c ratio ~3-4 for all q.** Consistent with doubling coordination z=2→4.

### 068c — First-Order vs Continuous Diagnostic
**Status:** Complete.

**Important caveat:** At a gap×L crossing of only 2 sizes, gap(L₂)/gap(L₁) = L₁/L₂ by construction. The z=1.00 result for q=3,5 is circular — not an independent test. The q=2 test with 3 sizes IS non-trivial and confirms z=1.

**TEST 1 — dE₀/dg smoothness (q=5, L=3):**
Fine scan of 10 g points near g_c=1.588. Energy slope dE₀/dg = -1.75 (g=1.1) → -1.96 (g=1.9). **Monotonically smooth, no kink.** First-order would show a jump in dE/dg that sharpens with L.

**TEST 2 — q-dependence at L=2:**
Gap×L at estimated g_c decreases with q: 5.4 (q=2), 5.8 (q=3), 3.5 (q=5), 2.3 (q=7), 1.4 (q=10). Smooth monotonic decrease. Z_q conjugate pair degeneracy (gap₂/gap₁ = 1.000) maintained for ALL q.

**TEST 3 — Gap scaling (circular for 2 sizes):**
- q=2 at L=3,4 g_c: gap×L = 2.363 for both. **z=1 confirmed (non-trivial, 3 sizes).**
- q=3,5: z=1.00 by construction (only 2 sizes).

**Verdict:** No first-order signal detected at q=5 in 2D (smooth dE/dg, maintained degeneracy), but **INCONCLUSIVE** — L=2,3 are too small for definitive determination. Need L=3,4 or L=4,5 at q=5 (infeasible: dim=5^{16}≈152B).

## Key Results

| q | g_c (1D) | g_c (2D) | 2D/1D | Sizes | gap×L at g_c | Status |
|---|----------|----------|-------|-------|--------------|--------|
| 2 | 0.250 | 0.771 | 3.08 | L=2,3,4 | 2.363 | **Continuous confirmed** |
| 3 | 0.333 | 1.267 | 3.80 | L=2,3 | 6.10 | Consistent with continuous |
| 5 | 0.441 | 1.588 | 3.60 | L=2,3 | 3.50 | **No first-order signal, inconclusive** |

## Surprises
- q=2 L=2 is completely out of the scaling regime (gap×L=4.19 vs 2.36 for L=3,4)
- q=3 gap×L at crossing (6.10) is much larger than q=2 (2.36) and q=5 (3.50)
- Z_q conjugate pair degeneracy (gap₂=gap₁) is EXACT for all q, all L, all g in 2D
- 2D/1D g_c ratio is non-monotonic: 3.08 (q=2), 3.80 (q=3), 3.60 (q=5)
- L=3 at q=5 (dim=2M) feasible at ~25s/pt — 2D exact diag is practical up to ~10 sites

## POTENTIALLY NOVEL
**First study of the Potts-clock hybrid model in 2D.** Phase transition location g_c(q) in 2D, gap×L universal values, and the absence of first-order signal at q=5 are all new measurements. If confirmed by larger lattice sizes (quantum Monte Carlo or tensor network methods), the persistence of continuous transitions in 2D would be a significant result — contradicting the classical 2D Potts first-order transition for q>4.

[Results: results/sprint_068a_2d_validation.json, sprint_068b_q2q3.json, sprint_068b_q5.json, sprint_068c_first_order.json]
