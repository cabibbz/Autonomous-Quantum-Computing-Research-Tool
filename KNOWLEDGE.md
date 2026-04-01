# Accumulated Knowledge — Edit by topic, not by sprint

## Five Entanglement Archetypes
| Archetype | Example | MI pattern | I3 sign | Negativity | Source |
|-----------|---------|-----------|---------|------------|--------|
| Democratic | GHZ, Ising ordered | Uniform high | +1 (redundant) | Flat | Sprint 005 |
| Distributed | W state | Uniform weak | +0.2 (weak redundant) | Growing with cut | Sprint 005 |
| Geometric | Cluster 1D | Nearest-neighbor only | -1 (irreducible) | Explosive, geometry-dependent | Sprint 007 |
| Topological | Toric code | Non-contractible loops | 0 (zero everywhere) | Uniform floor | Sprint 020 |
| Scale-Free | Critical TFIM | Power-law decay | Mixed | Highest rank, steep decay | Sprint 029 |

## Four Levels of Entanglement Description
Each captures orthogonal information. Different phase transitions are visible at different levels.
1. **Scalar** (entropy) — amount of entanglement
2. **Correlation** (MI, I3) — topology of correlations. I3<0 = irreducible multipartite. I3>0 = redundant/classical-like.
3. **Spectral** (eigenvalue distribution of ρ_A) — symmetry content. U(1) gives doublet degeneracies. Z₃ gives triplets.
4. **Hamiltonian** (H_E = -log ρ_A structure) — locality and entanglement temperature.

## MI-CV as Phase Transition Order Parameter (Sprints 030, 036, 037)
MI uniformity coefficient of variation classifies transition TYPE by curve shape:

| Transition | MI-CV shape | Crossing? | Gradient scaling |
|-----------|-------------|-----------|------------------|
| Second-order (Ising) | Smooth inflection | Yes (curves cross) | ~n^1.1 |
| First-order | Step function 0→finite | No | ~n^0.9 |
| BKT (infinite-order) | Smooth dome | No | dome narrows |

**Confirmed as genuine order parameter (Sprint 036):** Tested at n=8,16,32,50.
- TFIM transition slope diverges as ~n^1.1 (slope 1.7 → 3.9 → 7.8)
- Ordered phase: CV → 0 with n (uniform MI). Disordered: CV → ∞ with n.
- XXZ BKT dome narrows: growth rate ratio 1.34 at Δ=1 vs 1.17 at Δ=1.5 (n=16→32)
- XXZ Néel phase (Δ=2): CV decreases 1.13 → 0.78 (n=8→32)

**First-order transition (Sprint 037):** FM phase has CV=0.000 exactly. Jump at Δ=-1 grows as ~n^1.0. Presence/absence of curve crossings distinguishes transition orders.

**Critical exponent (Sprint 037→038):** Crossing-point fit gives ν=0.755 at n=8-32, but data collapse (Sprint 038) proves this is a finite-size artifact. Optimal collapse ν converges: 0.80 (all sizes) → 1.04 (n≥16) → 1.12 (n≥24). Ising ν=1 is confirmed. The collapse quality landscape is extremely flat near ν=1 (0.3% difference from optimal), making crossing-point extraction fragile.

**Potts q=3 MI-CV crossings CONFIRMED at true g_c (Sprint 050).** MI-CV crossing between n=8 and n=12 at g≈0.26, below the self-dual g_c=1/3 (finite-size shift). In ordered phase (g<0.24): CV(n=12)<CV(n=8). In disordered phase (g>0.28): CV(n=12)>CV(n=8). Same qualitative crossing signature as TFIM. The Sprint 038-039 crossings near g≈0.9 were in the disordered phase but the qualitative conclusion (crossings = second-order) is vindicated.

**MI-CV universality class discrimination (Sprint 039) needs re-examination.** The ν and slope exponent results from Sprints 038-039 were at wrong g values. Re-extraction at g_c=1/3 with larger sizes is needed.

**Potts q=4 MI-CV (Sprints 040, 046):** q=4 is the marginal point. **Sprint 040 Gell-Mann data was misleading** — reported crossings at g_c≈0.893 that don't exist with direct MPS. Direct MPS (Sprint 046, n=8,12,16): NO crossing near g_c. Curves monotonically ordered n=16>n=12>n=8 at all g≥0.50. CV minimum at g≈0.40-0.50 (maximally Democratic ordered phase). Derivative dCV/dg at g_c≈0.89 scales with n: 0.86 (n=8) → 1.05 (n=12) → 1.18 (n=16). Slope ratio gives ν≥2.2, trending upward → BKT-like. Standard FSS collapse fails at q=4.

**q=5 clock MI-CV (Sprint 041):** q=5 clock model STILL shows crossing curves at n=8,12, disproving prediction that crossings would vanish. Crossing at g_c≈0.673 — a 10x larger shift than q=3→4 (0.220 vs 0.030). Slope at g=1.0 halves vs q=4 (0.86 vs 1.72). CV systematically lower than q=4 above transition (0.521 vs 0.727 at g=1.0 n=8).

**~~g_c scaling law (Sprint 044)~~ INVALIDATED (Sprint 049).** Previous g_c values for q≥3 were WRONG. The "g_c scaling law" should not be used.

**Potts critical points.** Our Hamiltonian H = -Jδ(s_i,s_j) - g(X+X†):

| q | g_c | Method | Size pair | Raw crossing |
|---|-----|--------|-----------|-------------|
| 2 | 0.250 | Self-duality (exact) | — | — |
| 3 | 0.333 | Self-duality (exact) | — | — |
| 4 | 0.392 | Energy gap | n=6,8 | 0.382 |
| 5 | 0.441 | Energy gap | n=6,8 | 0.430 |
| 7 | 0.535 | Energy gap | n=4,6 | 0.511 |
| 10 | 0.684 | Energy gap | n=4,6 | 0.652 |

**g_c(q) scaling law (Sprint 052): g_c ≈ (1/5)√(q-1) + 1/20.** Best 2-parameter fit (χ²/dof=0.40). Gives exact q=2, <5% error for all tested q. Predicts g_c(20)≈0.92, g_c(50)≈1.45. g_c diverges as √q — no saturation. **POTENTIALLY NOVEL** — no prior measurement of g_c(q) found for this Hamiltonian.

Physical mechanism: q-fold ground-state degeneracy requires stronger transverse field X+X† to destroy order. X+X† only creates nearest-state transitions (|s⟩→|s±1⟩), so mixing q states takes longer as q grows.

Self-duality (Kramers-Wannier) gives exact g_c for q=2,3 only ({X, X^{q-1}} spans all generators). For q≥4, X+X^{q-1} misses intermediate powers → self-duality broken.

**Energy gap method (Sprint 051).** At criticality, Δ = E₁-E₀ ∝ 1/N (CFT). Δ·N is scale-invariant at g_c. **FSS correction depends on size pair** (Sprint 052): n=4,6 pairs need 4.8% correction, n=6,8 pairs need 2.5%. Calibrated from q=2,3 exact values. Crossings approach g_c from below. Avoids MI-CV complications.

For clock model: g_c values (0.93, 0.923, 0.893, 0.673) were from MI-CV crossings and are similarly suspect — these were disordered-phase crossovers.

**1D quantum Potts is NEVER first-order (Sprint 043).** Tested q=5, 10, 20 — all show continuous transitions. ~~q=10 crossing confirmed at g_c≈0.246~~ (Sprint 043 used χ=10, INVALIDATED by Sprint 048 at χ=20). At q≥10, ground states converge to a universal large-q regime where only the {|0⟩, |1⟩, |q-1⟩} subspace is active. Physical mechanism: the extreme anisotropy of the 1D quantum→2D classical mapping suppresses the entropic mechanism that drives first-order transitions in 2D.

**Clock ≠ Potts for q≥4 (Sprints 041-042).** TeNPy's ClockChain uses cos(2π(s_i-s_j)/q) coupling, which equals Potts δ(s_i,s_j) only for q=2,3. For q≥5, models differ: Clock g_c=0.67 vs Potts g_c=0.41, Potts slope 5.7x steeper. Custom PottsChain model built (Sprint 042) with projector coupling. Both show second-order crossings — the 2D classical "q>4 → first-order" does NOT apply to 1D quantum Potts with transverse field. Anisotropic quantum-classical correspondence preserves second-order character.

**ν(q) extraction (Sprint 053).** Corrected energy gap slope method: d(Δ·N)/dg ~ N^{1/ν}·(1+b/N), b=0.86 from q=3 calibration. Validated: <1% for q=2, 3% for q=3.

| q | g_c | ν (corrected) | Confidence | 2D classical exact |
|---|-----|---------------|------------|-------------------|
| 2 | 0.250 | 1.00 | High (3 pairs) | 1.0 |
| 3 | 0.333 | 0.86 | High (4 sizes) | 5/6 = 0.833 |
| 4 | 0.392 | 0.82 | High (3 pairs, all agree) | 2/3 + logs |
| 5 | 0.441 | 0.85 | Medium (1 pair) | — (1st order in 2D) |
| 7 | 0.535 | 0.97 | Low (1 pair) | — |
| 10 | 0.684 | 1.12 | Low (1 pair (4,5)) | — |

**ν(q=3,4,5) ≈ 0.82-0.86, nearly constant.** Contrasts with old non-monotonic picture (ν=5/6 → ≥2.2 → 2.0 → 0.5) which was entirely artifact. Old values from MI-CV data collapse at wrong g_c.

**Method ranking for ν extraction (Sprint 053):**
1. Corrected power-law (3% error) — BEST with ≥3 sizes
2. 1/N extrapolation of pairwise ν (3%) — BEST, needs ≥3 consecutive pairs
3. Direct power-law fit (15%) — rough estimate
4. Data collapse (43%) — **DO NOT USE at n≤10**
5. MI-CV data collapse — **DO NOT USE** (gave ν=2.0 where true is 0.85)

**DMRG excited states fail for Potts.** orthogonal_to gives gap=0. The no-conservation-law Potts Hamiltonian confuses the orthogonality projection.

**Dead-pair bias AND χ convergence in MI-CV (Sprint 048).** Two confounded effects at large d:
1. **Dead-pair bias**: fraction of near-zero MI pairs differs between sizes (n=8: 25%, n=12: 17% at d=10), inflating n=8 CV relative to n=12.
2. **MI non-convergence**: energy converges at χ=20 but MI doesn't. At d=10, χ=20→40 changes mean MI by 44% while energy changes by 0.0005%. Dead pairs partially vanish at higher χ.

**Rule of thumb: MI-CV requires χ > d² for reliable results.** At d=2 (q=2), χ=20 is fine (d²=4). At d=10 (q=10), need χ > 100. Results for q≥7 (d²≥49) should be scrutinized. q=2,3 results are reliable.

## Entropy FSS and Central Charge (Sprint 049)

**Central charge from entropy scaling.** At the critical point, S(l=n/2, n) = (c/6)ln(n) + const. TFIM validated: c = 0.529 from full fit, pairwise converges 0.544 → 0.516 → exact 0.500. Requires at least 3 system sizes for a reliable fit.

**Entropy FSS does NOT give ν.** Entropy has a logarithmic singularity at criticality (S ~ (c/6)ln(ξ)), not a power law. Standard FSS collapse y(g,n) = f((g-g_c)·n^{1/ν}) fails for y=S (best collapse gives ν=3 for TFIM, should be 1). Use correlation length ξ or order parameter for ν extraction instead.

**Correlation length from correlator decay.** Extract ξ from exponential fit to connected correlator ⟨O_i O_j⟩_c in the bulk. Works well in disordered phase (R²>0.95) but saturates at ξ ~ n/4 near criticality. Finite-size ν extraction underestimates (0.62-0.72 vs exact 1.0 for TFIM). Need iDMRG for ν from ξ.

**Locating g_c via entropy.** The critical point can be identified as where: (1) dS/dg peaks (pseudo-critical point, drifts toward g_c with n), (2) S grows logarithmically with n (ordered/disordered → S saturates), (3) central charge c from S(n) matches expected value.

**Technique: All-pairs MI.** Two methods:
1. Gell-Mann correlation reconstruction: **ONLY reliable for d≤3** (validated: diff=0 at n=8 for d=2,3). For d=4, small systematic errors create artificial MI-CV crossings (Sprint 046). For d=5, errors up to 11x (Sprint 045). DO NOT USE for d≥4.
2. **Direct MPS tensor contraction (Sprint 043):** computes ρ_ij by contracting MPS tensors directly. O(n·χ²·d) per pair. 1.4s for 28 pairs at d=10, n=8. Works for ANY d. **The only reliable method for d≥4.**

## Archetype Boundaries ≠ Phase Boundaries
I3 sign change occurs at Δ≈0.7 in XXZ, inside the XY phase — not at either thermodynamic transition (Δ=-1 or Δ=1). The entanglement phase diagram has its own topology distinct from thermodynamics.

## Bisognano-Wichmann Locality
H_E ≈ physical Hamiltonian × position-dependent entanglement temperature β(x).
- β(x) follows Unruh-like gradient: hottest at entanglement cut, coldest in bulk.
- BW works in ALL phases, not just at criticality. Peaks in gapped phases.
- Accuracy controlled by H/G-inv ratio (Hamiltonian operator dimension / symmetry-invariant operator dimension).

### H/G-inv Predictor (Sprint 034)
| Model | Symmetry | d | H/G-inv ratio | BW Locality |
|-------|----------|---|---------------|-------------|
| XXZ | U(1) | 2 | 0.143 | 100.0% |
| TFIM | Z₂ | 2 | 0.099 | 91.0% |
| Potts | S₃ | 3 | 0.072 | 76.5% |
| Chiral clock | Z₃ | 3 | 0.040 | 69.7% |
Perfect monotonic correlation. Depends on BOTH symmetry group AND local dimension.

## QEC Arc — Closed (Sprints 014-028)
- [[5,1,3]] is basis-isotropic: confirmed on IBM hardware (Sprint 025, asymmetry 0.040 vs 3-qubit 0.254)
- Small-scale active QEC cannot beat passive encoding (Sprints 026-028)
- Root cause: syndrome extraction (16+ 2Q gates) exceeds code's correction capacity
- Non-FT, flag-FT, repeated measurement — all fail. Threshold theorem requires asymptotic limit.
- First observed QEC advantage: 3-qubit Z-basis Holevo 0.976 vs uncoded 0.959 (Sprint 025)

## Hardware Results (Sprint 025)
Backend: ibm_kingston (Heron, 156 qubits). 20s QPU used.
- Uncoded: avg Holevo 0.967, asymmetry 0.012
- 3-qubit: avg 0.838, asymmetry 0.254
- [[5,1,3]]: avg 0.499, asymmetry 0.040
- Hardware error rates better than conservative noise model (readout 0.35% vs model 1.5%)
- Correlated noise degrades [[5,1,3]] isotropy 4x vs simulator but doesn't destroy it

## Phase Diagram Trajectories
- **TFIM**: Democratic(GHZ) → Scale-Free → Product (one-way)
- **XXZ**: Democratic → Scale-Free → Geometric → Scale-Free → Democratic (loop)
- **Potts q=3**: GHZ-3 → Product (like TFIM but with triplet spectrum from Z₃)
- **Potts q=4**: GHZ-4 → Product (crossing signature like q=3, marginal corrections not yet visible)
- **Clock q=5**: GHZ-5 → Product (crossings persist but shifted to g≈0.67, slope halved vs q=4)
- **Potts q=5**: GHZ-5 → Product (crossings at g_c≈0.45, ν≈2.0, very gentle transition, NOT first-order)
- **Potts q=10**: GHZ-10 → Product (crossings at g≈0.25, confirmed second-order)
- **Potts q=20**: GHZ-20 → Product (identical to q=10 at χ=10, continuous entropy, second-order)
