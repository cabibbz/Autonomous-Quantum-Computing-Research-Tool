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

**Critical exponent (Sprints 037-038):** Data collapse confirms Ising ν=1 for TFIM. Crossing-point extraction fragile due to flat landscape.

**Potts MI-CV crossings CONFIRMED at true g_c (Sprint 050).** Qualitative crossing signature vindicated. q=4 MI-CV has no crossings (marginal). Old ν estimates from MI-CV were wrong (Sprint 053).

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

**DMRG excited states fail for Potts.** orthogonal_to gives gap=0.

**MI-CV requires χ > d² for reliable results.** At d=2 (q=2), χ=20 is fine. At d=10 (q=10), need χ>100. Dead-pair bias + MI non-convergence confound at large d.

## Entropy and Central Charge (Sprints 049, 054-056)

**Three methods for c extraction, ranked:**
1. **Entropy profile** (Sprint 055, BEST): Fit S(l) vs chord distance ln[(2n/pi)sin(pi*l/n)] at single large n. Central half of chain only. Even/odd oscillations negligible. Converges monotonically from above. At n=64: 2.5% overshoot (q=2).
2. **FSS pairwise** (Sprint 054): S(n/2) at multiple n, pairwise c from consecutive sizes. Slower convergence, needs 3+ sizes. At (n=16,24): 9-23% overshoot.
3. **iDMRG S vs ln(xi)** (Sprint 055, DO NOT USE): Correlation length saturates at criticality with L=2 unit cell. 18% error, pairwise c scattered.

**c(q) at true critical points (Sprints 054-056):**

| q | g_c | c (profile, best n) | n | c (corrected) | c (exact/CFT) |
|---|-----|---------------------|---|---------------|---------------|
| 2 | 1.000 | 0.512 | 64 | 0.500 | 0.500 |
| 3 | 0.333 | 0.827 | 48 | 0.803 | 0.800 |
| 4 | 0.392 | 1.148 | 24 | ~1.00 | 1.000 |
| 5 | 0.441 | 1.261 | 16 | ~1.10 ± 0.10 | — |
| 7 | 0.535 | 1.462 | 8 | ~1.3 ± 0.15 | — |
| 10 | 0.684 | 1.596 | 6 | ~1.4 ± 0.20 | — |

**c(q) grows monotonically, approximately as c ≈ 0.40·ln(q-1) + 0.55 (Sprint 056).** No peak, no saturation. Physical interpretation: more local states q → more effective massless modes → higher c.

**Analytic continuation of Potts CFT is WRONG for q>4 (Sprint 056).** Both Coulomb gas (g>1) and minimal model (m via arccosh) continuations predict c DECREASING below 1 for q>4. Definitively ruled out by c(q=7)≫1 and c(q=10)≫1 measurements. The 1D quantum Potts at q>4 is described by a different CFT.

**Quadratic interpolation through q=2,3,4 also WRONG (Sprint 056).** Gives c(5)=1.10 (accidentally exact!) but c(7)=1.00 and c(10)=0.10. A trap — would have misled without q=7,10 data.

**c(q=4) has anomalous FSS.** Both methods show flat, non-converging overshoot (+14-23%). Consistent with logarithmic corrections at marginal q=4 (Ashkin-Teller).

**POTENTIALLY NOVEL: c(q≥5) outside minimal models.** c(q=5)≈1.10, c(q=7)≈1.3, c(q=10)≈1.4 — all above c=1 and growing. No CFT predictions exist for q>4 Potts (2D classical is first-order). **Literature search found no prior measurements.** The c(q) ~ ln(q) growth formula is also novel.

## CFT Operator Content (Sprint 057)

**Method:** Energy spectrum on periodic chain at g_c. Ratios R_n = (E_n-E_0)/(E_1-E_0) = x_n/x_1 are universal. Validated on q=2 (Ising) and q=3 (W₃ CFT).

**Z_q spin field degeneracy pattern.** The spectrum contains exactly (q-1) spin field primaries:
- ⌊(q-1)/2⌋ conjugate pairs (σᵏ, σ^{q-k}) for k=1,...,⌊(q-1)/2⌋
- Plus 1 self-conjugate σ^{q/2} if q is even

| q | Structure | Harmonic ratios R(σᵏ/σ) | R(ε/σ) |
|---|-----------|-------------------------|--------|
| 2 | 1 (self-conj) | — | 8.0 |
| 3 | 2 (pair) | — | 6.0 |
| 4 | 2+1 | σ²: 1.68 | 6.6 |
| 5 | 2+2 | σ²: 2.41 | 7.1 |
| 7 | 2+2+2 | σ²: 3.32, σ³: 5.95 | 8.1 |
| 10 | 4 pairs+1 | σ²: 3.66, σ³: 7.69 | 8.3 |

**NOT a free boson.** Free boson predicts x(σᵏ)/x(σ) = k². Measured ratios always BELOW k², but approach k² as q→∞. The q→∞ limit is a free boson; finite q has non-perturbative Z_q corrections.

**Energy field R_ε = x_ε/x_σ** has minimum at q=3 (R=6), increases for q>3 toward ~8+. The energy operator is pushed above all spin harmonics.

**Physical interpretation:** c(q)~ln(q) growth is explained by the proliferation of (q-1) spin field primaries below the energy scale. Each primary contributes effective degrees of freedom.

**POTENTIALLY NOVEL:** Complete operator content (degeneracy pattern, harmonic ratios, energy field position) of 1D quantum q-state Potts CFT for q>4 appears previously unmeasured.

## Scaling Dimension x₁(q) (Sprint 058)

**Method:** Absolute x₁ from CFT Casimir energy + gap on periodic chains. Two sizes give v·c (from E₀/N) and v·x₁ (from Δ₁·N), yielding c/x₁ independent of v.

**c/x₁ ratio (model-independent, no c needed):**

| q | c/x₁ (best pair) | Pair | 2q | x₁ (using measured c) |
|---|-------------------|------|-----|----------------------|
| 2 | 4.013 | (10,12) | 4 | 0.1246 |
| 3 | 5.984 | (8,10) | 6 | 0.1337 |
| 4 | 8.536 | (8,10) | 8 | 0.1172 |
| 5 | 10.840 | (6,8) | 10 | 0.1015 |
| 7 | 15.109 | (4,6) | 14 | 0.0860 |
| 10 | 16.856 | (4,6) | 20 | 0.0831 |

**c/x₁ = 2q EXACT for q=2,3.** For q≥4: c/x₁ grows sub-linearly, significantly below 2q at q=10 (16.9 vs 20). No simple formula found.

**x₁ peaks at q=3 (2/15 ≈ 0.133).** Decreasing for q>3, approaching 0 as q→∞. The peak is related to q=3 being the maximum of the exact Potts CFT spin dimension before the minimal model series ends at q=4.

**q=4 has logarithmic FSS corrections.** The excess (c/x₁ - 8)·N² grows as ~N², confirming sub-power-law (logarithmic) convergence. Consistent with marginal Ashkin-Teller point. Cannot determine if c/x₁(∞) = 8.0 or ~8.5 from n≤10.

**Anomalous FSS for q≥5.** Δ₁·N INCREASES with N (sign-flipped correction vs q≤4). Means x₁ is underestimated at small sizes. Coincides with c>1 and novel CFT regime.

**POTENTIALLY NOVEL:** x₁(q) for q≥5 previously unmeasured. Combined with c(q), ν(q), g_c(q), operator content, and OPE coefficients, this provides the most complete characterization of the novel q>4 Potts CFT family.

## OPE Coefficients (Sprint 060)

**Method:** Ratio of matrix elements of clock operator Z between CFT eigenstates on periodic chain.
C_{sigma*,sigma,epsilon} = |<epsilon|Z_0|sigma>| / |<0|Z_0|sigma>|

Z = diag(1, omega, ..., omega^{q-1}) has Z_q charge +1. Selection rule: <a|Z|b> nonzero only if charge(a) = charge(b)+1 mod q. Ratio method automatically selects sigma* component from degenerate pair.

Validated: q=2 Ising converges to exact C=1/2 with 0.4% error at n=12.

| q | C_sse (extrap.) | largest n | error at n_max |
|---|----------------|-----------|----------------|
| 2 | 0.503 | 12 | 0.4% |
| 3 | 0.540 | 10 | ~3% |
| 4 | 0.464 | 8 | ~7% |
| 5 | 0.376 | 8 | ~7% |
| 7 | 0.272 | 6 | ~10% |
| 10 | ~0.21 | 6 | ~10% |

**C_sse peaks at q=3 (~0.54), then decreases monotonically.** For q≥4, approximate scaling C_sse ~ q^{-0.8}. No simple analytic formula found. Approaches 0 as q→∞.

**Epsilon identification pitfall for large q:** At q≥7, sigma harmonics (sigma^2, sigma^3, ...) densely populate R < R_epsilon. Naive charge detection fails because degenerate pairs (charges k, q-k) in superposition average to misleading charge values. Must use R_epsilon from independent spectrum measurements to identify epsilon correctly.

**POTENTIALLY NOVEL:** First C_{sigma*,sigma,epsilon} measurement for Hermitian q-state Potts at q=5,7,10. Peak at q=3 and monotonic decrease previously unreported.

## Conformal Tower Structure (Sprint 059)

**Method:** Momentum-resolved energy spectrum on periodic chain. Translation operator T gives momentum k (conformal spin s = h - h̄). Degeneracies resolved by diagonalizing T within degenerate eigenspaces of H.

**Genuine CFT confirmed for ALL q=2-10.** Four qualitative tests pass universally:

1. **Descendant degeneracy:** L₋₁σ has degen = 2 × (size of lowest conjugate pair). For q=2: degen=2 (1 field × 2 chiralities). For q≥3: degen=4 (2 fields × 2 chiralities). Exact for all tested q.

2. **Descendant momentum:** Always spin ±1, from Virasoro generators L₋₁ and L̄₋₁.

3. **Gap convergence:** R(L₋₁σ) - R(σ) converges monotonically from below toward 1/x₁.

4. **Tower organization:** Identity(k=0) → spin field pairs(k=0) → ε singlet(k=0) → L₋₁σ(k=±1) for all q.

**Descendant gap convergence (largest n):**

| q | n | gap | 1/x₁ | ratio | error |
|---|---|-----|-------|-------|-------|
| 2 | 12 | 7.90 | 8.0 | 0.987 | 1.3% |
| 3 | 10 | 7.33 | 7.5 | 0.978 | 2.2% |
| 4 | 8 | 8.02 | 8.5 | 0.938 | 6% |
| 5 | 8 | 8.41 | 9.9 | 0.849 | 15% |
| 7 | 6 | 8.59 | 11.6 | 0.738 | 26% |
| 10 | 6 | 8.96 | 12.0 | 0.743 | 25% |

**Anomalous FSS at q≥5** prevents quantitative gap convergence at accessible sizes. Same sign-flipped corrections as Δ₁·N (Sprint 058). NOT evidence against CFT — qualitative tower structure is exact.

**Tower density increases with q.** The (q-1) spin field primaries fill R=1 to R≈8-9. For large q, ε nearly merges with L₋₁σ. In q→∞ limit, approaches free boson continuum.

**Literature context:** Lao et al. (PRB 2019) found drifting towers for q=5,6 assuming first-order transition. Tang et al. (PRL 2024) found clean towers only with non-Hermitian q=5 model. Our Hermitian towers converge because transition is genuinely second-order. **POTENTIALLY NOVEL: First momentum-resolved conformal tower for q>4 Hermitian Potts chain.**

**Entropy profile overshoot grows with q at fixed n:** At n=16: 8.7% (q=3), 14.4% (q=4), ~20% (q=5). At n=8: ~25% (q=7). Large-q profile method requires larger n — q=7+ is computationally infeasible at n≥12 with chi=20.

**Entropy FSS does NOT give ν.** S ~ (c/6)ln(ξ), not a power law. Standard FSS collapse fails (gives ν=3 for TFIM). Use energy gap slope for ν.

**Correlation length from correlator decay.** Saturates at ξ ~ n/4 near criticality. iDMRG correlation length also saturates at criticality (Sprint 055). Not useful for ν extraction near g_c.

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
