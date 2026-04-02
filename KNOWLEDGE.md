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

## ⚠ Model Identity — READ THIS (External Review, April 2026)

Our Hamiltonian H = -Jδ(s_i,s_j) - g(X+X†) is a **Potts-clock hybrid**: Potts coupling + clock transverse field. It is **not the standard quantum Potts model** from the literature.

| Model | Coupling | Transverse field | Symmetry | q>4 transition | Floating phase |
|-------|----------|-----------------|----------|---------------|----------------|
| S_q Potts | δ(s_i,s_j) | Σ_{k=1}^{q-1} X^k | S_q | **First-order** | No |
| Z_q clock | cos(2π(s_i-s_j)/q) | X + X† | Z_q | **Two BKT** | **Yes** (Sprint 067c) |
| **Our hybrid** | δ(s_i,s_j) | X + X† | Z_q | Continuous 2nd order | **No** (Sprint 067) |

For q=2,3 all three are equivalent. For q≥4 they differ.

**Sprint 065 CONFIRMED: hybrid ≠ clock universality class.** Three independent probes at q=5:
- c/x₁: hybrid 10.77 vs clock 9.43 (12% diff, slowly shrinking but nonzero)
- ν: hybrid 0.83 (finite, power-law) vs clock ~2+ (diverging, BKT)
- Clock q=7 has NO gap crossing (BKT); hybrid has clear crossing (second-order)

**Key literature (search before claiming novelty):**
- **Gorbenko, Rychkov & Zan (JHEP 2018, SciPost 2018):** Complex CFT for q>4 S_q Potts. Predicts complex c (e.g., c ≈ 1.138 ± 0.021i for q=5).
- **Ma & He (PRB 99, 195130, 2019):** Measured effective c from entropy for q=5,6,7 on Hermitian S_q chain.
- **Tang et al. (PRL 133, 076504, 2024):** Non-Hermitian q=5 S_q Potts. 11 scaling dimensions, 9 OPE coefficients.
- **Sun, Luo & Chen (arXiv:2006.11361, 2020):** Z_q clock has two BKT transitions for q>4. Continuous transitions known for clock-type fields.
- **Jacobsen & Wiese (PRL 133, 077101, 2024):** All S_q Potts exponents for q>4 via analytic continuation.

**Sprint 076 CONFIRMED: hybrid ≠ S_q Potts universality class.** Direct head-to-head at q=5:
- Degeneracy: S_q has 4-fold (S₅ permutation), hybrid has 2+2 split (Z₅ conjugate pairs, x₃=2.41)
- g_c: S_q 0.200, hybrid 0.438 (field strength ratio 6.47x → g_c ratio 2.19x)
- FSS corrections: S_q 10x larger (Δ·N drifts 4.8% vs 0.46% over n=4-8)
- z_eff: S_q 1.071, hybrid 1.007 (both power-law at n≤8)

**Surprise: S_q Potts looks CFT-like at n≤8 for ALL q=5,7,10 (Sprints 076-077).** Clean conformal towers, exact degeneracies, Δ·N stable to <5%. The "walking" correlation length ξ* exceeds n=8 even at q=10. First-order nature requires much larger systems. Complex CFT (Gorbenko et al.) describes S_q behavior at moderate sizes.

**S_q Potts critical points (Sprint 076-077):**

| q | g_c (S_q) | g_c (hybrid) | g_c ratio | S_q 1st degen | Hybrid 1st degen |
|---|-----------|-------------|-----------|---------------|-----------------|
| 3 | 0.333 | 0.333 | 1.00 | identical | identical |
| 5 | 0.200 | 0.441 | 2.21 | 4-fold (S₅) | 2-fold (Z₅) |
| 7 | 0.144 | 0.535 | 3.72 | 6-fold (S₇) | 2-fold (Z₇) |
| 10 | 0.101 | 0.684 | 6.77 | (q-1)-fold | 2-fold |

**S_q 1st excited degeneracy = q-1 exactly** (full S_q permutation merges all spin fields). Hybrid always has 2-fold (Z_q conjugate pairs). Sharpest discriminator between models.

**R_ε (energy field) differs 2-3x between models.** At q=7: S_q R_ε=3.48, hybrid R_ε=8.07. Different operator content = different universality classes.

**Δ·N drift comparison (Sprint 077c):** S_q Δ·N ≈ 0.6 nearly q-independent (q=5-10). Hybrid Δ·N decreases: 0.48 (q=5) → 0.36 (q=7) → 0.25 (q=10).

**Genuinely novel (for our hybrid model):** g_c(q) formula, c·x₁ ≈ 0.112, systematic CFT data q=2-30, distinct universality class from both S_q Potts AND Z_q clock (Sprints 065, 076). A new universality class with continuous power-law transitions (finite ν) that is neither first-order nor BKT.

**Retracted claims:** ~~"No CFT predictions exist for q>4"~~ (complex CFTs do). ~~"Analytic continuation is WRONG"~~ (gives complex values, verified). ~~"1D quantum Potts NEVER first-order"~~ (true for our Z_q-field model, but S_q Potts IS first-order for q>4).

**Resolved questions (Sprint 065):** ~~(1) Does hybrid flow to clock universality at large n?~~ **NO** — different ν (power-law vs BKT). ~~(3) One transition or two (like clock's two BKT)?~~ Clock has BKT; hybrid has single power-law transition.

**Resolved (Sprint 067): No floating phase in hybrid model.** Wide gap scan (g=0.02-3.0), DMRG entropy (n=16,24), and clock comparison all confirm: hybrid has ONE transition. Clock q=5 has floating phase spanning g=[0.30, 0.92] (Δg=0.62); hybrid is critical only at g_c±0.06. The δ-function coupling suppresses intermediate phase. Fourth difference between hybrid and clock universality classes.

**Resolved (Sprint 066):** (2) ~~Truly continuous or weakly first-order?~~ **No first-order signal detected.** Three diagnostics: Δ·N decelerates at q=10 (n=4,5,6), gap minimum is size-independent (~0.030), anomalous FSS saturates at ~5.5% for q≥7. Consistent with continuous transition with sign-flipped corrections for q≥5. Cannot rule out ξ >> 12 weakly first-order, but no positive evidence.

**In future sprints: call this "Potts-clock hybrid" or "Z_q Potts chain", not "quantum Potts."**

**2D status (Sprints 068-071):** L=2 is out of scaling regime for all q. Full 2D (torus) exhausted at q=5: L=3 only usable size. Sprint 071 extended to Ly=2 cylinder (ladder) geometry: g_c(cyl, q=2)=0.451, g_c(cyl, q=5)=0.714. Order parameter smooth for both — no first-order signal on cylinder either. DMRG impractical for q≥5 cylinder (d=5 per site → chi≥50 needed, too slow). Exact diag gap×Lx crossing on cylinder reliable for Lx=3-7 (q=2) and Lx=3-4 (q=5).

## 2D Extension (Sprint 068)

**First study of the hybrid model on 2D square lattices.** H = -Σ_{<ij>} δ(s_i,s_j) - g Σ_i (X+X†) on periodic L×L torus. Method: gap×L crossing (z=1 at Lorentz-invariant critical point).

| q | g_c (2D) | Sizes | gap×L at g_c | Verdict |
|---|----------|-------|--------------|---------|
| 2 | 0.771 | L=2,3,4 | 2.363 | **Continuous confirmed** (L=3,4 collapse to 4 decimals) |
| 3 | 1.267 | L=2,3 | 6.10 | Consistent with continuous |
| 5 | 1.588 | L=2,3 | 3.50 | No first-order signal (inconclusive, only 2 sizes) |

**2D/1D g_c ratio ≈ 3-4 for all q.** Non-monotonic: 3.08 (q=2), 3.80 (q=3), 3.60 (q=5). Consistent with doubling coordination number z=2→4.

**Z_q conjugate pair degeneracy (gap₂=gap₁) is EXACT in 2D** for all q, all L, all g tested. The symmetry structure carries over from 1D.

**L=2 (2×2) is out of scaling regime for q=2.** Gap×L = 4.19 at the L=3,4 g_c, far from converged 2.36. L≥3 needed for 2D Ising. This FSS correction likely worsens for larger q.

**Key limitation:** L=2,3 cannot distinguish continuous from first-order at q=5. The gap×L crossing ratio is trivially 2/3 with only 2 sizes. dE/dg is smooth but L=3 may be too small. Need QMC or tensor networks for definitive answer.

## Cylinder Geometry (Sprints 071-075)

**Cylinder geometry:** Open x-direction, periodic y-direction. Enables gap×Lx crossing with exact diag for moderate sizes. Bridges 1D and 2D.

**Ly=2 cylinder (z=3) — complete dataset:**

| q | g_c(1D) | g_c(cyl, Ly=2) | g_c(2D) | cyl/1D | Crossing pairs |
|---|---------|----------------|---------|--------|----------------|
| 2 | 0.250 | 0.451 | 0.771 | 1.80 | 3 |
| 3 | 0.333 | 0.565 | 1.267 | 1.69 | 3 |
| 5 | 0.441 | 0.714 | 1.588 | 1.62 | 1 |

**Cyl/1D ratio monotonically decreasing in q:** 1.80→1.69→1.62.

**Ly convergence for q=2 (complete to Ly=5, Sprint 075):**

| Ly | g_c | Progress to 2D | Increment |
|----|-----|---------------|-----------|
| 1 | 0.250 | 0% | — |
| 2 | 0.451 | 38.6% | +0.201 |
| 3 | 0.655 | 77.7% | +0.204 |
| 4 | 0.688 | 84.0% | +0.033 |
| 5 | 0.701 | 86.6% | +0.013 |
| ∞ | 0.771 | 100% | — |

**Convergence is POWER-LAW ~1/Ly², NOT exponential (Sprint 075c).** Best fit: g_c(Ly) = 0.771 - 1.285/Ly^2.03 (RMS=0.016). Exponential fit RMS=0.025. The Ly=1→3 regime is nearly linear (+0.201, +0.204), with abrupt deceleration at Ly≥4 (+0.033, +0.013). Exponential model from Sprint 073 predicted Ly=5: 0.745, actual: 0.701 (5.9% error).

**Ly convergence for q=3 (to Ly=3, Ly=4 infeasible):**

| Ly | g_c | Progress to 2D |
|----|-----|---------------|
| 1 | 0.333 | 0% |
| 2 | 0.565 | 24.8% |
| 3 | 0.797 | 49.7% |
| ∞ | 1.267 | 100% |

**Ly convergence for q=5 (to Ly=3, Sprint 074):**

| Ly | g_c | Progress to 2D |
|----|-----|---------------|
| 1 | 0.441 | 0% |
| 2 | 0.714 | 23.8% |
| 3 | 0.974 | 46.5% |
| ∞ | 1.588 | 100% |

Note: q=5 Ly=3 from (Lx=2,3) crossing only — Lx=2 has known FSS issues. Lower confidence than q=2,3 values.

**Convergence is q-dependent, saturating for q≥3 (Sprints 073-074).** At Ly=3: q=2 at 77.7%, q=3 at 49.7%, q=5 at 46.5%. The big slowdown is q=2→q=3.

**Cylinder entanglement entropy (Sprint 075b).** c_eff from S ~ (c/6)·ln(Lx) at g_c: Ly=1 (0.61), Ly=2 (0.75), Ly=3 (0.82). c_eff growth is FSS overshoot — 2D Ising also has c=0.5. S per bond cut DECREASING with Ly: 0.285→0.157→0.107→0.089 — area-law behavior emerging.

**Computational limits:** q=2 Ly=5 Lx=5 (dim=33M) takes 541s/pt. q=3 Ly=4 Lx=4 (dim=43M) 838s/pt. Cylinder convergence with exact diag capped at Ly=5 for q=2, Ly=3 for q≥3.

**All transitions smooth on cylinders.** Order parameter continuous for q=2,3,5 on Ly=2 cylinder. Max slope: q=3 (1.99) > q=5 (1.88) > q=2 (1.21).

**DMRG impractical for q≥5 cylinder.** d=5 per site → chi=20 has massive truncation errors.

**POTENTIALLY NOVEL:** If the hybrid remains continuous in 2D at q>4, it would contradict the standard Potts (first-order) and clock (BKT) behavior. The 1/Ly² convergence law and q-dependent convergence rate may be novel characterizations of the 1D→2D crossover.

## 2D Entanglement Entropy (Sprint 069)

**Ordered phase entropy = ln(q) EXACTLY** for all q=2,3,5 and all L=2,3,4 tested. The Z_q-symmetric ground state in the ordered phase is a GHZ-like superposition |ψ⟩ = (1/√q) Σ_s |s,s,...,s⟩, giving S = ln(q) for any non-trivial bipartition. Size-independent, geometry-independent.

**Critical entropy follows area law.** S/boundary at w=1 strip, L=3: q=2 (0.065), q=3 (0.039), q=5 (0.053). Non-monotonic in q. L=2 is out of scaling for all q.

**Entropy peak is NOT at g_c for accessible sizes.** At L≤4, the GHZ entropy ln(q) dominates over the critical area-law entropy α·2L. Crossover at L ~ ln(q)/(2α) ≈ 10-16. This is qualitatively different from 1D, where S peaks at g_c for modest n.

**Sharp entropy drop at g ≈ 0.6·g_c** (not at g_c itself). For q=5 L=3: S drops from 1.62 to 0.33 between g=0.68 and g=1.57. This marks destruction of the ordered-phase superposition, not the critical point.

**Complement symmetry:** S(w) = S(L-w) exact on torus for all q, L, g.

## 2D Transition Diagnostics (Sprint 070)

**q=2 confirmed continuous in 2D** by three independent probes (L=3→4):
- d²E₀/dg² peak scales as L^0.16 (consistent with α=0, log divergence for 2D Ising)
- χ_F/N (fidelity susceptibility from overlap) scales as L^0.94 (consistent with ν=1)
- dE₀/dg smooth across g_c with 1.2% change L=3→4 (no latent heat)

**q=5 INCONCLUSIVE — accessible sizes too small.** Only L=2,3 available; L=2 is out of scaling (10-50x off in all quantities). Cannot extract meaningful L-scaling from single reliable size (L=3). Key observations:
- dE₀/dg smooth across g_c (no latent heat at L=3)
- F_min = 0.985 at L=3 (lower than q=2: 0.999, q=3: 0.998 at same L)
- d²E₀/dg² peak at L=3 is SMALLER for q=5 (1.58) than q=2 (3.00) — less singular
- Peak positions are far below g_c for all q at accessible L

**L=2 pathologically out of scaling for ALL diagnostics.** Not just gap (Sprint 068) — d²E/dg² peaks are 10-50x smaller, χ_F is 20x smaller, overlap is 1.000. L=2→3 scaling exponents (5-7) are artifacts.

**Eigenstate-sum χ_F fails at dim > 5000.** Truncation to k=10 states out of millions severely underestimates χ_F. Ground state overlap method is the correct approach for large systems.

**Bottleneck: q=5 L=4 has dim ≈ 10^11.** Exact diag impossible. Need cylinder DMRG (Ly=2) or QMC to resolve the 2D transition question.

## MI-CV as Phase Transition Order Parameter (Sprints 030, 036, 037)
MI-CV classifies transition type: second-order (crossing), first-order (step), BKT (dome). Confirmed at n=8-50. Potts crossings vindicated at true g_c. DO NOT use MI-CV for ν extraction (Sprint 053).

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

**Energy gap method (Sprint 051).** Δ·N scale-invariant at g_c. FSS correction: n=4,6 pairs +4.8%, n=6,8 pairs +2.5% (calibrated from q=2,3). Crossings approach from below.

**Hybrid has continuous transitions for all tested q** (q=5,10,20). Power-law second-order with finite ν (Sprint 065). Sprint 066 found NO first-order signal at q=10 up to n=6: gap minimum size-independent, Δ·N decelerating. Clock model has BKT (ν→∞). S_q Potts is first-order for q>4.

**Clock ≠ Hybrid for q≥4 — CONFIRMED (Sprints 041-042, 063, 065, 067).** Four independent differences at q=5:
1. ν: hybrid 0.83 (finite) vs clock →∞ (BKT) [Sprint 065]
2. c/x₁: hybrid 10.77 vs clock 9.43 (12%) [Sprint 065]
3. g_c: hybrid 0.441 vs clock 0.52 (18%) [Sprint 063]
4. **Floating phase: clock has it (g=[0.30,0.92]), hybrid does NOT** [Sprint 067c]
NOT FSS artifacts — qualitatively different transition types and phase structures.

**ν(q) extraction (Sprint 053, 065).** Energy gap slope: d(Δ·N)/dg ~ N^{1/ν}·(1+b/N), b=0.86 from q=3 calibration.

| q | g_c | ν (corrected) | Confidence | 2D classical exact |
|---|-----|---------------|------------|-------------------|
| 2 | 0.250 | 1.00 | High (3 pairs) | 1.0 |
| 3 | 0.333 | 0.86 | High (4 sizes) | 5/6 = 0.833 |
| 4 | 0.392 | 0.82 | High (3 pairs, all agree) | 2/3 + logs |
| 5 | 0.441 | 0.85 | Medium (1 pair) | — (1st order in 2D) |
| 7 | 0.535 | 0.97 | Low (1 pair) | — |
| 10 | 0.684 | 1.12 | Low (1 pair (4,5)) | — |

**ν(q=3,4,5) ≈ 0.82-0.86, nearly constant (hybrid).** Sprint 065 confirmed ν≈0.83 at q=5 with n=4,6,8.

**Method ranking for ν extraction:**
1. Corrected power-law (3% error) — BEST with ≥3 sizes
2. 1/N extrapolation of pairwise ν (3%)
3. Direct power-law (15%), Data collapse (DO NOT USE at n≤10), MI-CV collapse (DO NOT USE)

**DMRG excited states fail for Potts** (orthogonal_to gives gap=0). MI-CV requires χ > d².

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

**Coulomb gas analytic continuation gives COMPLEX c for q>4** (see Model Identity section). Gorbenko-Rychkov-Zan predict Re(c) > 1, consistent with our measurements. The naive real-valued continuation (g>1) is wrong, but the complex continuation is correct. Ma & He (2019) already measured effective c(q=5,6,7) on the S_q Hermitian chain.

**Quadratic interpolation through q=2,3,4 also WRONG (Sprint 056).** Gives c(5)=1.10 (accidentally exact!) but c(7)=1.00 and c(10)=0.10.

**c(q=4) has anomalous FSS.** Both methods show flat, non-converging overshoot (+14-23%). Consistent with logarithmic corrections at marginal q=4 (Ashkin-Teller).

**c(q≥5) for our Potts-clock hybrid model.** c>1 is expected (both standard Potts Re(c) and clock model c exceed 1). Our specific c values and the c ~ ln(q) fit formula may be novel FOR THIS MODEL, but c>1 itself is not surprising. Compare to Ma & He's S_q values and clock model values to establish what's model-specific.

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

**NOT a free boson at finite q, but converges as O(q^{-3}) (Sprint 061).** Free boson predicts R(σᵏ)/R(σ) = k². Measured ratios always BELOW k², approaching with corrections that scale as q^{-3.0} (k=2), q^{-3.2} (k=3,4). At q=30: σ² at 99.85% of k²=4. Requires Z_q charge resolution to separate σ² (charge 2) from ε (charge 0).

**Energy field R_ε grows with q at n=4.** Minimum at q=3 (R=6). For q≥10 at n=4: R_eps ≈ 0.88·q (linear growth). Sigma² and epsilon are ALWAYS in different charge sectors and DIVERGE, not merge. At larger n (q=10, n=6), R_eps=8.3 — may saturate in thermodynamic limit.

**Physical interpretation:** c(q)~ln(q) growth from proliferation of (q-1) spin field primaries. Free boson is the q→∞ limit on a decompactifying circle: R ~ √(ln q) explains c ~ ln(q) and x₁ → 0.

**Operator content is novel FOR OUR HYBRID MODEL.** Tang et al. (2024) measured operator content for S_q Potts (non-Hermitian). Our data is for the Potts-clock hybrid — a different model. Compare to check if universality classes match.

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

**x₁ peaks at q=3 (2/15 ≈ 0.133).** Decreasing for q>3, NOT saturating. From n=4 descendant gap: x₁(15)≈0.071, x₁(20)≈0.055, x₁(30)≈0.038. Approaches 0 as q→∞ ~ q^{-0.3} (best-n fit) or faster.

**q=4 has logarithmic FSS corrections.** The excess (c/x₁ - 8)·N² grows as ~N², confirming sub-power-law (logarithmic) convergence. Consistent with marginal Ashkin-Teller point. Cannot determine if c/x₁(∞) = 8.0 or ~8.5 from n≤10.

**Anomalous FSS for q≥5 (Sprint 058, 066).** Δ₁·N INCREASES with N (sign-flipped correction vs q≤4). At q=10 (n=4,5,6): increases decelerate (+3.1%, +2.2%), consistent with convergence from below. The anomaly saturates at ~5.5% (n=4→6) for q≥7. NOT a first-order signal — gap minimum stays ~0.030 (size-independent). Coincides with c>1 and novel CFT regime.

**x₁(q) best fit (Sprint 062a):** x₁ ≈ 0.206·q^(-0.449), RMS=0.010. Neither 1/ln(q) nor power of ln(q-1) fits well.

**c·x₁ ≈ 0.112 for Potts q≥3 (Sprints 062-063).** Product c(q)·x₁(q) = 0.112 ± 0.005 for q=3-15 Potts. NOT exactly 1/9 (q=3 exact is 8/75=0.1067). q=2 Ising is an outlier (c·x₁=1/16). Potts-specific: clock q=5 gives c·x₁≈0.146 (30% higher). The near-constancy is approximate, useful for predictions but not an exact identity.

**c(q=15) = 1.549 (DMRG n=8, Sprint 063a).** Gives c·x₁ = 0.110, supporting c·x₁ constancy up to q=15.

**Single compact boson RULED OUT (Sprint 062a).** c>1 for q≥5 is incompatible with any single free boson (always c=1). Multiple effective DOF required.

**Novel for hybrid model:** x₁(q), c·x₁ ≈ 0.112, and the full CFT characterization are new for the Potts-clock hybrid Hamiltonian. Tang et al. measured x for S_q Potts at q=5 (complex-valued). Our real x₁ values are for a different model. The key question is whether the hybrid universality class is truly distinct — Sprint 063 suggests yes (different c·x₁ from clock).

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

**C_sse peaks at q=3 (~0.54), then decreases monotonically.** Power law: C_sse ~ 1.5·q^{-0.87} (best-n, q≤10) or ~2.3·q^{-1.19} (n=4, q≤30). Extended to q=30 at n=4 via charge resolution: C_sse(30)≈0.038. Approaches 0 as q→∞ (epsilon decouples).

**Epsilon identification pitfall for large q:** At q≥7, sigma harmonics (sigma^2, sigma^3, ...) densely populate R < R_epsilon. Naive charge detection fails because degenerate pairs (charges k, q-k) in superposition average to misleading charge values. Must use R_epsilon from independent spectrum measurements to identify epsilon correctly.

**Novel for hybrid model:** First C_sse measurement for the Potts-clock hybrid at q=5,7,10. Tang et al. (2024) measured OPE for non-Hermitian S_q Potts at q=5 (complex values). Our real-valued C_sse for a different model is new.

## Twisted BC & Spin Stiffness (Sprint 062)

q=2,3 are Luttinger liquids (ρ_s·L converges); q≥4 are NOT (10-22% drift). Twist field = σ field exactly at q=2. Twist ratios converge toward spectrum ratios from below.

## Conformal Tower Structure (Sprint 059)

**Genuine CFT confirmed for ALL q=2-10.** Four tests: descendant degeneracy exact, momentum = spin ±1, gap converges from below, tower organization universal. Anomalous FSS at q≥5 (15-25% gap error) — qualitative structure perfect despite quantitative offset. Tower density → free boson continuum as q→∞.

**Entropy profile overshoot grows with q at fixed n:** At n=16: 8.7% (q=3), 14.4% (q=4), ~20% (q=5). At n=8: ~25% (q=7). Large-q profile method requires larger n — q=7+ is computationally infeasible at n≥12 with chi=20.

**Entropy FSS does NOT give ν.** S ~ (c/6)ln(ξ), not a power law. Standard FSS collapse fails (gives ν=3 for TFIM). Use energy gap slope for ν.

**Correlation length from correlator decay.** Saturates at ξ ~ n/4 near criticality. iDMRG correlation length also saturates at criticality (Sprint 055). Not useful for ν extraction near g_c.

**Locating g_c via entropy.** The critical point can be identified as where: (1) dS/dg peaks (pseudo-critical point, drifts toward g_c with n), (2) S grows logarithmically with n (ordered/disordered → S saturates), (3) central charge c from S(n) matches expected value.

**Technique: All-pairs MI.** Two methods:
1. Gell-Mann correlation reconstruction: **ONLY reliable for d≤3** (validated: diff=0 at n=8 for d=2,3). For d=4, small systematic errors create artificial MI-CV crossings (Sprint 046). For d=5, errors up to 11x (Sprint 045). DO NOT USE for d≥4.
2. **Direct MPS tensor contraction (Sprint 043):** computes ρ_ij by contracting MPS tensors directly. O(n·χ²·d) per pair. 1.4s for 28 pairs at d=10, n=8. Works for ANY d. **The only reliable method for d≥4.**

## Archetype Boundaries ≠ Phase Boundaries
I3 sign change occurs at Δ≈0.7 in XXZ, inside the XY phase — not at either thermodynamic transition (Δ=-1 or Δ=1). The entanglement phase diagram has its own topology distinct from thermodynamics.

## Bisognano-Wichmann Locality (Sprints 032-034, 064)
BW: H_E ≈ H_phys × β(x), Unruh-like gradient. Peaks in ordered phase. Controlled by H/G-inv ratio.
At n_A=4: q=2 (91%), q=3 (76.5%), q=4 (56%), q=5 (42.3%). ~15% drop per unit q.
Z_q-invariant fraction = exactly 1/q. H/G-inv ~ 5/q^5. Peak shifts to ordered phase for q≥4.

## QEC Arc — Closed (Sprints 014-028)
[[5,1,3]] is basis-isotropic (confirmed on IBM hardware Sprint 025). Active QEC never beats passive at small scale.

## Hardware (Sprint 025)
ibm_kingston, 20s QPU. [[5,1,3]] asymmetry 0.040 vs 3-qubit 0.254. 580s QPU remaining.

## Phase Diagram Trajectories
- **TFIM**: Democratic(GHZ) → Scale-Free → Product (one-way)
- **XXZ**: Democratic → Scale-Free → Geometric → Scale-Free → Democratic (loop)
- **Hybrid q=3**: GHZ-3 → Product (like TFIM but with triplet spectrum from Z₃)
- **Hybrid q=5**: GHZ-5 → Product (continuous transition at g_c≈0.44)
- **Hybrid q=10,20**: Continuous transitions confirmed
