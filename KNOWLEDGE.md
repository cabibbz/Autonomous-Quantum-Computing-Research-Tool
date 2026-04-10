# Accumulated Knowledge — Edit by topic, not by sprint

## Open Items — Check each sprint, remove when done
- **~~Test energy-entropy hierarchy in a DIFFERENT model~~** DONE (Sprint 104). Tested in J1-J2 chain. Result: **NOT universal.** Hierarchy direction is model-dependent. Potts walking gives O(1) entropy deviation — unique to walking mechanism. J1-J2 gives only O(1%) differences. The Casimir-Re(c) finding is walking-specific, not a general CFT principle.
- **~~Harden chi_F scaling at q=5~~** DONE (Sprint 103, updated Sprint 116). alpha(q) mapped q=2-25. All polynomial/power-law forms ruled out. **Log+loglog alpha(q) = 2.62 ln(q) - 1.77 ln(ln(q)) - 1.26 marginally AIC-best** (dAIC=1.4 over pure log) with 10 data points (Sprint 116). Pure log: 1.86 ln(q) - 0.96.
- **~~Test χ_F in different model~~** DONE (Sprint 105). J1-J2 BKT: invisible (α→0). MG first-order: saturates. Walking super-scaling is unique.
- **~~Understand χ_F mechanism~~** DONE (Sprint 106). α = β_me + 2z_m - 1. Single multiplet dominates. Spectral gap symmetry-forbidden.
- **~~Harden χ_F mechanism~~** DONE (Sprint 107). 5 sizes at q=5, cross-validated. CONFIRMED NOVEL.
- **~~Log corrections to α(q)~~** CORRECTED (Sprint 109). Sprint 108's 1/ln(N) extrapolation was WRONG for q=2,3. Power-law 1/N² corrections recover exact ν: q=2→α_∞=1.00 (exact 1.0), q=3→α_∞=1.40 (exact 7/5). Walking (q≥5) zero corrections confirmed. q=4 BKT genuinely slow (neither model converges at accessible sizes).
- **KNOWLEDGE.md is over budget (~560 lines vs ~200 target).** Compress old sections into one-line summaries.
- **Hardware validation** — 580s QPU unused for 89 sprints. Strongest prediction: q=2 Ising χ_F at g_c, or Heisenberg chain c_eff on 5-10 qubits.

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

## ⚠ Model Identity — READ THIS (Audit April 2026)

### TWO MODELS WERE STUDIED — NOT ONE

**Sprints 033-075:** Studied the **Potts-clock hybrid** H = -Jδ(s_i,s_j) - g(X+X†). Z_q symmetry.
**Sprints 076-118:** Switched to the **standard S_q Potts model** H = -Jδ(s_i,s_j) - gΣ_{k=1}^{q-1} X^k. S_q symmetry. All six claimed novel findings (098, 102-118) are on this model.

The switch happened at Sprint 076 (introduced S_q for comparison) and was never reversed. Code audit confirms: all experiment scripts from 076+ use `build_sq_potts` with `for k in range(1, q)` or TeNPy `SqField = ones(q,q) - eye(q)`.

| Model | Coupling | Transverse field | Symmetry | Sprints | q>4 transition |
|-------|----------|-----------------|----------|---------|---------------|
| S_q Potts | δ(s_i,s_j) | Σ_{k=1}^{q-1} X^k | S_q | **076-118** | Weakly 1st-order (walking) |
| Z_q clock | cos(2π(s_i-s_j)/q) | X + X† | Z_q | 034,065,067 | Two BKT |
| Potts-clock hybrid | δ(s_i,s_j) | X + X† | Z_q | **033-075** | Continuous 2nd order |

For q=2,3 all three are equivalent (up to field rescaling at q=2). For q≥4 they differ.

### Implications for novelty claims

The S_q Potts model is the SAME model studied by Gorbenko-Rychkov-Zan, Ma & He, Tang et al., and Jacobsen & Wiese. Therefore:
- Agreement between Casimir energy and Re(c) is **expected** (GRZ predictions are for this model)
- g_c = 1/q is a **known exact result** (Kramers-Wannier self-duality)
- The model is NOT novel — but our **probes** (χ_F spectral decomposition, Casimir vs entropy comparison) may be

### What IS still potentially novel on the S_q model
- **χ_F spectral decomposition** (Sprint 107): selection rule + single-multiplet dominance mechanism
- **χ_F scaling data across q** (Sprints 102-117): systematic measurement, appears to be new data
- **Casimir vs entropy systematic comparison** (Sprint 098): useful quantitative observation

### ⚠ Spectral chi_F: Methodology Status (Sprints 125-126)
Spectral Lehmann sum was missing factor 2 (Sprint 125). Even corrected, spectral method captures only the dominant state — non-dominant states at high eigenvalue indices cause **systematic negative alpha bias of 0.005–0.038** (Sprint 126). Bias grows with system size. **Use exact chi_F (finite-difference) for all exponent claims.** The selection rule is standard Z_q representation theory. The formula alpha = beta_me + 2·z_m − 1 is a tautological identity.

### chi_F Exponents: Authoritative Values (Sprints 127-128, exact chi_F, GPU-extended sizes)

| q | Hybrid alpha | S_q alpha | # sizes (S_q) | Pairwise drift (S_q) |
|---|-------------|-----------|---------------|---------------------|
| 3 | 1.481+/-0.014 | 1.468+/-0.012 | 6 (n=4-14) | Decreasing (1.57->1.44) |
| **4** | **1.549+/-0.012** | **1.794+/-0.011** | 6 (n=4-11) | Oscillating ~1.79 |
| 5 | 1.352+/-0.043 | **2.139+/-0.019** | 5 (n=4-10) | Increasing |
| **6** | **1.186+/-0.038** | **2.375+/-0.006** | 6 (n=4-9) | Increasing (2.35->2.40) |
| 7 | 0.971+/-0.058 | 2.584+/-0.015 | 4 (n=4-8) | Noisy |

Hybrid alpha(q) non-monotonic: peaks at q~4 (1.55), drops to ~1.0 at q=7. Extrapolation gives alpha_inf=0 for q>=5 (sub-power-law, possibly logarithmic). S_q alpha(q) monotonically increasing, pairwise trending upward for q>=5.

**z_m crosses 1 at q_cross = 3.58 +/- 0.04 for hybrid** (Sprint 121a). z_m > 1 = walking, z_m < 1 = continuous. **S_q q=4 pairwise alpha oscillates around 1.79** -- no clear drift direction at accessible sizes.

**alpha(q) for S_q:** alpha(q) = 1.337*ln(q) - 0.023, chi2/dof=5.0 (Sprint 128, 5 data points q=3-7). Quadratic-in-ln(q) marginally better (chi2/dof=3.7). q=6 sits exactly on the log curve.

**g_c(hybrid):** q=2->0.250, q=3->0.333, q=4->0.393, q=5->0.438, q=6->0.474, q=7->0.535, q=10->0.684.

**Asymptotic extrapolation (Sprint 128):** alpha_eff(N) = alpha_inf + c/N^p. Validated at q=3 (alpha_inf=1.412 vs exact 1.400). S_q q=4 too noisy. S_q q>=5 alpha still growing with N. Hybrid q>=5 collapses to alpha_inf=0.

**Literature on S_q q=4:** Marginal operator gives log corrections: chi_F ~ L^2(ln L)^{-p} (Salas-Sokal 1997, Balog et al. 2007). Predicted p=3/2 (specific heat sector) but NO prior chi_F log correction measurement exists. At n=4-11 (periodic BC), alpha=2 with log corrections is the WORST fit (dAIC=42 vs power+1/N^2). Best fit: alpha=1.757 with 1/N^2 correction. Hybrid q=4 alpha=1.460 with 1/N^2 correction -- clearly different universality class.

**DMRG extension (Sprint 124):** Open-BC chi_F at n=6-20 (8 sizes). Pairwise alpha drifts UPWARD: 1.505->1.523. Power+1/N^2 corrected alpha_open=1.524+/-0.002. Log-corrected alpha=2 still worst fit (R^2=0.9997 vs 0.999999). **But drift direction is upward** -- consistent with eventual convergence to higher value (periodic 1.77 or log-corrected 2.0). Asymptotic regime needs L>>20 (likely L>100). iDMRG overlap method FAILED for S_q Potts (non-abelian symmetry).

### Hybrid model findings (Sprints 033-075, 119-121) — separate body of work
Sprint 065 confirmed hybrid ≠ clock, Sprint 076 confirmed hybrid ≠ S_q Potts. Sprints 119-121: chi_F spectral decomposition on hybrid confirms continuous transitions for q≥5 and walking→continuous boundary at q_cross=3.58. The hybrid model at q≥4 is genuinely a different universality class from S_q Potts.

**Key literature (search before claiming novelty):**
- **Gorbenko, Rychkov & Zan (JHEP 2018, SciPost 2018):** Complex CFT for q>4 S_q Potts.
- **Ma & He (PRB 99, 195130, 2019):** Measured c_eff for q=5,6,7 on Hermitian S_q chain.
- **Tang et al. (PRL 133, 076504, 2024):** Non-Hermitian q=5 S_q Potts.
- **Jacobsen & Wiese (PRL 133, 077101, 2024):** All S_q Potts exponents via analytic continuation.
- **Sun, Luo & Chen (arXiv:2006.11361, 2020):** Z_q clock has two BKT transitions for q>4.

**Retracted:** ~~"No CFT predictions for q>4"~~, ~~"Analytic continuation WRONG"~~, ~~"Potts NEVER first-order"~~, ~~"Our model is a novel Potts-clock hybrid" (for sprints 076+)~~ — see archive for details.

## Detailed Findings

All detailed findings from sprints 029-115 (walking regime, entanglement spectrum, BW analysis, 2D extension, hybrid model characterization, MI-CV, entropy methods, CFT content, OPE, SREE, χ_F) are in **KNOWLEDGE_ARCHIVE.md**. Read that file when you need:
- Methodology details (how to extract c, x₁, ν, C_sse)
- Past numerical values or data tables
- What's been tried and ruled out
- Hybrid model vs clock vs S_q Potts comparisons
