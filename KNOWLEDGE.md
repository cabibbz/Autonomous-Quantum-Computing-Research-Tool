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

### χ_F Spectral Decomposition: Universal Mechanism, Model-Specific Exponents (Sprints 119-120)
Single-multiplet dominance (frac=1.000) and spectral gap selection rule are **universal** — identical in both S_q Potts and hybrid Potts-clock. The decomposition alpha = beta_me + 2*z_m - 1 is exact in both. But the exponents diverge:

| q | Hybrid alpha | S_q alpha | Hybrid z_m | S_q z_m |
|---|-------------|-----------|-----------|---------|
| 3 | 1.47 | 1.40 | 1.03 | 1.03 |
| **4** | **1.55** | **1.77** | **1.00** | **~1.05** |
| 5 | 1.41 | 2.09 | 0.91 | 1.31 |
| 7 | 0.95 | 2.65 | 0.77 | ~1.4 |
| 10 | 0.55 | ~3.2 | 0.69 | — |

**z_m crosses 1 at q_cross = 3.58 ± 0.04 for hybrid** (Sprint 121a, 5 fit forms all converge). z_m > 1 → walking/super-scaling, z_m < 1 → continuous/ordinary. Hybrid q=4 is barely walking (z_m=1.004), not exactly marginal. Alpha peak also at q≈3.58. **S_q q=4: alpha=1.777±0.003, z_m=1.092±0.001** (Sprint 121b, 8 sizes n=4-11) — firmly walking, z_m slowly drifting toward 1 (1.097→1.079). 19% alpha divergence between models at largest sizes.

**g_c(hybrid):** q=2→0.250, q=3→0.333, q=4→0.393, q=5→0.438, q=7→0.535, q=10→0.684 (Sprint 120a).

**Literature on S_q q=4:** Marginal operator gives log corrections: chi_F ~ L^2(ln L)^{-p} (Salas-Sokal 1997, Balog et al. 2007). Predicted p=3/2 (specific heat sector) but NO prior chi_F log correction measurement exists. At n=4-11, alpha=2 with log corrections is the WORST fit (dAIC=42 vs power+1/N^2). Best fit: alpha=1.757 with 1/N^2 correction. Asymptotic regime needs L>>11 (Lv et al. 2019 needed L=1024). Hybrid q=4 alpha=1.460 with 1/N^2 correction — clearly different universality class.

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
