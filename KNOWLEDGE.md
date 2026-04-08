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


**Genuinely novel (confirmed):** Casimir tracks Re(c) (Sprint 098), chi_F spectral mechanism (Sprint 107), alpha(q) scaling (Sprints 103-116, 10 pts q=5-25). Walking-specific, not universal (Sprint 104-105).

**Retracted:** ~~"No CFT predictions for q>4"~~, ~~"Analytic continuation WRONG"~~, ~~"Potts NEVER first-order"~~ — see archive for details.

## Detailed Findings

All detailed findings from sprints 029-115 (walking regime, entanglement spectrum, BW analysis, 2D extension, hybrid model characterization, MI-CV, entropy methods, CFT content, OPE, SREE, χ_F) are in **KNOWLEDGE_ARCHIVE.md**. Read that file when you need:
- Methodology details (how to extract c, x₁, ν, C_sse)
- Past numerical values or data tables
- What's been tried and ruled out
- Hybrid model vs clock vs S_q Potts comparisons
