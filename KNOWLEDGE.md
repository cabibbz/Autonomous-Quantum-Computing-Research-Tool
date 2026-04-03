# Accumulated Knowledge — Edit by topic, not by sprint

## Open Items — Check each sprint, remove when done
- **~~Test energy-entropy hierarchy in a DIFFERENT model~~** DONE (Sprint 104). Tested in J1-J2 chain. Result: **NOT universal.** Hierarchy direction is model-dependent. Potts walking gives O(1) entropy deviation — unique to walking mechanism. J1-J2 gives only O(1%) differences. The Casimir-Re(c) finding is walking-specific, not a general CFT principle.
- **~~Harden χ_F scaling at q=5~~** DONE (Sprint 103, updated Sprint 110). α(q) mapped q=2-9. Formula revised: α ≈ 0.26q + 0.81 (walking-only, q=5-9). Old 0.315q+0.469 biased by q=4 BKT inclusion.
- **~~Test χ_F in different model~~** DONE (Sprint 105). J1-J2 BKT: invisible (α→0). MG first-order: saturates. Walking super-scaling is unique.
- **~~Understand χ_F mechanism~~** DONE (Sprint 106). α = β_me + 2z_m - 1. Single multiplet dominates. Spectral gap symmetry-forbidden.
- **~~Harden χ_F mechanism~~** DONE (Sprint 107). 5 sizes at q=5, cross-validated. CONFIRMED NOVEL.
- **~~Log corrections to α(q)~~** CORRECTED (Sprint 109). Sprint 108's 1/ln(N) extrapolation was WRONG for q=2,3. Power-law 1/N² corrections recover exact ν: q=2→α_∞=1.00 (exact 1.0), q=3→α_∞=1.40 (exact 7/5). Walking (q≥5) zero corrections confirmed. q=4 BKT genuinely slow (neither model converges at accessible sizes).
- **KNOWLEDGE.md is over budget (~560 lines vs ~200 target).** Compress old sections into one-line summaries.
- **Hardware validation** — 580s QPU unused for 83 sprints. Strongest prediction: q=2 Ising χ_F at g_c, or Heisenberg chain c_eff on 5-10 qubits.

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

**S_q Potts looks CFT-like at n≤12 for q=5 (Sprints 076-078).** Clean conformal towers, exact degeneracies, Δ·N stable to <5%. DMRG gap×N INCREASES n=8→12 (1.88→2.00, open BC) — no walking breakdown. Walking correlation length ξ* > 12 sites.

**Complex CFT central charge formula (verified Sprint 079c):** √Q = 2·cos(π/p), c = 1 - 6/[p(p-1)]. For q>4, p is complex: α = arccosh(√q/2), p = iπ/α. Re(c) = 1 + 6α²/(π² + α²). Verified: q=2→0.5, q=3→0.8, q=5→1.138 (all exact).

**c_eff(q=5) ≈ 1.15, matches complex CFT Re(c) ≈ 1.138 to ~1% (Sprints 078c, 079c).** q=5 is the unique "sweet spot" where complex CFT quantitatively predicts c_eff at ALL accessible sizes (n=6-24).

**Walking regime boundary mapped (Sprints 079-081).** c_eff/Re(c) ratio decays exponentially for q>5:

| q | Re(c) | c_eff(best) | n_best | c/Re(c) | dc/d(ln n) | Walking? |
|---|-------|-------------|--------|---------|------------|----------|
| 3 | 0.800 | 0.893 | 24 | 1.12 | — | real CFT (+12% FSS) |
| 5 | 1.138 | 1.152 | 12 | 1.01 | +0.014 | YES — sweet spot |
| 6 | 1.253 | 1.115 | 12 | 0.89 | **-0.048** | MARGINAL BREAKING |
| 7 | 1.351 | 1.059 | 12 | 0.78 | -0.091 | BREAKING |
| 8 | 1.438 | 1.062 | 7 | 0.74 | — | BREAKING |
| 9 | 1.515 | 1.012 | 6 | 0.67 | — | NO |
| 10 | 1.584 | 0.946 | 6 | 0.60 | — | NO |

**Fit: c_eff/Re(c) = 1.004 × exp(−0.105(q−5)).** q=5 is the exact walking threshold.

**Key discriminator: drift rate dc/d(ln n) (Sprint 081).** q=5 c_eff RISES (+0.014), q=6 DROPS (-0.048), q=7 drops faster (-0.091). q=6 drift is 2× closer to q=7 than q=5 → marginal breaking, NOT a second walking case. Walking boundary is a CROSSOVER, not a sharp transition. q=6 stability at n=6-8 was misleading (too short a range).

**q=6 walking extrapolated breakdown at n≈38 (Sprint 081c).** c_eff crosses 0.85×Re(c) at n≈38, just outside DMRG reach with chi=30.

**gap×N INCREASES with q** (2.01 at q=6, 2.10 at q=6 n=10, 2.17 at q=8) even as walking breaks down. Breakdown is in entropy, not the gap.

**c_eff converges to ~1.0-1.15 for ALL q≥5** at moderate sizes — the complex CFT q-dependence only manifests at n << ξ*(q).

**POTENTIALLY NOVEL (Sprints 079-080):** Complete c_eff/Re(c) curve q=3-10 with exponential decay law. First measurement at q=6,8,9. Walking boundary at q=5 exact.

**Spin-spin correlator x_σ(q) ≈ 0.13 nearly universal (Sprint 082).** Periodic chain correlator G(r) = C × [chord(r)]^{-2x_σ} fits to 0.04% accuracy for ALL q=2-8. Validated: q=2 gives 0.123 (exact 0.125, 1.5% FSS), q=3 gives 0.132 (exact 0.133, 1.3%).

| q | x_σ(corr) | v(q) | gap×N |
|---|-----------|------|-------|
| 2 | 0.123 | 1.02 | 0.786 |
| 3 | 0.132 | 0.89 | 0.733 |
| 4 | 0.135 | 0.81 | 0.686 |
| 5 | 0.136 | 0.75 | 0.639 |
| 6 | 0.135 | 0.71 | 0.601 |
| 7 | 0.132 | 0.68 | 0.560 |
| 8 | 0.130 | 0.66 | 0.536 |

**Velocity v(q) = gap×N/(2π·x_σ) decreases monotonically with q.** The decreasing gap×N is primarily velocity reduction, NOT x_σ change. x_σ is nearly constant; v(q) carries all q-dependence.

**Walking breakdown is an entropy phenomenon AT FINITE SIZE (Sprints 082-084, corrected Sprint 087).** Four independent probes:
1. **Correlators (x_σ)**: perfectly conformal for ALL q=2-8 (Sprint 082)
2. **Casimir energy (E₀)**: governed by Re(c) for ALL q=2-8, even where c_eff deviates 40% (Sprint 083)
3. **Entanglement entropy (c_eff)**: deviates from Re(c) for q>5 (Sprints 079-081)
4. **Fidelity susceptibility (χ_F)**: scaling exponent α crosses 2.0 at walking boundary, then increases linearly with q (Sprints 102-103, revised Sprint 110). **Walking-specific** (Sprint 105): BKT gives α→0, MG first-order saturates, only walking gives persistent α>2. **Mechanism identified (Sprint 106):** α = β_me + 2z_m - 1, where z_m is the multiplet gap exponent and β_me is the matrix element growth exponent. Both linear in q. **Updated fits (Sprint 110, walking-only q=5-9):** z_m(q) = 0.082q + 0.741, β_me(q) = 0.098q + 0.333, α(q) ≈ 0.26q + 0.81. Single-multiplet dominance (frac=1.000) confirmed through q=9.
5. **Entanglement spectrum multiplet dominance**: M/[(q-1)/q] crosses 1.0 at q≈4 (Sprints 089-090)

**Microscopic mechanism: entropy concentration in (q-1)-fold multiplet (Sprint 084).** Entanglement spectrum at g_c=1/q shows:
- First excited entanglement level has **(q-1)-fold degeneracy** for all q (S_q symmetry)
- Entanglement gap Δξ ≈ 0.47 + 0.78·ln(q), INCREASES with q
- Level 1 absorbs progressively more entropy: 48% (q=2) → 63% (q=5) → 69% (q=8)
- Tail entropy (levels ≥ 2) correlates with c_eff/Re(c), Pearson r = 0.80

| q | deg_1 | %S(lev 0) | %S(lev 1) | %S(tail) | c_eff/Re(c) | Δξ | M/[(q-1)/q] |
|---|-------|-----------|-----------|----------|-------------|------|-------------|
| 2 | 1 | 33.1% | 47.8% | 19.1% | 1.00 | 1.03 | 1.35 |
| 3 | 2 | 27.3% | 55.9% | 16.8% | 1.12 | 1.32 | 1.09 |
| 4 | 3 | 24.2% | 73.3% | 2.5% | 1.00 | 2.54 | **1.00** |
| 5 | 4 | 22.1% | 62.7% | 15.2% | 1.01 | 1.69 | 0.96 |
| 6 | 5 | 21.0% | 65.8% | 13.2% | 0.89 | 1.87 | 0.91 |
| 7 | 6 | 19.8% | 66.8% | 13.5% | 0.78 | 1.98 | 0.93 |
| 8 | 7 | 19.2% | 69.1% | 11.7% | 0.74 | 2.13 | 0.90 |

**Why only entropy sees walking breakdown at small n:** Energy/gap/correlator observables depend on lowest entanglement levels (ground + first excited), which are perfectly conformal. Entropy sums over ALL levels → sensitive to weight redistribution in the tail. Entanglement gap grows with q at fixed n, but DECREASES with n at fixed q.

**Tail weight grows as UNIVERSAL power law w_tail ~ n^2 for ALL critical q (Sprints 088-090).** DMRG midchain Schmidt spectrum, log-log fit for n≥8:
- q=2: b=1.98, q=3: b=2.01, q=4: b=2.02, q=5: b=2.01, q=7: b=2.07 (all ± 0.14)
- **Mean b = 2.018 ± 0.029** — identical for all q=2-7 within error bars.
- ~~q=7 b≈3.0~~ was pre-asymptotic (only 4 pts n=6-12). With n=6-16 (6 pts) and n≥8: b=2.07 (Sprint 089a).
- Entanglement gap: Δξ ~ -0.43·ln(n) universal for q=2-5, -0.50 for q=7.

**Energy-entropy decoupling is TEMPORARY.** Tail weight reaches 10% at n≈130 (q=2), n≈109 (q=3), n≈108 (q=5). Walking regime creates a **hierarchy of crossover scales**: entropy breaks first (n ~ 10-20), then spectral observables (n ~ 100+). The RATE of tail growth is the same for all q — walking-specific effects are in STATIC PARTITION, not dynamics.

**Weight flows multiplet → tail, not ground → tail.** %S(lev0) saturates while %S(lev1) monotonically decreases, feeding the tail. S_lev0_inf decreasing monotonically: 0.333 (q=2), 0.281 (q=3), 0.251 (q=4), 0.230 (q=5), 0.204 (q=7). All fit S_inf + α/n with R² > 0.999.

**Multiplet dominance M = %S(lev1)/[%S(lev0)+%S(lev1)] is a monotonic walking discriminator (Sprints 089-090).** M increases with q: 0.676 (q=2), 0.724 (q=3), 0.752 (q=4), 0.771 (q=5), 0.797 (q=7) at n=16. The ratio M/[(q-1)/q] crosses 1.0 at q≈4 — **confirmed by q=4 measurement (Sprint 090): M/[(q-1)/q] = 1.003 at n=16.** Real CFT has M > (q-1)/q (multiplet over-represents), complex CFT has M < (q-1)/q (multiplet under-represents).

**q_cross converges to q≈4.0 in the thermodynamic limit (Sprint 090c).** At n=8: q_cross=4.49, n=12: 4.24, n=16: 4.07, n=20: 3.97, n=24: 3.92. Extrapolates to q≈4.0, coinciding exactly with the real-to-complex CFT boundary.

**Democracy index (Sprints 089-090):** At q=2, ground state is depleted relative to democratic share. At q=5, multiplet is depleted. q=4 is the crossover point where democracy index ≈ 1.0.

**Entanglement Hamiltonian H_E = -log(ρ_A) across walking boundary (Sprint 091).** BW fidelity measured for S_q Potts at g_c=1/q, periodic chain. At fixed nA=4:

| q | Potts NN locality | non-Potts fraction | BW fidelity | BW alpha |
|---|---|---|---|---|
| 2 | 99.97% | 0.029% | 0.9997 | 2.39 |
| 3 | 99.92% | 0.081% | 0.9994 | 2.40 |
| 4 | 99.00% | 0.955% | 0.9942 | 3.23 |
| 5 | 97.08% | 2.729% | 0.9832 | 3.20 |

Non-Potts fraction grows exponentially ~exp(1.6·q). **Biggest jump at q=3→4 (11.7×)**, coinciding with real-to-complex CFT boundary. NNN and 3-body Potts operators add <0.2% — residual is genuinely non-Potts operators. BW locality dominated by nA (subsystem size), not q: q=2 drops from 99.98% (nA=3) to 62.6% (nA=7). q=5 non-Potts grows 1.7× faster per nA than q=2. BW alpha jumps from ~2.4 (q=2,3) to ~3.2 (q≥4) — another q=4 discontinuity.

**Operator content of BW corrections (Sprint 092).** Clock-shift (Weyl-Heisenberg) basis decomposition of H_E. Non-BW operators classified as D(density/clock), F(field/shift), M(mixed/XZ-type):
- **Dominant non-BW: MM (mixed×mixed) and DFD/DMD (3-body density-field chains).** MM generalizes Pauli YY — represents entanglement between "position" and "momentum" sectors (φ·π correlators). DFD is density-field-density mediated correlation.
- **At nA=3, non-BW is flat across q** (0.005-0.009% for q=2-5). q-dependent growth only manifests at nA≥4.
- **q=2 vs q=3 at nA=4: same operator types, similar weight** (0.008% both). The 11.7× jump at q=3→4 (Sprint 091) is an AMPLITUDE effect, not new operator types.
- **For q=2 (Pauli basis), non-BW grows from 0.005% (nA=3) to 6.3% (nA=6).** At large nA: 3-body ZXZ/ZZZ dominate, then 4-body ZXXZ appears. Body order and range of dominant corrections increase with nA.

**BW entanglement temperature profile (Sprint 093).** BW predicts H_E = Σ β(x) h_x with sin-envelope β(x). Verified:
- **BW envelope shape is correct** for both Potts coupling and field operators (R² > 0.99 at nA≥4).
- **BW Frobenius R² degrades with q.** nA=4: q=2 R²=0.9995, q=3 0.9992, q=4 0.9984, q=5 0.965. nA=3: monotonic from 0.9997 (q=2) to 0.9986 (q=7).
- **Walking amplifies BW deviation with subsystem size.** 1-R² amplification from nA=3→4: 1.8× (q=2), 1.5× (q=3), 2.0× (q=4), **34× (q=5)**. The 34× at q=5 vs 1.8× at q=2 is the clearest walking amplification measured.
- **BW residual is BULK-concentrated for real CFT, UNIFORM for walking.** Boundary enrichment (fraction at boundaries / expected if uniform): q=2 nA=4: 0.42×, q=3: 0.47×, q=5: 0.98×. Real CFT's BW corrections concentrate in the bulk of A; walking generates corrections uniformly.
- **BW scale α increases with q:** 6.5 (q=2) → 10.3 (q=7). Higher local dimension → higher entanglement temperature.

**BW breakdown vs subsystem size nA — threshold behavior (Sprints 094-095).** Mapped 1-R²(nA) for q=2 (nA=3-7), q=3 (3-6), q=4 (3-5), q=5 (3-4):
- **NOT a power law.** Data shows two-regime threshold: slow growth for nA ≤ nA*, then catastrophic collapse.
- **Real threshold at nA*=6 for q=2 (corrected Sprint 095).** Sprint 094's nA*=5 was inflated 18× by equal-bipartition penalty. At n≫nA, nA=5 gives 1-R²=5.6e-4 (same order as nA=4). nA=6 collapses at ALL n tested (1-R²≈0.17 even at n=18).
- **Equal-bipartition anomaly (Sprint 095).** At nA/n=0.5, BW accuracy is 2-18× worse than at nA/n<0.4. Penalty peaks at threshold nA: nA=5 q=2 has 18× penalty, nA=3 has 2×. This is a systematic bias in all Sprint 094 data.
- **BW corrections are UV (lattice-scale), not IR.** 1-R² depends on nA alone, not nA/n. Increasing chain length n at fixed nA does NOT improve BW. This means corrections come from lattice UV structure, not finite-size effects.
- **B(q) = 0.48q + 1.09 still valid** for exponential rate within each regime. Walking amplification at nA=3 (q=5/q=2): 10-15×, ratio-independent.

**BW threshold mechanism (Sprint 096).** The threshold is NOT in the entanglement spectrum (eigenvalues of ρ_A are smooth across the threshold) but in the OPERATOR CONTENT (eigenvectors) of H_E:
- **Entanglement spectrum is smooth:** tail weight changes 8% from nA=5→6, but 1-R² jumps 212× (q=2 n=14). Confirmed for q=3,5 too.
- **BW residual dominated by max-range operators:** at nA=6, 96% of residual in operators spanning ≥4 sites.
- **Adding longer-range Potts operators barely helps:** free-fit with all Potts ops (21 params) gives 1-R²=0.163, vs BW 0.171. Only 5% improvement.
- **The threshold is non-Potts operator content.** Non-Potts fraction: 6.7e-4 (nA=5) → 0.163 (nA=6) — a 244× jump.
- **BW envelope is optimal for Potts subspace.** Free NN coefficients only 2× better than BW envelope at small nA.
- **Implication:** BW accuracy is limited by operator algebra completeness, not entanglement structure or operator range. To improve beyond BW, need mixed (clock×shift) operators.

**H_E operator compactness — BW breakdown is fundamental (Sprint 097).** Full operator basis decomposition (Pauli for q=2, clock-shift for q=3):
- **Pre-threshold: H_E is maximally compact.** Exactly 2nA-1 BW operators capture >99.9%. Participation ratio 4-11. Nothing to add.
- **At threshold: compactness transition.** Non-BW content erupts (0.07%→18% for q=2, 0.08%→22% for q=3). PR jumps to 12 (total), 357 (non-BW). Need 1052 operators for 99% R² at q=2 nA=6.
- **Non-BW content is DIFFUSE.** 743 operators for 90% of non-BW weight. No compact "BW + corrections" ansatz exists.
- **XZ (mixed) operators dominate:** 44% of non-BW at q=2 threshold. High-body (4-6), long-range (3-5) operators.
- **Universal across q=2,3.** Non-Potts fraction larger for q=3 (22.4%) than q=2 (16.3%), consistent with exp(1.6q) growth.
- **BW is the optimal compact H_E approximation.** The BW story is now complete: it works perfectly until operator algebra mismatch becomes O(1), then fails fundamentally.

**Rényi entropy decomposition of walking (Sprints 085-086).** Size-pair extraction c_α = 6·ΔS/((1+1/α)·ln(N₂/N₁)) on periodic BC. Optimal α shifts gradually with q:

| q | pair | c₁/Rec | c₂/Rec | c₃/Rec | c_∞/Rec | best α |
|---|------|--------|--------|--------|---------|--------|
| 2 | (10,14) | 0.996 | 1.055 | 1.114 | 1.101 | 0.5 (1.000) |
| 3 | (8,10) | 0.994 | 1.081 | 1.131 | 1.080 | 0.5 (1.001) |
| 5 | (6,8) | 1.004 | 1.114 | 1.127 | 1.040 | 1 (1.004) |
| 6 | (6,8) | 0.997 | 1.102 | 1.098 | 1.006 | 1 (0.997) |
| 7 | (6,7) | 0.817 | 0.904 | 0.889 | 0.809 | 2 (0.904) |
| 8 | (6,7) | 0.799 | 0.877 | 0.852 | 0.772 | 2 (0.877) |

**~~α=3 as optimal probe for broken walking~~ RETRACTED (Sprint 086).** Sprint 085's single-size extraction c_α = 12·S/((1+1/α)·ln(N/π)) is contaminated by non-universal constant c'_α which varies with α. Size-pair extraction (cancels c'_α) shows α=1-2 is best for all q. The "α=3 recovery" was c'₃ coincidentally compensating the walking deviation.

**No Rényi index recovers Re(c) for broken walking (q≥7).** Best size-pair c_α/Re(c) ≈ 0.90 at q=7, 0.88 at q=8. Walking breakdown is a genuine multi-eigenvalue phenomenon.

**DMRG confirms walking breakdown in Rényi c₁ (Sprint 086b).** q=7 size-pair c₁/Re(c) drops from 1.05 (n=6,8) to 0.83 (n=10,12). q=5 drops only 5% over same range.

**Open-BC CC profile fit unreliable for Rényi extraction** (R² < 0.9 at n≤12, Sprint 086a). Size-pair extraction from midchain entropy is more robust.

**c_∞ does NOT recover Re(c) at large q** — disproves the hypothesis that min-entropy (probing only λ_max) is insensitive to walking breakdown.

**POTENTIALLY NOVEL (Sprints 084-086):** First entanglement spectrum decomposition across walking boundary. Entropy concentration mechanism. First systematic Rényi c_α(q,α) mapping via size pairs. Demonstration that optimal Rényi index shifts 0.5→1→2 across walking boundary.

**✅ CONFIRMED NOVEL: Fidelity susceptibility α(q) across walking boundary (Sprints 102-103).** χ_F_max ~ N^α. Complete α(q) curve:

| q | α_meas | α_∞ (log-corrected) | #sizes | regime |
|---|--------|--------------------:|--------|--------|
| 2 | 0.980 | — | 5 | real CFT (exact ν=1.0) |
| 3 | 1.43 | 1.22 | 7 | real CFT (exact ν=5/6), strong log correction |
| 4 | 1.77 | ~2.0* | 8 | BKT (extremely slow convergence) |
| 5 | 2.09 | 2.08 | 6 | walking (ZERO correction) |
| 6 | 2.37 | — | 2 | broken walking |
| 7 | 2.65 | — | 3 | broken walking |

*q=4 α_∞=2.0 from exact ν=2/3, but inaccessible at finite size due to BKT log corrections.

**α(q) = 0.315·q + 0.469 for q≥5 ONLY** (walking regime, zero FSS corrections). For q≤3, finite-size α is inflated by power-law 1/N² corrections (Sprint 109, correcting Sprint 108). Power-law fits recover exact ν to sub-percent: q=2 α_∞=1.001 (exact 1.0), q=3 α_∞=1.405 (exact 1.4). q=4 BKT: neither power-law nor log model converges — would need n>100. Walking pairwise α at q=5 converges upward: (6,7)→2.076, (7,8)→2.082, (8,9)→2.090, (9,10)→2.099.

**Gap ratio (multiplet/spectral) decreases with q** (Sprint 108c): 6.2 (q=3) → 5.2 (q=4) → 4.6 (q=5). The χ_F-dominant singlet approaches the symmetry-forbidden spectral gap as q increases.

**✅ CONFIRMED NOVEL: χ_F spectral mechanism (Sprints 106-107).** χ_F is 100% dominated by the (q-1)-fold degenerate S_q multiplet for q≥3. The spectral gap (level 1) has ZERO matrix element with H_field — symmetry-forbidden. α decomposes exactly:

α = β_me + 2z_m - 1, where gap_m ~ N^{-z_m} and |⟨mult|H_field|0⟩|² ~ N^{β_me}

| q | sizes | z_m | β_me | α(decomp) | α(known) |
|---|-------|-----|------|-----------|----------|
| 2 | 5 (n=6-14) | 0.989 | 0.066 | 1.044 | 1.099 |
| 3 | 3 (n=6-10) | 1.031 | 0.381 | 1.442 | 1.414 |
| 5 | 5 (n=6-10) | 1.155 | 0.776 | 2.085 | 2.044 |
| 7 | 3 (n=6-8) | 1.312 | 1.025 | 2.649 | 2.674 |

Linear fits: z_m(q) = 0.065q + 0.845, β_me(q) = 0.188q − 0.238. Reconstructed α(q) = 0.318q + 0.452 matches direct fit exactly. Cross-validated: spectral/finite-diff ratio confirms 98% captured for q=5 (single multiplet) vs 87% for q=2 (multi-state). Energy and entanglement multiplet gaps share S_q symmetry origin but scale independently (z_energy ≈ 1.0-1.3 vs z_ent ≈ 0.2).

**✅ CONFIRMED NOVEL: Casimir energy obeys complex CFT Re(c) (Sprints 083, 098).** E₀/N = ε_∞ - πv·Re(c)/(6N²) + O(1/N⁴). Hardened with GPU-extended sizes (Sprint 098):

| q | sizes | c_Cas/Rec(pw_last) | c_eff/Rec | c_Cas/Rec(extrap∞) |
|---|-------|-------------------|-----------|-------------------|
| 2 | 6-14 (5pts) | 0.986 | 1.000 | 0.984 |
| 3 | 4-12 (5pts) | 0.988 | 1.116 | 0.983 |
| 5 | 4-10 (4pts) | 1.010 | 1.012 | 1.008 |
| 7 | 4-8 (5pts) | **1.000** | 0.784 | 1.008 |
| 8 | 4-7 (4pts) | 0.992 | 0.739 | 1.008 |

**Pairwise-last c/Rec: mean=0.995, std=0.009. Casimir is 16× more consistent with Re(c) than entropy (std=0.145).**

1/N⁴ corrections: |d|/q ≈ 0.044 (q-independent, universal lattice correction). 2-param fits biased 2-4%, but pairwise at largest sizes and extrapolated ∞ values all within 1-2% of Re(c).

**Scope:** Holds at all accessible sizes (n≤14 for q=2, n≤8 for q=7). Walking breakdown is exclusively an entropy phenomenon — energy sees Re(c) everywhere entropy fails.

**Energy-entropy hierarchy is NOT universal (Sprint 104).** Tested in J1-J2 spin-1/2 chain at three regimes:

| Model | c | c_eff(18,20) | c_Cas(18,20) | Winner |
|-------|---|-------------|-------------|--------|
| XX (Δ=0) | 1.000 | 0.994 | 1.006 | tied |
| Heisenberg (Δ=1) | 1.000 | 0.999 | 1.014 | entropy 21× |
| J1-J2 BKT (J2=0.2412) | 1.000 | 0.996 | 1.000 | Casimir 14× |

Hierarchy direction depends on correction type: (1) marginal operator logs pollute Casimir → entropy wins (Heisenberg), (2) BKT log³ corrections pollute entropy → Casimir wins. **Max J1-J2 entropy deviation = 0.6% vs Potts walking = 26% — the Potts effect is 41× larger.** The O(1) entropy deviation with accurate Casimir is unique to the walking mechanism (entanglement spectrum reorganization).

**Im(c) oscillations undetectable at accessible sizes (Sprints 099-100).** Complex CFT predicts oscillatory corrections ~cos(ω·ln(N))/N⁴ with period ~e^300 for q=5. Literature confirms: no paper has detected Im(c) from Casimir energy of a real Hamiltonian. Only non-Hermitian deformation (Tang et al. 2024) accesses Im(c) directly.

**Open-BC Casimir extraction FAILS for q≥5 (Sprint 100).** DMRG gives open-BC energies. The open-BC formula E₀ = ε_∞·N + e_s - πvc/(24N) has O(1/N) boundary corrections that scale dramatically with q: c/Re(c) = 0.97 (q=2, N≤24), 0.64 (q=5, N≤12), 0.36 (q=7, N≤8). Periodic BC is 10-100× more accurate at same N (1/N² vs 1/N corrections). DMRG cannot extend the Casimir-Re(c) result beyond exact diag sizes.

**Open BC correlators give WRONG exponents (Sprint 082a).** Raw power-law fit on DMRG open chain inflates η by ~5×. NEVER use open-BC raw power law for x extraction — use periodic chain with chord distance.

**S_q Potts g_c = 1/q exactly (Sprint 078a, self-duality).** Kramers-Wannier duality gives g_c = 1/q. Verified numerically: all crossings approach 1/q from above with ~1/n² corrections. q=2 deviation 0.05% at (10,12), q=5 deviation 0.4% at (6,8).

**S_q Potts critical points (Sprints 076-078):**

| q | g_c (S_q) = 1/q | g_c (hybrid) | g_c ratio | S_q 1st degen | Hybrid 1st degen |
|---|-----------------|-------------|-----------|---------------|-----------------|
| 2 | 1/2 = 0.500 | 0.250 | 2.00 | identical (Ising) | identical |
| 3 | 1/3 = 0.333 | 0.333 | 1.00 | identical | identical |
| 5 | 1/5 = 0.200 | 0.441 | 2.21 | 4-fold (S₅) | 2-fold (Z₅) |
| 7 | 1/7 = 0.143 | 0.535 | 3.75 | 6-fold (S₇) | 2-fold (Z₇) |
| 10 | 1/10 = 0.100 | 0.684 | 6.84 | (q-1)-fold | 2-fold |

**q=2: S_q ≠ hybrid (Sprint 078a).** S_q field = Σ X^k = X (one operator). Hybrid field = X+X† = 2X. Factor-of-2 field difference gives S_q g_c = 0.500, hybrid g_c = 0.250. S_q = hybrid ONLY at q=3.

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

## Symmetry-Resolved Entanglement Entropy (Sprint 101)

**S_q Potts at g_c=1/q, periodic BC.** SREE decomposes ρ_A into Z_q charge sectors via P_α = (1/q) Σ_k ω^{-αk} G_A^k.

**S_number/S_total is q-independent at fixed geometry (Sprint 101c).** At n=4 nA=2: ratio=0.953 for ALL q=2-7. At n=6 nA=3: 0.908. At n=8 nA=4: 0.875. Decreasing with n (→0 thermodynamically since S_n=O(1), S_t=O(log n)). Implies S_number ∝ c(q) at fixed size.

**Charge-0 always enriched: p(0)*q >> 1.** Increases monotonically with q: 1.55 (q=2) to 5.15 (q=10) at n=6. p(0) → 1/2 as q→∞ (not 1/q). Convergence to equipartition requires n >> 100.

**Not a walking discriminator** at accessible sizes — all q-dependence is smooth and monotonic. No feature at q=4 or q=5.

## Fidelity Susceptibility (Sprint 102)

**χ_F(g) = (2/N)(1 - |⟨ψ(g)|ψ(g+δg)⟩|)/δg².** Measures ground state sensitivity to coupling change. At continuous transition: χ_F_max ~ N^{2/ν-1}. At first-order: χ_F_max ~ N².

**Validates real CFT exponents for q=2,3.** q=2: ν=1.009 (exact 1.0, <1% error). q=3: ν=0.841 (exact 0.833, <1% error).

**Walking gives anomalous first-order-like scaling (q=5).** Scaling exponent α=2.09 (from n=6,8 only), giving ν_eff=0.648 instead of gap-derived 0.83. α≈2 matches first-order prediction for 1D. This is 22% discrepancy, far beyond the 1% FSS error at q=2,3.

| q | ν_exact | α_expected | α_measured | ν_fidelity | deviation |
|---|---------|------------|------------|------------|-----------|
| 2 | 1.000 | 1.00 | 0.98 | 1.009 | -2% |
| 3 | 0.833 | 1.40 | 1.38 | 0.841 | -1% |
| 4 | 0.667 | 2.00 | 1.69 | 0.743 | -16% |
| 5 | 0.83* | 1.41 | **2.09** | 0.648 | **+48%** |

*gap-derived ν

**χ_F grows exponentially with q at fixed n.** At n=6: 0.73 (q=2), 4.09 (q=3), 13.98 (q=4), 36.29 (q=5), 166.88 (q=7). Fit: ~exp(1.06q), growth rate 2.88× per unit q.

**No finite-size peak shift for q≥5.** Peak sits exactly at g_c=1/q. For q=2,3 the peak approaches g_c from below as ~1/N. Zero shift consistent with self-duality fixing g_c exactly + sharp transition.

**CAVEAT:** q=5 result based on only 2 sizes (n=6,8). Need more sizes to confirm. q=4 deviation may reflect logarithmic corrections at the marginal case.
