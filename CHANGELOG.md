# Changelog — Quantum Explorer

## QPU Budget
- Used: 20s of 600s (Sprint 025)
- Remaining: 580s

## Compressed Sprint History (Sprints 001-024)

- **001** — Bell states & CHSH. S=2.834 matches Tsirelson bound. GHZ scaling timed out (partial_trace limit).
- **002** — CHSH landscape (19.9% of angles violate), noise kills at p≈9.5%, cluster states gain entropy under qubit loss.
- **003** — Cluster state qubit loss: position-dependent, two-qubit loss unlocks 3 bits of hidden entanglement.
- **004** — MI and I3 distinguish all archetypes. GHZ has I3=+1 (redundant), Cluster has I3=-1 (irreducible).
- **005** — Discord: GHZ=0 (classical pairwise), W=0.28 (most quantum), Cluster=0. Three archetypes confirmed.
- **006** — Concurrence, negativity, monogamy. Cluster has 7x GHZ entanglement at certain bipartitions.
- **007** — 2D cluster: topological entanglement, position-independent qubit loss, 3x baseline entropy.
- **008** — Phi/IIT: 2D cluster 100% Phi retention under loss. Phi measures structure not quantumness.
- **009** — Structured noise fingerprints: each state has unique 6D signature. Cluster weakest to dephasing.
- **010** — GME witnesses, body-order hierarchy, 4-measurement classifier (182x compression vs tomography).
- **011** — Entanglement dynamics: CZ gates can destroy MI. 2D cluster annihilates ALL pairwise MI.
- **012** — Scrambling & OTOCs: GHZ has OTOC=1 despite MI=15. Entanglement ≠ scrambling.
- **013** — Hayden-Preskill: scrambling democratizes recovery. Page curve visible. Light cone in 1D.
- **014** — QEC = static scrambling. [[5,1,3]] Page curve matches Hayden-Preskill. Code distance = body-order.
- **015** — Code families: distance ≠ topology. Democratic/Selective/Hierarchical code archetypes.
- **016** — Structured noise breaks topology-blindness. "Keep and correct" vs "filter and project" strategies.
- **017** — Threshold as phase transition. Holevo is the order parameter. Concatenation sharpens transition.
- **018** — Syndrome information: lossy compression. Fault tolerance is emergent collective phenomenon.
- **019** — Channel capacity: depolarizing threshold p≈0.20. Codes far from capacity. Phase damping non-monotonic.
- **020** — Toric code: zero I3, local indistinguishability, bimodal Page curve. Fourth archetype (Topological).
- **021** — Combined T1+T2: basis isotropy is everything. [[5,1,3]] asymmetry 0.005 vs 3-qubit 0.680.
- **022** — Noise-adapted codes: isotropy always wins for quantum computation at small scale.
- **023** — Concatenated bias-tailoring: figure of merit determines winner. Min-basis Holevo → [[5,1,3]] 81%.
- **024** — Surface code at d=3: never wins. Boundary tax breaks isotropy. Advantage is architectural, not informational.

- **025** — Real hardware QPU test (20s ibm_kingston). [[5,1,3]] isotropy confirmed: asymmetry 0.040 vs 3-qubit 0.254.
- **026** — Active syndrome extraction: non-FT always hurts at all error rates.
- **027** — Flag-FT syndrome: solves the wrong problem, still never beats passive.
- **028** — Repeated syndrome: 48 2Q gate overhead kills. QEC active correction arc closed.
- **029** — TFIM phase transition: Scale-Free archetype discovered. MI-CV as order parameter.
- **030** — XXZ archetype loop. MI-CV classifies transition type (crossing/step/dome).
- **031** — Entanglement spectrum: three levels of description (scalar/correlation/spectral).
- **032** — Bisognano-Wichmann: symmetry controls BW accuracy. Fourth level (Hamiltonian).
- **033** — Potts BW: local dimension matters more than group size.
- **034** — BW operator algebra: H/G-inv ratio perfectly predicts BW locality.

Full details for all compressed sprints are in sprints/sprint_NNN.md.

---

## Detailed Sprint Log (Recent)

- **035-043** — BW scaling, MI-CV FSS/classification, Potts phase transitions. Many results invalidated by Sprint 049 (wrong g_c).
- **044** — g_c scaling law g_c∝(q-3)^{-0.85}. Blind prediction q=7 accurate to 1.6%. (INVALIDATED Sprint 049)
- **045** — ν(q=5)≈2 from MI-CV data collapse. Gell-Mann MI WRONG for d≥5. (ALL ν VALUES WRONG, Sprint 053)
- **046** — ν(q=4)≥2.2, BKT-like. Gell-Mann misleading at d=4. (WRONG, Sprint 053)
- **047** — ν(q=7)≈0.5. q=4 confirmed unique BKT peak. (WRONG, Sprint 053)
- **048** — Dead-pair bias invalidates raw MI-CV at d≥10. χ convergence ≠ MI convergence.

- **049** — q=3 Potts g_c WRONG (was at 1.0, true ≈ 0.33). 10 sprints of Potts results invalidated. Entropy FSS gives c not ν.
- **050** — True Potts critical points via self-duality (q=2: g_c=0.25, q=3: g_c=1/3 exact). MI-CV crossings vindicated at true g_c.
- **051** — Energy gap method: g_c INCREASES with q (0.25→0.33→0.39→0.44→0.52). Reverses 10 sprints of wrong results.
- **052** — g_c(q) ≈ (1/5)√(q-1) + 1/20. q=10 measured (g_c≈0.684). POTENTIALLY NOVEL.
- **053** — ν(q) at true critical points: ν(q=3,4,5) ≈ 0.82-0.86, nearly constant. All old ν values WRONG.

Full details for all compressed sprints are in sprints/sprint_NNN.md.

- **054-055** — c(q) from entropy scaling. c grows with q: 0.5→0.8→1.0→>1. Entropy profile best method. iDMRG fails at criticality.
- **056** — c(q) ≈ 0.40·ln(q-1) + 0.55. Analytic continuation and quadratic interpolation RULED OUT for q>4.
- **057** — CFT operator content: Z_q spin field degeneracy, harmonic ratios sub-quadratic, free boson is q→∞ limit. POTENTIALLY NOVEL.
- **058** — x₁ peaks at q=3. c/x₁ = 2q for q=2,3 only. Δ₁·N sign flip at q≥5 (anomalous FSS).
- **059** — Conformal towers confirm genuine CFT for ALL q=2-10. Descendant degeneracy=4 exact for q≥3.
- **060** — OPE C_{σ*,σ,ε} peaks at q=3 (~0.54), decreases ~q^{-0.8}. Z_q charge detection fails for conjugate pairs.

Full details for all compressed sprints are in sprints/sprint_NNN.md.

- **061** — Large-q limit: harmonic ratios → k² as O(q^{-3}). x₁ decreases to 0.038 at q=30. R_ε grows linearly with q. POTENTIALLY NOVEL.
- **062** — c·x₁ ≈ 0.112 for q≥3 (near 1/9). Single boson RULED OUT (c>1). Twist field = σ for q=2 exactly.
- **063** — c·x₁ NOT exactly 1/9 (q=3 exact: 8/75). Potts-specific, not Z_q-universal. Clock c > Potts c at same q.
- **064** — BW locality degrades monotonically with q (~15% drop per unit q). G-inv fraction = exactly 1/q. Peak shifts to ordered phase for q≥4.

- **065** — Hybrid ≠ clock universality: 3 probes at q=5 (c/x₁, ν, gap crossing). Clock BKT, hybrid power-law. NEW universality class. POTENTIALLY NOVEL.
- **066** — No first-order signal at q=10 (n≤6). Δ·N decelerates, gap minimum size-independent. Continuous confirmed.
- **067** — Hybrid has ONE transition, no floating phase. Clock has wide floating phase (Δg=0.62). Fourth hybrid-clock difference.

Full details for all compressed sprints are in sprints/sprint_NNN.md.

### Sprint 068 — 2D Hybrid Model: First Study, No First-Order Signal at q=5
**Status:** Complete (3 experiments).

**First study of Potts-clock hybrid in 2D.** Built Hamiltonian on periodic L×L square lattices. Gap×L crossing method validated on q=2 (2D Ising, L=2,3,4).

**2D critical points measured.**
| q | g_c(2D) | 2D/1D ratio | Sizes |
|---|---------|-------------|-------|
| 2 | 0.771 | 3.08 | L=2,3,4 |
| 3 | 1.267 | 3.80 | L=2,3 |
| 5 | 1.588 | 3.60 | L=2,3 |

**q=2 continuous transition confirmed (068b).** L=3 and L=4 gap×L agree to 4 decimal places (2.3632 vs 2.3631). L=2 completely out of scaling regime (gap×L=4.19).

**q=5: no first-order signal detected (068c).** dE₀/dg smooth across g_c (no kink). Z₅ conjugate pair degeneracy exact in 2D. But only L=2,3 available — **inconclusive** for definitive determination. Need larger lattices.

**Gap×L crossing ratio is circular with 2 sizes.** gap(L₂)/gap(L₁) = L₁/L₂ at the crossing BY DEFINITION. Not an independent z=1 test. Three sizes needed (only available for q=2).

**POTENTIALLY NOVEL:** First measurements of 2D g_c for the hybrid model. If continuous behavior at q=5 is confirmed by larger lattices (QMC), the hybrid would be unique: continuous in 2D where both standard Potts (first-order) and clock (BKT) show different behavior.

**Surprises:**
- L=2 is out of scaling regime for q=2 2D Ising (gap×L off by 77%)
- gap₂/gap₁ = 1.000 for ALL q,L,g in 2D — Z_q degeneracy is exact
- 2D/1D ratio is non-monotonic: 3.08, 3.80, 3.60 for q=2,3,5
- Gap×L at g_c varies widely: 2.36 (q=2), 6.10 (q=3), 3.50 (q=5)
- L=3 at q=5 (dim=2M) runs at 25s/pt — practical 2D exact diag

[Full report: sprints/sprint_068.md]

### Sprint 069 — 2D Entanglement Entropy: Area Law and Ordered-Phase Dominance
**Status:** Complete (3 experiments).

**Ordered phase entropy = ln(q) EXACTLY.** For all q=2,3,5 and all L=2,3,4 tested. The Z_q-symmetric ground state is a GHZ-like superposition with S = ln(q) independent of system size and bipartition geometry.

**Critical entropy follows area law.** S/boundary at w=1 strip, L=3: q=2 (0.065), q=3 (0.039), q=5 (0.053). Non-monotonic in q. L=2 out of scaling for all q.

**Entropy peak is NOT at g_c.** At accessible sizes (L≤4), the GHZ entropy ln(q) dominates over the critical area-law entropy. For q=2 L=4: peak at g=0.58 vs g_c=0.771. For q=5 L=3: peak at g=0.68 vs g_c=1.588. Crossover at L ~ ln(q)/(2α) ≈ 10-16.

**Sharp entropy drop at g ≈ 0.6·g_c** for q=5 L=3. S drops from 1.62 to 0.33 between g=0.68 and g=1.57.

**Surprises:**
- Entropy peak in ordered phase, NOT at criticality — qualitatively different from 1D
- S_ordered = ln(q) exactly, independent of L and bipartition
- Area-law coefficient non-monotonic in q (q=3 lowest at 0.039)
- Sharp entropy crossover well below g_c

[Full report: sprints/sprint_069.md]

### Sprint 070 — 2D Energy Derivatives & Fidelity: q=2 Continuous, q=5 Inconclusive
**Status:** Complete (3 experiments).

**q=2 2D Ising continuous transition confirmed (L=3→4).** d²E₀/dg² peak scales as L^0.16 (consistent with α=0). Fidelity susceptibility χ_F/N scales as L^0.94 (consistent with ν=1). dE₀/dg smooth — no latent heat.

**q=5 2D transition nature INCONCLUSIVE.** Only L=2,3 accessible; L=2 pathologically out of scaling (all quantities 10-50x too small). dE/dg smooth at L=3 (no latent heat signal). F_min=0.985 (lower than q=2,3 but not diagnostic). d²E peak smaller at q=5 than q=2 — less singular, not more.

**Eigenstate-sum χ_F FAILS for dim > 5000.** Truncation to k=10 states captures negligible fraction. Ground state overlap method works correctly for any dim.

**L=2 out of scaling for ALL diagnostics** — not just gap. L=2→3 apparent exponents (5-7) are artifacts.

**Fundamental bottleneck:** q=5 L=4 dim ≈ 10^11 — infeasible with exact diag. Cylinder DMRG or QMC needed.

**Surprises:**
- d²E/dg² and χ_F peaks are in ordered phase (far below g_c) for all q at small L
- q=5 d²E peak (1.58) SMALLER than q=2 (3.00) at L=3 — opposite to first-order expectation
- Overlap-derived χ_F matches eigenstate sum only at full diag (dim < 5000)
- dE/dg converges smoothly for all q — no hint of latent heat even at q=5

[Full report: sprints/sprint_070.md]

### Sprint 071 — 2D Cylinder Geometry: Ly=2 Ladder for q=2,5
**Status:** Complete (3 experiments).

**Ly=2 cylinder (ladder) extends 2D study beyond torus bottleneck.** Open x-direction, periodic y-direction. Coordination z=3. Enables exact diag gap×Lx crossings for moderate sizes.

**g_c(cylinder) measured for q=2 and q=5.**
| q | g_c(1D) | g_c(cyl) | g_c(2D) | cyl/1D |
|---|---------|----------|---------|--------|
| 2 | 0.250 | 0.451 | 0.771 | 1.80 |
| 5 | 0.441 | 0.714 | 1.588 | 1.62 |

**q=2 crossings converge (3 pairs): 0.447→0.451→0.454.** q=5 has single pair (3,4)=0.714.

**No first-order signal for q=5 on cylinder (071c).** Order parameter ⟨δ(s_i,s_j)⟩ smooth across g_c. Max slope steeper for q=5 (1.88) vs q=2 (1.21) — more critical, not less. Extends Sprint 066 finding from 1D to ladder geometry.

**DMRG impractical for q≥5 cylinder (071b).** d=5 per site → chi=20 has massive truncation errors. Multi-hour runtimes even at chi=30. Exact diag is the viable approach for Lx≤4-5.

**DMRG on q=2 cylinder works (071a).** g_c≈0.45 from entropy peak (drifts with Lx). c_eff=0.19 truncation-limited.

**Surprises:**
- DMRG entropy peak drifts strongly on cylinders — unreliable g_c proxy
- q=5 DMRG completely impractical at chi=20-30
- Cyl/1D ratio appeared non-monotonic with 2 points (1.80 vs 1.62) — Sprint 072 shows monotonically decreasing (1.80→1.69→1.62)
- q=5 order parameter slope STEEPER than q=2 — opposite to first-order expectation

[Full report: sprints/sprint_071.md]

### Sprint 072 — q=3 Cylinder & Ly=3 Extension
**Status:** Complete (3 experiments).

**q=3 fills the Ly=2 cylinder gap (072a).** g_c(q=3, Ly=2 cyl) = 0.565 from 3 crossing pairs (Lx=3-6, dim up to 531k): 0.555→0.567→0.573. Monotonic convergence from below. Cyl/1D ratio = 1.69.

**Cyl/1D ratio monotonically decreasing (072b).** 1.80 (q=2) → 1.69 (q=3) → 1.62 (q=5). NOT non-monotonic as previously noted for 2D/1D. Higher q needs proportionally less coupling boost on cylinder.

**q=3 order parameter smooth (072b).** Max slope 1.99 (Lx=5,6) — steepest of the three q values tested. No first-order signal.

**Ly=3 cylinder for q=2 (072c).** g_c = 0.655 from 4 crossing pairs (Lx=3-7, dim up to 2M): 0.645→0.655→0.659→0.662. 63.7% of the way from Ly=2 (0.451) to 2D (0.771). Gap×Lx at crossing converges toward ~2.4.

**Surprises:**
- Cyl/1D ratio is monotonically decreasing, not non-monotonic (corrects Sprint 071 note)
- q=3 has steepest order parameter slope (1.99 > q=5's 1.88 > q=2's 1.21)
- 2D/cyl ratio peaks at q=3 (2.24) — the 2D jump from cylinder is largest for q=3
- Ly=3 crossings converge tighter than Ly=2 (spread 0.017 in 4 pairs)

[Full report: sprints/sprint_072.md]

### Sprint 073 — Ly Convergence: q=3 Ly=3, q=2 Ly=4, Cross-q Analysis
**Status:** Complete (3 experiments).

**q=3 Ly=3 cylinder (073a).** g_c = 0.797 from single crossing pair (Lx=3,4). Only 49.7% of the way to 2D (vs q=2's 77.7% at Ly=3). Lx=5 (dim=14.3M) at 251s/pt — infeasible for full scan.

**q=2 Ly=4 cylinder (073b).** g_c = 0.688 from 2 crossing pairs: (3,4)=0.683, (4,5)=0.693. 84.0% to 2D. Convergence decelerating: Ly=3→4 jump (+0.033) is 6x smaller than Ly=2→3 (+0.204). Exponential fit predicts Ly=5: g_c=0.745.

**Cross-q convergence analysis (073c).** q=3 converges to 2D **1.9x SLOWER** than q=2. Exponential decay lengths: B(q=2)=1.60, B(q=3)=3.03. q=3 has both larger 2D/1D ratio (3.80 vs 3.08) AND slower exponential approach. Ly convergence is NOT universal across q — different models require different cylinder widths to approximate 2D behavior.

**Surprises:**
- q=3 at Ly=3 only halfway to 2D — much slower than q=2
- q=2 Ly=3→4 jump surprisingly small (+0.033) — near saturation
- q=2 convergence ratio non-monotonic (0.61, 0.36, 0.72) — not pure exponential
- 2D/1D g_c ratio correlates with convergence slowness: bigger gap = slower approach

[Full report: sprints/sprint_073.md]

### Sprint 074 — Cylinder Limits: q=3 Ly=4 Infeasible, q=5 Ly=3 g_c=0.974, Convergence Saturates
**Status:** Complete (3 experiments).

**q=3 Ly=4 cylinder INFEASIBLE (074a).** Lx=4 (dim=43M) takes 838s per g-point — no scan possible. Only Lx=3 available (no crossing). Exponential fit prediction g_c ≈ 0.917 cannot be tested.

**q=5 Ly=3 cylinder (074b).** g_c = 0.974 from (Lx=2, Lx=3) crossing. 46.5% to 2D. Lx=4 (dim=244M) infeasible. Caveat: Lx=2 known to be out of scaling regime.

**Convergence saturates for q≥3 (074c).** At Ly=3: q=2 at 77.7%, q=3 at 49.7%, q=5 at 46.5%. The big slowdown is q=2→q=3 (B: 1.19→2.49). q=3→q=5 is marginal (B: 2.49→2.83). q=3 Ly increments exactly constant (0.232) — possibly linear, not exponential.

**Surprises:**
- q=3 Ly=4 hits exact diag wall at dim=43M (838s/pt)
- q=3 and q=5 converge at nearly identical rates (~47-50% at Ly=3)
- Convergence slowdown saturates — q≥3 may all approach 2D at similar Ly
- q=3 increments exactly 0.232 for both Ly=1→2 and Ly=2→3 (linear!)
- Cyl/1D ratio decreasing at both Ly=2 AND Ly=3 — universal trend

[Full report: sprints/sprint_074.md]

### Sprint 075 — q=2 Ly=5 Cylinder & Cylinder Entanglement Entropy
**Status:** Complete (3 experiments).

**q=2 Ly=5 cylinder (075a).** g_c = 0.701 from (Lx=3, Lx=4) crossing. 86.6% to 2D. Ly=4→5 increment only +0.013 (Ly=3→4 was +0.033). Lx=5 (dim=33M) takes 541s/pt — infeasible. Exponential prediction (0.745) OVERSHOOTS by 5.9%.

**Cylinder entanglement entropy (075b).** First measurement of S on cylinder geometries at criticality. c_eff from open-BC fit: Ly=1 (0.61), Ly=2 (0.75), Ly=3 (0.82). Growth is FSS overshoot — 2D Ising also has c=0.5. S per bond cut DECREASING: 0.285→0.157→0.107→0.089. Area-law behavior emerging.

**Power-law convergence, NOT exponential (075c).** g_c(Ly) = 0.771 - 1.285/Ly^2.03 fits better (RMS=0.016) than exponential (RMS=0.025). Ly=1→3 increments are nearly linear (+0.201, +0.204), abrupt deceleration at Ly≥4. This supersedes the exponential model from Sprint 073.

**Surprises:**
- Exponential prediction overshoots by 5.9% — convergence slower than expected
- Ly=1→2 and Ly=2→3 increments nearly identical (+0.201, +0.204) — linear regime
- Abrupt deceleration at Ly=3→4 (6x smaller increment)
- c_eff at Ly=1 already 22% above exact c=0.5 — significant small-Lx overshoot
- S per bond cut monotonically decreasing — approach to area law

[Full report: sprints/sprint_075.md]

### Sprint 076 — S_q Potts vs Hybrid: Direct Universality Class Comparison at q=5
**Status:** Complete (3 experiments).

**First head-to-head comparison of S_q Potts vs Potts-clock hybrid at q=5.** Built standard S_q Potts Hamiltonian (field = Σ_{k=1}^{q-1} X^k, full S_q symmetry) and compared to hybrid (field = X+X†, Z_q symmetry).

**S_q Potts g_c(q=5) = 0.200 (076a).** From gap×N crossings at n=4,6,8. Verified S_q ≡ hybrid at q=3 (10-digit match). S_q field strength (q-1)=4 per unit g, hybrid 2cos(2π/5)≈0.618. Ratio 6.47x → g_c ratio 2.19x.

**Qualitatively different conformal towers (076b).** S_q has 4-fold degenerate first excited state (S₅ symmetry). Hybrid has 2+2 split (Z₅ conjugate pairs, x₃=2.41). x₃ is 141% different between models. Both show clean CFT convergence at n≤8.

**S_q Potts FSS corrections 10x larger than hybrid (076c).** Δ·N drift n=4→8: S_q 4.8% vs hybrid 0.46%. Power-law fits better than exponential for both (25x for S_q, 50x for hybrid). q=7 S_q gap collapses 69% from n=4→6 (possible first-order signal at larger q).

**Surprises:**
- S_q Potts at q=5 looks like genuine CFT at n≤8 — no first-order signature visible
- Degeneracy structure is sharpest discriminator (exact integers, no fitting needed)
- Δ₂/Δ₁ = 1.0000 exactly for BOTH models — conjugate pair degeneracy universal
- q=3 verification perfect: confirms code correctness

[Full report: sprints/sprint_076.md]
