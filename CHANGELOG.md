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

### Sprint 056 — c(q) Formula: Logarithmic Growth Confirmed, Analytic Continuation Ruled Out
**Status:** Complete (4 experiments).

**Analytic continuation of Potts CFT RULED OUT for q>4.** Both Coulomb gas and minimal model continuations predict c DECREASING below 1 for q>4. Measured c(7)≈1.46 raw and c(10)≈1.60 raw are far above 1.0. The 1D quantum Potts at q>4 is NOT described by the standard Potts CFT.

**Quadratic interpolation also RULED OUT.** Fits exact q=2,3,4 and accidentally gives c(5)=1.10, but predicts c(7)=1.00 and c(10)=0.10 — contradicted by measurements.

**c(q) grows monotonically, approximately as ln(q).** Best fit: c ≈ 0.40·ln(q-1) + 0.55 (RMS=0.061). No peak, no saturation.

**New c measurements:** c(q=7) profile 1.462 at n=8 (DMRG) and 1.537 at n=6 (exact diag). c(q=10) profile 1.596 at n=6 (exact diag, 10^6 dim Hilbert space, 42s).

**Surprises:**
- Analytic continuation gives c DECREASING for q>4 — opposite to physical data
- Quadratic through q=2,3,4 accidentally gives c(5)=1.10 — a numerical coincidence trap
- Raw n=6 profile for q=10 (1.596) coincidentally matches pre-correction log prediction (1.59)
- q=4 overshoot 8.7% worse than calibration model — log corrections propagate

[Full report: sprints/sprint_056.md]

### Sprint 057 — Operator Content & Scaling Dimensions: What CFT Describes q>4 Potts?
**Status:** Complete (6 experiments).

**CFT spectrum method validated.** Energy gap ratios on periodic chains give scaling dimension ratios. q=2 (Ising): R_ε/R_σ → 8.0 at n=12, matching x_ε/x_σ = 1/(1/8) = 8 exactly. q=3 (W₃): R → 6.0 for energy, 10 for disorder, matching CFT. Absolute x₁ extraction: 0.3% error for q=2,3.

**Z_q spin field degeneracy pattern DISCOVERED.** Exactly (q-1) spin field primaries organized as ⌊(q-1)/2⌋ conjugate pairs (+ 1 self-conjugate if q even). Confirmed for q=2,3,4,5,7,10. This is the operator content of a previously uncharacterized CFT for q>4.

**Harmonic ratios x(σᵏ)/x(σ) are sub-quadratic.** NOT a free boson (k² prediction fails by 20-139%). Ratios approach k² as q→∞ — the free boson is the large-q limit. No closed-form formula found.

**Energy field R_ε has minimum at q=3 (=6), increases for q>3.** The energy operator is pushed above all spin harmonics.

**POTENTIALLY NOVEL:** Complete operator content of 1D quantum Potts CFT for q>4 appears previously unmeasured. The degeneracy structure, harmonic ratios, and physical mechanism (q-1 primaries → c~ln(q)) are new.

**Surprises:**
- q=4 σ² at R=1.68, not 4 (free boson) or 2 (sin²) — no simple formula
- q=10 σ⁵ is a clean singlet at n=4, confirming even-q self-conjugation
- Energy field ratio R_ε non-monotonic (minimum at q=3)
- Harmonic ratios approach k² from below — free boson is asymptotic limit

[Full report: sprints/sprint_057.md]

### Sprint 058 — x₁(q) Scaling Dimensions: Peak at q=3, Sub-linear c/x₁ Growth
**Status:** Complete (5 experiments).

**New computations:** q=4 at n=10 (dim=10⁶, 15s), q=5 at n=8 (dim=390k, 9s), q=10 at n=6 (dim=10⁶, 19s). Extended Sprint 057 data to larger sizes for precise x₁ extraction.

**c/x₁ = 2q EXACT for q=2,3, but NOT for q≥4.** Model-independent ratio c/x₁ extracted from Casimir energy + gap: 4.01 (q=2), 5.98 (q=3), 8.54 (q=4), 10.84 (q=5), 15.11 (q=7), 16.86 (q=10). Sub-linear growth — q=10 gives 16.9, far below 2q=20. No simple formula fits.

**x₁ peaks at q=3.** x₁(q): 0.125, 0.133, 0.117, 0.101, 0.086, 0.083 for q=2,3,4,5,7,10. Non-monotonic: rises from q=2→3, falls for q>3 toward 0.

**Logarithmic FSS corrections at q=4.** Excess (c/x₁−8)·N² grows as ~N² (not constant), confirming sub-power-law correction. Consistent with marginal Ashkin-Teller point.

**Δ₁·N sign flip at q≥5.** For q≤4, Δ₁·N decreases with N (normal). For q≥5, it INCREASES — anomalous finite-size behavior unique to the novel CFT regime.

**Surprises:**
- c/x₁(q=10) = 16.9, definitively ruling out c/x₁ = 2q for large q
- Δ₁·N sign flip at q=5 — qualitative change in FSS corrections
- q=4 log corrections make 8.0 vs 8.5 distinction impossible at n≤10
- x₁(q=7) ≈ x₁(q=10) ≈ 0.08 — possible saturation at large q

[Full report: sprints/sprint_058.md]

### Sprint 059 — Conformal Tower Analysis: Genuine CFT Confirmed for q>4 Potts
**Status:** Complete (4 experiments).

**Momentum-resolved conformal towers measured for q=2,3,4,5,7,10.** Translation operator T diagonalized within degenerate eigenspaces to extract conformal spin. Method validated on q=2 (Ising) and q=3 (W₃ CFT) with <2.5% descendant gap error.

**Genuine CFT structure confirmed for ALL q.** Four qualitative tests pass universally: (1) descendant degeneracy = 4 (q≥3) or 2 (q=2), exact; (2) descendant momentum = spin ±1; (3) gap converges from below; (4) tower organization Identity→σ pairs→ε→L₋₁σ.

**Descendant gap matches 1/x₁ for q=2,3 (~2% error).** For q≥5, anomalous FSS gives 15-25% gap error at accessible sizes — same sign-flipped corrections as Sprint 058. Qualitative structure perfect despite quantitative gap offset.

**POTENTIALLY NOVEL:** First momentum-resolved conformal tower for q>4 Hermitian Potts chain. Previous work (Lao et al. PRB 2019) found drifting towers assuming first-order; Tang et al. (PRL 2024) needed non-Hermitian model. Our towers converge because transition is genuinely second-order.

**Surprises:**
- Tower structure qualitatively identical for ALL q=2-10 — no transition at q=4
- Descendant degeneracy 4 is exact, not approximate, for all q≥3
- ε field nearly merges with L₋₁σ at large q (R_ε ≈ 9 vs R_desc ≈ 10 for q=10)
- q=10 n=6 (dim=10⁶) tower extraction in 66s — exact diag still viable

[Full report: sprints/sprint_059.md]

### Sprint 060 — OPE Coefficients: C_{sigma*,sigma,epsilon} Peaks at q=3
**Status:** Complete (5 experiments).

**OPE extraction method validated.** Ratio C = |<epsilon|Z_0|sigma>|/|<0|Z_0|sigma>| using clock operator Z (charge +1). Z_q selection rules ensure only charge-0 states (identity, epsilon) couple. q=2 Ising: converges to exact 1/2 with 0.4% error at n=12.

**C_{sigma*,sigma,epsilon}(q) measured for q=2-10.** Peaks at q=3 (~0.54), then DECREASES monotonically: 0.46 (q=4), 0.38 (q=5), 0.27 (q=7), ~0.21 (q=10). Approximate scaling C ~ q^{-0.8} for q≥4. No simple analytic formula.

**Z_q charge detection fails for degenerate conjugate pairs.** Charges (k, q-k) in superposition average to misleading values. For q=10, sigma^2 (charges 2,8) falsely detected as charge 0. Must use independent R_epsilon for identification.

**POTENTIALLY NOVEL:** First C_{sigma*,sigma,epsilon} for Hermitian q>4 Potts. Peak at q=3 and monotonic decrease unreported in literature.

**Surprises:**
- C_sse peaks at q=3 (~0.54), NOT at q=2 (0.50)
- No simple analytic formula for C_sse(q) despite simple formulas for g_c(q) and c(q)
- Z_q charge superposition problem is universal: ALL conjugate pairs fail naive detection
- sigma^k–sigma cross-channel couplings are nonzero (Z maps charge k→k+1)

[Full report: sprints/sprint_060.md]

### Sprint 061 — Large-q Limit: Free Boson Convergence at q=15-30
**Status:** Complete (4 experiments).

**CFT data extended to q=30 at n=4.** Exact diag with periodic BC, dim up to 810k (q=30 n=4). Z_q charge resolution via symmetry generator G=X₁⊗...⊗Xₙ cleanly separates epsilon (charge 0) from sigma² (charge 2).

**Harmonic ratios converge to free boson k² as O(q^{-3}).** At q=30: R(σ²)/4=0.9985 (0.15% deficit), R(σ³)/9=0.9967, R(σ⁴)/16=0.9942, R(σ⁵)/25=0.9911. Convergence rate q^{-3.0} to q^{-3.2} depending on harmonic — much faster than naive O(1/q).

**x₁ does NOT saturate.** From n=4 descendant gap: x₁(15)≈0.071, x₁(20)≈0.055, x₁(25)≈0.045, x₁(30)≈0.038. Continues toward 0.

**Epsilon and sigma² diverge at large q.** R(σ²) → 4 (free boson). R(ε) grows as ~0.88·q at n=4. They are ALWAYS in different Z_q charge sectors.

**C_sse continues decreasing:** 0.093 (q=15), 0.063 (q=20), 0.047 (q=25), 0.038 (q=30). Power law ~q^{-1.2} at n=4.

**n=4,5 energy gap crossings unreliable for large q.** q=10 calibration shows 21% FSS correction needed (vs 5% for n=4,6). Odd-n sublattice effect. Formula g_c predictions used instead.

**POTENTIALLY NOVEL:** First charge-resolved conformal spectrum of Hermitian q-state Potts for q=15-30. The O(q^{-3}) convergence rate to free boson and the R_epsilon~q growth are new measurements.

**Surprises:**
- Harmonic ratio convergence is O(q^{-3}), not O(1/q) — no known prediction
- R_epsilon grows linearly with q at n=4 while R_sigma² saturates at 4
- x₁ continues decreasing to 0.038 at q=30, no saturation
- n=4,5 pairs need 21% FSS correction (vs 5% for n=4,6) — odd-n effect

[Full report: sprints/sprint_061.md]

### Sprint 062 — Compactification Radius & Twisted BC: c·x₁ ≈ 1/9 Discovery
**Status:** Complete (3 experiments).

**c·x₁ ≈ 1/9 for q≥3.** Product c(q)·x₁(q) = 0.112 ± 0.005 for all q=3-10, close to 1/9=0.111. q=2 Ising is an outlier (c·x₁=1/16). Model-independent constraint on the CFT. **POTENTIALLY NOVEL.**

**Single compact boson RULED OUT.** c>1 for q≥5 incompatible with c=1 for any single boson. Multiple DOF required. x₁ ≈ 0.206·q^(-0.449) (best power-law fit).

**q=2,3 are Luttinger liquids; q≥4 are NOT.** Twisted BC spin stiffness ρ_s·L converges for q=2,3 (drift <1.5%) but 10-22% for q≥5. ΔE(2)/ΔE(1) far from free boson prediction of 4: q=4→1.42, q=5→1.78, q=7→2.22.

**q=2 twist field = σ field EXACTLY.** ΔE(twist)/Δ₁ = 1.000 at all sizes. q=3: ratio ≈ 1.01.

**Twist ratios converge toward spectrum ratios from below.** 15-33% deficit at accessible sizes, shrinking with n.

**Surprises:**
- c·x₁ ≈ 1/9 for q≥3 — no theoretical prediction
- q=2 twist field has zero FSS correction (exact at n=4)
- ρ_s·L drift grows dramatically with q (1% → 22%)
- q=10 twist sector ratios at n=4 saturate near k=4

[Full report: sprints/sprint_062.md]

### Sprint 063 — Testing c·x₁: Independent c(q=15), Clock Model, Exact Analysis
**Status:** Complete (3 experiments).

**c(q=15) = 1.549 (DMRG n=8).** Gives c·x₁ = 0.110, within 1% of target 0.111. Pattern holds up to q=15. q=20 DMRG infeasible (>600s).

**c·x₁ is Potts-specific, NOT Z_q-universal.** Clock q=5 gives c·x₁ ≈ 0.146 (30% above Potts 0.112). The cos(2π·Δs/q) coupling creates different CFT data than δ(s_i,s_j).

**c·x₁ is NOT exactly 1/9.** Exact q=3 value is 8/75 = 0.10667 (4% below 1/9). Approximate constancy is a coincidence of the particular q-dependences.

**Clock g_c(q=5) = 0.52, not 0.67.** Energy gap method corrects old MI-CV estimate. Clock c(q=5) ≈ 1.17 > Potts 1.10.

**Surprises:**
- c·x₁ holds to 1% at q=15 despite being approximate
- Clock c > Potts c at same q — cos coupling is "more critical"
- Minimal model c·x₁ → 0 as m→∞ (opposite to Potts trend)

[Full report: sprints/sprint_063.md]

### Sprint 064 — Entanglement Hamiltonian at q>4 Potts: BW Locality Degrades Monotonically
**Status:** Complete (3 experiments).

**BW locality measured for q=4,5,7,10 Potts at true critical points.** Extends Sprints 032-034 to the novel CFT regime. At fixed n_A=4: q=2 (91%) > q=3 (76.5%) > q=4 (56%) > q=5 (42.3%) — monotonic decrease, ~15% drop per unit q.

**H/G-inv ratio predicts ordering perfectly.** Z_q-invariant operator fraction = exactly 1/q — a clean algebraic identity. H/G-inv scales as ~5/q^5, plummeting because symmetry-allowed operators grow as q^{2n_A-1} while physical Hamiltonian terms are fixed.

**Unruh-like gradient persists for ALL q.** Bond and site couplings in H_E decrease monotonically toward entanglement cut. But this captures only 36-50% of H_E at q=4,5 (vs 76-91% at q=2,3).

**Peak BW locality shifts into ordered phase for q≥4.** At q=2,3 peak is near g_c; at q≥4 peak is at g≈0.30 regardless of g_c value.

**Surprises:**
- G-inv fraction = exactly 1/q — clean formula, previously uncomputed
- Linear envelope wins for q≥4 (sin_inv for q≤3) — envelope character changes at q=4 boundary
- ~15% drop per unit q suggests BW → 0% near q~8 at n_A=4
- BW degradation correlated with c(q)~ln(q): richer CFTs have less local H_E

[Full report: sprints/sprint_064.md]

### Sprint 065 — Hybrid vs Clock Universality: DIFFERENT Classes Confirmed
**Status:** Complete (3 experiments).

**Hybrid ≠ clock universality — confirmed by THREE independent probes at q=5.** c/x₁: hybrid 10.77 vs clock 9.43 (12% diff, shrinking but nonzero). ν: hybrid 0.83 (finite, power-law) vs clock ~2+ (diverging → BKT). Clock q=7 has NO gap crossing (BKT); hybrid has clear crossing.

**Clock q=5 has BKT transition.** Gap slopes barely grow: 2.29→2.67→2.95 for n=4,6,8. Hybrid slopes clearly grow: 3.99→6.16→8.43. BKT confirmed as expected from Sun et al. 2020.

**Potts-clock hybrid defines a NEW universality class for q≥5.** Not first-order (like S_q Potts), not BKT (like Z_q clock). Continuous, second-order, power-law transition with finite ν ≈ 0.83. **POTENTIALLY NOVEL.**

**Clock DMRG ~10x slower than hybrid** due to dense cos coupling (q² terms vs q terms in MPO). Limits clock characterization at large q.

**Surprises:**
- ν comparison is the sharpest discriminator: qualitatively different (finite vs diverging)
- Clock q=7 gap monotonically increasing — no crossing possible
- c/x₁ difference slowly shrinks (14%→12%) suggesting SOME shared structure
- Clock q=7 g_c ≈ 0.30 (if applicable), far below hybrid 0.535

[Full report: sprints/sprint_065.md]

### Sprint 066 — Weakly First-Order Test: No First-Order Signal at q=10
**Status:** Complete (3 experiments).

**No first-order signal detected at accessible sizes (n≤6).** Three diagnostics all consistent with continuous transition:

**Δ·N at g_c DECELERATES (066a).** q=10 at n=4,5,6: Δ₁·N = 0.242, 0.250, 0.255. Changes: +3.1%, +2.2%. Decelerating = convergent, not divergent. First-order would accelerate.

**Gap minimum is SIZE-INDEPENDENT (066b).** Δ_min ≈ 0.030 for n=4,5,6. Not exponentially shrinking (first-order signature). Pseudo-critical point drifts from g=0.50→0.60 toward g_c=0.684.

**Cross-q comparison (066c).** q=3 (known continuous): Δ·N decreases from 0.746 to 0.722 (n=4→12). q≥5: Δ·N increases (anomalous). Normalized change at n=6 saturates at ~5.5% for q≥7. Sign-flipped corrections, not first-order.

**Closes the last major open question.** The hybrid model defines a genuinely continuous universality class for q=2-10, distinct from first-order (S_q) and BKT (Z_q clock).

**Δ₂/Δ₁ = 1.0000 exactly** at all sizes — Z₁₀ conjugate pair degeneracy perfectly maintained.

**Surprises:**
- Gap minimum barely changes with N (0.0304, 0.0300, 0.0308) — no avoided crossing
- q=3 from above, q≥5 from below — qualitative sign flip in FSS corrections
- Anomaly saturates between q=7 and q=10 (both ~5.5%)
- Pseudo-critical point drift is linear in n, not 1/N^{1/ν}

[Full report: sprints/sprint_066.md]

### Sprint 067 — Two-Transition Scan: Hybrid Has No Floating Phase
**Status:** Complete (3 experiments).

**Hybrid model has exactly ONE transition — no floating phase.** Three independent probes confirm:

**Wide gap scan (067a).** Scanned g=[0.02, 3.0] at q=5 (n=4,6,8) and q=7 (n=4,6). No second gap minimum. Gap increases monotonically from g_c. Gap*N → 2·sin(2π/q) at large g.

**Entropy scan (067b).** DMRG n=16,24 at q=5. c_eff=1.22 at g_c=0.44, drops to 0.003 at g=0.55. No c≈1 region above g_c. Sharp critical-to-gapped transition.

**Clock comparison (067c).** Same method on clock model: **floating phase clearly detected** at g=[0.30, 0.92] (Δg=0.62). Hybrid gapless only at g=[0.38, 0.50] (Δg=0.12, just g_c). Validates method — absence in hybrid is real, not a size limitation.

**Fourth difference between hybrid and clock universality:** clock has floating phase, hybrid does not. The Potts δ-function coupling suppresses the intermediate Luttinger-liquid phase that the cosine coupling allows.

**Surprises:**
- Clock floating phase is WIDE (Δg≈0.62) — easily detectable even at n=4,6
- Hybrid c_eff drops from 1.22 to 0.003 between g=0.44 and g=0.55 (sharp!)
- No local gap minima in either model — floating phase creates flat region, not a second dip
- Large-g gap*N = 2·sin(2π/q): exact formula for paramagnetic limit
- Clock n=8 scan 7x slower than hybrid (886s vs 128s) due to dense coupling matrix

[Full report: sprints/sprint_067.md]

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
