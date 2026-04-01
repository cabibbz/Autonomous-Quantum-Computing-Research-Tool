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

### Sprint 054 — Central Charge c(q) at True Critical Points
**Status:** Complete (5 experiments).

**c(q) extracted from entropy scaling S=(c/6)ln(n)+const.** Method validated: c(q=2) pairwise converges 0.544→0.516 toward exact 0.500. c(q=3) pairwise converges 0.934→0.884 toward exact 0.800. Chi convergence proven: entropy identical at chi=20 and chi=80 for q=3 n=16. DMRG matches exact diag to machine precision (ΔE=2e-14, ΔS=0 at q=3 n=8).

**c(q=4) ≈ 1.23 raw, FLAT — not converging at n≤24.** Pairwise c: 1.231, 1.222, 1.229 — barely varies. Overshoot +23% vs CFT c=1. Qualitatively different from q=2,3 which converge clearly. Consistent with logarithmic corrections at marginal q=4 (Ashkin-Teller).

**c(q=5) ≈ 1.37 raw, converging downward.** Pairwise c: 1.395→1.335. Even with ~20% overshoot correction (calibrated from q=2,3), c(q=5) ≈ 1.1-1.2 — **above c=1, outside minimal model series.** No CFT prediction exists. **POTENTIALLY NOVEL.**

**c(q) increases monotonically with q:** 0.5 → 0.8 → 1.0 → (>1). Central charge grows as Potts state space grows.

**Surprises:**
- c(q=4) overshoot is FLAT — logarithmic corrections prevent convergence at n≤24
- Chi convergence trivial (chi=20 suffices) — NOT the overshoot source
- c(q=5) > 1 even after generous correction — outside Potts minimal models
- Overshoot ratio increases with q: 9%→11%→23% at (n=16,24)

[Full report: sprints/sprint_054.md]

### Sprint 055 — Entropy Profile Method & c(q=5) > 1 Confirmed
**Status:** Complete (7 experiments, 1 abandoned).

**iDMRG fails at criticality.** S vs ln(xi) with infinite MPS gives c=0.41 for TFIM (18% error). Correlation length saturates at xi~370 even at chi=160. Pairwise c values scatter wildly. **Do not use iDMRG for c extraction.**

**Entropy profile method validated and superior.** S(l) vs chord distance ln[(2n/pi)sin(pi*l/n)] at single large n. TFIM: c = 0.524→0.512 at n=32→64 (exact 0.500). q=3: c = 0.870→0.827 at n=16→48 (exact 0.800). 2-5x lower overshoot than FSS at same sizes.

**c(q=5) > 1 CONFIRMED by two independent methods.** Profile: c=1.261 (n=16). FSS (Sprint 054): c=1.335 (n=12,16). Both overshoot-corrected to c ≈ 1.10 ± 0.10. **Even with 25% correction, c > 1.0.**

**c(q=4) flat overshoot confirmed.** Profile: c=1.144 (n=16), 1.148 (n=24) — barely changes. Independent confirmation of logarithmic corrections.

**Surprises:**
- iDMRG correlation length saturates at ~370 at chi=160 — L=2 unit cell fundamentally limited
- Profile overshoot at n=16 grows strongly with q: 8.7% (q=3) → 14.4% (q=4) → ~20% (q=5)
- q=7 DMRG at n=12 exceeds 120s — d=7 infeasible for profile method
- Even/odd oscillations negligible in entropy profile (<0.001)

[Full report: sprints/sprint_055.md]

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
