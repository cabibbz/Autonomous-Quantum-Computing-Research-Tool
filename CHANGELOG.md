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

### Sprint 049 — Correlation Length & Entropy FSS: q=3 Potts Critical Point WRONG
**Status:** Complete (6 experiments).

**TFIM central charge c validated.** Entropy S(g_c=1.0, n) at n=16-64 gives c = 0.529, with pairwise estimates converging: 0.544 → 0.532 → 0.523 → 0.516 toward exact c=0.500.

**Correlation length from correlator decay works but is finite-size limited.** ξ matches exact ξ=1/ln(g) far from g_c (5% at g=2) but saturates at ξ~n/4. Naive ν extraction gives 0.62 (n=40), 0.72 (n=80) — well below ν=1.

**Entropy FSS does NOT give ν.** Logarithmic singularity S ~ ln(ξ) defeats power-law collapse. Best collapse at ν=3, not ν=1.

**q=3 Potts critical point at g_c ≈ 0.33, NOT 1.0.** Exact diag confirms S=0.047 at g=1.0 (deep disordered phase). True phase transition between g≈0.25-0.35 where entropy grows logarithmically with n and central charge estimate crosses c=4/5.

**10 sprints of Potts results (038-048) were at wrong g.** MI-CV crossings near g≈1.0 are a disordered-phase crossover. All Potts g_c values, the g_c scaling law, and ν(q) results need re-examination.

**Surprises:**
- q=3 at g=1.0 has S=0.047 — confirmed by exact diag, not a DMRG failure
- GHZ entropy S=ln(3) persists to g≈0.20 in ordered phase
- TFIChain Z₂ conservation is essential — PottsChain has NO conservation
- Most expensive mistake: 10 sprints of Potts physics at the wrong critical point

[Full report: sprints/sprint_049.md]

### Sprint 050 — True Potts Critical Points via Self-Duality & MI-CV Vindication
**Status:** Complete (3 experiments).

**Self-duality gives EXACT g_c for q=2,3.** Our Hamiltonian H = -Jδ - g(X+X†) has Kramers-Wannier self-duality when X+X† spans all non-trivial generators of Z_q. This holds for q=2 (X+X†=2X, all 1 generator) and q=3 (X+X†=X+X²=Γ, all 2 generators). For q≥4, X+X† misses intermediate powers → self-duality broken.
- q=2: g_c = J/4 = 0.25 (different convention from TFIChain g_c=1.0)
- q=3: g_c = J/3 = 0.333 (EXACT)

**MI-CV crossings CONFIRMED at true q=3 g_c.** n=8 and n=12 MI-CV cross at g≈0.26 (finite-size shifted below g_c=1/3). Ordered phase: CV(n=12)<CV(n=8). Disordered: reversed. The crossing signature from Sprints 038-039 is qualitatively correct.

**q=4 Potts: g_c ≈ 0.34 from pseudo-critical (n=8).** No self-dual prediction. DMRG at n=12 too slow (~120-230s/point at χ=40). Central charge c(6,8) peaks at 3.31 — massive overshoot, same as q=3.

**Surprises:**
- Self-duality breaks at q≥4 — fundamental asymmetry of the X+X† transverse field
- Our Potts q=2 has g_c=0.25, NOT 1.0 — different normalization from TFIChain
- MI-CV crossings survive at true g_c — prior qualitative conclusions rescued
- Central charge extraction from n≤12 overshoots by 3× for q=3 Potts (c=2.56 vs 0.8)
- q=4 pseudo-critical drifts DOWNWARD with n (0.35→0.34), opposite to q=3

[Full report: sprints/sprint_050.md]

### Sprint 051 — Energy Gap Method: g_c INCREASES with q, Reversing Old Picture
**Status:** Complete (5 experiments).

**Energy gap Δ·N crossing validated and applied.** At criticality, Δ = E₁-E₀ ∝ 1/N, so Δ·N is scale-invariant at g_c. Validated on q=2 (crossings→0.246 vs exact 0.250) and q=3 (crossings→0.325 vs exact 0.333). Systematic convergence from below, ~2-3% error.

**g_c(q) INCREASES with q:**
| q | g_c |
|---|-----|
| 2 | 0.250 (exact) |
| 3 | 0.333 (exact) |
| 4 | ~0.39 |
| 5 | ~0.44 |
| 7 | ~0.52 |

**Physical mechanism:** q-fold ground-state degeneracy requires stronger transverse field X+X† to destroy order. X+X† only creates nearest-state transitions, so mixing time grows with q.

**Reverses 10 sprints of wrong results.** Old trend (g_c decreasing) was from MI-CV disordered-phase crossovers. q=5 old estimate (0.41) was closest to truth (0.44) — accidentally near true g_c.

**Surprises:**
- g_c INCREASES with q — opposite to invalidated scaling law
- q=7 g_c ≈ 0.52, not 0.26 — off by factor of 2
- Energy gap converges faster than entropy or MI-CV (~2% error vs ~30%)
- q=5 old ν≈2.0 may be approximately correct (measured near true g_c)

[Full report: sprints/sprint_051.md]

### Sprint 052 — g_c(q) Scaling Law: √(q-1) Formula + q=10 Verification
**Status:** Complete (4 experiments).

**g_c(q) ≈ (1/5)√(q-1) + 1/20.** Best 2-parameter fit to 6 data points (q=2-10), χ²/dof=0.40. Coefficients are simple fractions (a=1/5, b=1/20). Gives exact q=2, <5% error for all tested q. **POTENTIALLY NOVEL.**

**q=10 measured for first time: g_c ≈ 0.684.** Energy gap Δ·N crossing at n=4,6 gives raw 0.652, corrected to 0.684 with 4.8% FSS correction. 5-point prediction (0.616) was 10% too low — logarithmic growth too slow, √(q-1) is correct.

**FSS corrections recalibrated by size pair.** n=4,6 pairs need 4.8% correction (not 2.5%). q=7 corrected upward from 0.524 to 0.535.

**Surprises:**
- g_c ∝ √(q-1) — simplest formula wins, free exponent 0.516 ≈ 1/2
- Integer-coefficient formula (4√(q-1)+1)/20 fits all 6 points
- g_c(q) diverges — no saturation at large q
- 5-point prediction underestimated q=10 by 10%

[Full report: sprints/sprint_052.md]

### Sprint 053 — ν(q) at True Critical Points: Method Validated, Old ν Values WRONG
**Status:** Complete (6 experiments).

**Corrected energy gap slope method validated.** d(Δ·N)/dg ~ N^{1/ν}·(1+b/N) with b=0.86 calibrated from q=3 (ν=5/6 exact). Achieves <1% error for q=2 (ν=1.00), 3% for q=3 (ν=0.86 vs 5/6). Data collapse is WORST method (43% error at n≤10).

**ν(q=3,4,5) ≈ 0.82-0.86, nearly constant.** q=2: 1.00 (exact 1.0). q=3: 0.86 (exact 5/6). q=4: 0.82 (3 size pairs all agree). q=5: 0.85 (1 pair). q=7: 0.97 (1 pair, low confidence). q=10: 1.12 (1 pair, suspect).

**ALL old ν values for q≥4 WRONG.** Old q=4 ν≥2.2, q=5 ν≈2.0, q=7 ν≈0.5 — all artifacts of MI-CV data collapse at wrong critical points. The "BKT peak" at q=4 and "non-monotonic ν(q)" were entirely artifactual.

**Surprises:**
- ν(q=4) = 0.82, NOT ≥2.2 — no BKT behavior
- ν(q=5) = 0.85, NOT 2.0 — MI-CV collapse overestimates ν by ~2.5×
- Data collapse is worst ν method (43% error) — corrected power-law is best (3%)
- ν nearly constant ~0.85 for q=3,4,5 — near-universal exponent

[Full report: sprints/sprint_053.md]

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
