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

Full details for all compressed sprints are in sprints/sprint_NNN.md.

---

## Detailed Sprint Log (Recent)

### Sprint 025 — Real Hardware QPU Test
**Status:** Complete. **QPU: 20s on ibm_kingston (Heron, 156q).**

[[5,1,3]] isotropy confirmed on real hardware: asymmetry 0.040 (6.4x more isotropic than 3-qubit's 0.254). Hardware beats simulator noise model. First QEC advantage: 3-qubit Z-basis 0.976 vs uncoded 0.959. Correlated noise degrades isotropy 4x vs simulator but doesn't destroy it.

[Full report: sprints/sprint_025.md]

### Sprint 026 — Active Syndrome Extraction: Non-FT Always Hurts
**Status:** Complete.

Active correction worse than passive at ALL error rates (p2q=0.0001 to 0.02). Root cause: CX gates propagate ancilla errors to weight-2 data errors, exceeding code's correction capacity.

[Full report: sprints/sprint_026.md]

### Sprint 027 — Flag-FT Syndrome: Solves the Wrong Problem
**Status:** Complete.

Flag qubits detect all weight-2 propagation (8/8, 0 false flags). But flag-FT still never beats passive — single-round extraction fundamentally insufficient. Threshold theorem requires FT gadgets + repeated measurement + time-aware decoding.

[Full report: sprints/sprint_027.md]

### Sprint 028 — Repeated Syndrome: Gate Overhead Kills
**Status:** Complete.

3-round majority improves syndrome quality but adds 48 2Q gates. Active correction NEVER beats passive at ANY noise level or round count. QEC active correction arc complete: [[5,1,3]] is too small.

[Full report: sprints/sprint_028.md]

### Sprint 029 — TFIM Phase Transition: Fifth Archetype
**Status:** Complete.

Ordered TFIM = GHZ. Critical point = new Scale-Free archetype (equidistant from all four prior archetypes, power-law correlations). MI uniformity CV is a new order parameter (0.03 → 0.39 → 1.16).

[Full report: sprints/sprint_029.md]

### Sprint 030 — XXZ Model: Archetype Loop & MI-CV Transition Classifier
**Status:** Complete.

XXZ traces a loop in archetype space: Democratic→Scale-Free→Geometric→Scale-Free→Democratic. MI-CV shape classifies transition type: jump=first-order, inflection=Ising, dome=BKT. Archetype boundaries ≠ phase boundaries (I3 sign change at Δ≈0.7 inside XY phase).

[Full report: sprints/sprint_030.md]

### Sprint 031 — Entanglement Spectrum: Three Levels of Description
**Status:** Complete.

Three complementary levels: scalar (entropy=amount), correlation (MI/I3=topology), spectral (eigenvalues=symmetry). CFT doublet structure at TFIM criticality. GHZ and Cluster 1D spectrally identical at half-cut. I3 sign change invisible in spectrum — orthogonal information.

[Full report: sprints/sprint_031.md]

### Sprint 032 — Bisognano-Wichmann: Symmetry Controls BW Accuracy
**Status:** Complete.

TFIM (Z₂): 91% BW locality. XXZ (U(1)): 100% BW locality. BW peaks in gapped phases, NOT at criticality. Entanglement temperature: Unruh-like gradient (hot at cut, cold in bulk). Fourth level of entanglement description: Hamiltonian (locality + temperature).

[Full report: sprints/sprint_032.md]

### Sprint 033 — Potts BW: Local Dimension Matters More Than Group Size
**Status:** Complete.

S₃ (order 6, d=3): 76.5% BW locality — LOWER than Z₂ (91%, d=2). Overturns "bigger group → better BW" prediction. BW accuracy depends on BOTH G and d. Potts ordered phase = GHZ-3 with triplet spectrum.

[Full report: sprints/sprint_033.md]

### Sprint 034 — BW Operator Algebra: H/G-inv Ratio as Predictor
**Status:** Complete.

Explicit operator counting: U(1)=0.143 > Z₂=0.099 > S₃=0.072 > Z₃=0.040 (H/G-inv ratio). Perfect monotonic correlation with BW locality: 100% > 91% > 76.5% > 69.7%. Z₃ clock ≡ S₃ Potts at q=3. Chiral clock genuinely breaks S₃→Z₃, BW drops to 69.7%.

[Full report: sprints/sprint_034.md]

### Sprint 035 — BW Size Scaling: Pauli Fraction Fails, Spectrum Succeeds
**Status:** Complete (3 experiments).

**CRITICAL FINDING: Sprint 032's "91% TFIM locality" is a finite-size illusion.** Same Pauli fraction metric gives 47% at n=12 and 8% at n=20. The metric breaks because ||H_E||² grows exponentially while TFIM operator count grows linearly.

**But BW FORM is correct at spectrum level.** Entanglement spectrum R²=1.0 at ALL tested sizes (n=8→16) when TFIM-type coefficients are fitted independently. The physically relevant eigenvalues live in the TFIM operator subspace. Overfitting caveat: only 3-4 significant eigenvalues vs 7-15 fitting parameters.

**The BW envelope is wrong.** Standard sin_inv envelope gives R²≈0.15 across all sizes. The β(i) profile at finite size with open boundaries is far from the CFT prediction.

**Full density matrix reconstruction fails catastrophically.** Projected entropy: 0.76 vs 0.52 exact at n=8, blowing up to 9.7 vs 0.64 at n=20 (1400% error). Non-TFIM terms collectively encode low-rank structure.

**MI-CV sharpens with system size.** TFIM transition slope increases (n=8→16). Ordered phase stable (CV≈0.05). XXZ BKT dome narrows.

**Surprises:**
- Spectrum R²=1.0 while Pauli fraction=8% — metrics disagree by 12×
- S_projected/S_exact diverges from 1.5× to 15× with system size
- BW envelope is far more wrong than BW form

[Full report: sprints/sprint_035.md]

### Sprint 036 — MI-CV Finite-Size Scaling: Genuine Divergence Confirmed
**Status:** Complete (3 experiments).

**MI-CV is a genuine order parameter.** Transition slope at TFIM h/J=1.0 diverges as ~n^1.1: slope 1.72 (n=8) → 3.95 (n=16) → ~7.8 (n=32). Ordered phase CV → 0, disordered CV → ∞ with system size. Crossing point at h/J ≈ 0.93 (7% finite-size shift).

**XXZ BKT dome narrows with size.** CV growth ratio at Δ=1.0 is 1.34 vs 1.17 at Δ=1.5 (n=16→32). Néel phase (Δ=2) shows CV convergence: 1.13 → 0.78. BKT boundary between diverging/converging CV sharpens toward Δ=1.

**New technique: correlation-function MI reconstruction.** All-pairs MI from 9 Pauli correlation matrices computed O(n²·χ²) from MPS. Enables n=50 all-pairs MI. Validated exact at n=8 (diff=0.000000).

**Surprises:**
- NN-only MI-CV is useless — boundary effects dominate, obscure phase transition
- Slope ~ n^1.1 suggests logarithmic corrections (expected for Ising universality)
- XXZ Néel CV drops to 0.78 at n=32 — ordered phases become increasingly Democratic

[Full report: sprints/sprint_036.md]

### Sprint 037 — Critical Exponent & First-Order Transition: MI-CV Classification Complete
**Status:** Complete (5 experiments).

**TFIM critical exponent from crossing points.** Four crossings at n=8-32: h_c shifts from 0.895 to 0.972. Power-law fit gives ν=0.755, but two-point estimates converge toward Ising ν=1 (0.74→0.88). Corrections to scaling dominate at n≤32.

**XXZ first-order transition at Δ=-1.** FM phase has CV=0.000 exactly. Step-function onset at Δ=-1 with jump growing as ~n^1.0. Gradient scales ~n^0.9. No curve crossings (unlike Ising).

**MI-CV transition classification complete.** Three transition types have three distinct signatures: crossing curves (Ising), step function (first-order), smooth dome (BKT). Presence/absence of crossings is the cleanest first-order diagnostic.

**Surprises:**
- FM phase CV is exactly zero, not just small — purest Democratic state found
- First-order gradient (~n^0.9) SLOWER than Ising (~n^1.1) — first-order sharpens width, not slope
- No crossing at first-order is a qualitative signature, simpler than quantitative gradient analysis

[Full report: sprints/sprint_037.md]

### Sprint 038 — Data Collapse & Potts MI-CV
**Status:** Complete (4 experiments).

**Data collapse confirms Ising universality.** Consolidated TFIM MI-CV data (n=8-50) and tested scaling collapse CV(h,n) = F((h-h_c)·n^{1/ν}). Optimal ν converges toward Ising exact (ν=1): all sizes give ν=0.80, large sizes (n≥16) give ν=1.04, largest pair (n≥24) gives ν=1.12. The finite-size fit ν=0.755 (Sprint 037) is an artifact — it gives WORSE collapse than ν=1.

**q=3 Potts MI-CV shows crossing curves.** Extended all-pairs MI to qutrits using Gell-Mann basis (8 generators). Potts transition at g/J=1 shows same crossing signature as TFIM: n=8 CV > n=12 CV in ordered phase, reversed in disordered phase. Crossing at g≈0.85-0.90. Potts transition is 14% steeper than TFIM at same n.

**Surprises:**
- ν=0.755 gives WORSE collapse than ν=1 — crossing-point ν extraction is fragile
- Quality landscape extremely flat near ν=1 (0.3% from optimal) — explains why crossing points gave misleading ν
- Potts CV is systematically higher than TFIM above transition — larger d amplifies MI inhomogeneity

[Full report: sprints/sprint_038.md]

### Sprint 039 — Potts Data Collapse: MI-CV Distinguishes Universality Classes
**Status:** Complete (3 experiments).

**Potts n=16 MI-CV confirms crossing at three sizes.** CV decreases with n in ordered phase (g=0.8: 0.341→0.271→0.236) and increases in disordered phase (g=1.0: 0.826→0.961→1.051). Crossing point drifts: g_c=0.923 (n=8,12) → 0.955 (n=12,16), converging toward self-dual g=1.0.

**Data collapse distinguishes Potts from Ising.** Joint optimization of (ν, g_c): Potts ν=5/6 gives 14% better collapse than Ising ν=1 (quality ratio 0.862 at n=12-16). Optimal ν converges from 0.81 (all sizes) to 0.87 (large sizes), approaching 5/6≈0.833. CAUTION: Fixed g_c=1.0 gives wrong answer (Ising wins) due to finite-size shift.

**Slope exponent provides second discriminator.** Potts slope ~ n^1.36 vs TFIM slope ~ n^1.1, consistent with 1/ν scaling (Potts 1/ν=1.2 > Ising 1/ν=1.0).

**Surprises:**
- Fixed g_c=1.0 gives WRONG universality class — joint optimization essential
- Slope exponent 1.36 overshoots 1/ν=1.2 — same finite-size inflation as TFIM
- n=24 qutrit DMRG without symmetry exceeds 60s/point — needs Z₃ conserve

[Full report: sprints/sprint_039.md]

### Sprint 040 — q=4 Potts MI-CV: Marginal Transition Shows Crossing, Not Dome
**Status:** Complete (3 experiments).

**q=4 Potts shows crossing curves, not dome.** Despite q=4 being the marginal boundary (2D Potts goes first-order at q>4), MI-CV at n=8,12 shows the same crossing signature as q=3 and TFIM. Crossing at g_c≈0.893, further below self-dual than q=3 (0.923).

**q=4 slope LOWER than q=3.** At g=1.0: q=4 slope 1.72 (n=8), 2.71 (n=12) vs q=3 slope 2.27, 3.98. Logarithmic corrections at q=4 suppress finite-size slope.

**q=4 CV systematically lower than q=3 above transition.** Ratio q4/q3 ≈ 0.87-0.95 in disordered phase. Larger d distributes correlations more evenly.

**Surprises:**
- Marginal q=4 does NOT show dome — crossings dominate at n≤12
- q=4 slope LOWER than q=3 — logarithmic corrections suppress slope
- Larger d gives MORE uniform MI (lower CV) above transition

[Full report: sprints/sprint_040.md]

### Sprint 042 — True q=5 Potts MI-CV: Second-Order, Not First-Order
**Status:** Complete (5 experiments).

**Custom Potts model built with Kronecker-delta coupling.** True q=5 Potts shows CROSSING curves at g_c≈0.41, disproving the first-order prediction. Potts transition 5.7x steeper than clock (slope 4.5-6.1 vs 0.8-1.0) but same qualitative second-order signature. The 2D classical "q>4 → first-order" rule does NOT apply to 1D quantum Potts with transverse field.

[Full report: sprints/sprint_042.md]

### Sprint 041 — q=5 Clock MI-CV: Crossings Persist, Shifted Dramatically
**Status:** Complete (4 experiments).

**q=5 clock model STILL shows crossing curves.** All three predictions (BKT dome, first-order step, no crossings) were wrong. MI-CV at n=8,12 shows crossings at g_c≈0.673, confirmed by interpolation between g=0.5 (CV decreases with n) and g=0.8 (CV increases).

**Crossing point shifts 10x more than q=3→4.** g_c trend: 0.93 (q=2) → 0.923 (q=3) → 0.893 (q=4) → 0.673 (q=5). The Δg=0.220 jump from q=4→5 is an order of magnitude larger than q=3→4 (0.030).

**Slope halves from q=4 to q=5.** At g=1.0 n=8: q=5 slope 0.86 vs q=4 slope 1.72. Consistent with gentler BKT-like transition in clock model.

**Clock ≠ Potts for q≥4.** ClockChain uses cos(2π(s_i-s_j)/q) coupling, which equals Potts only for q=2,3. True first-order test requires custom Potts model with Kronecker-delta coupling.

**Surprises:**
- ALL three predictions wrong — crossings persist at q=5
- g_c shift is 10x larger from q=4→5 than q=3→4
- Slope exactly halves (0.86 vs 1.72) — suspiciously clean
- Larger d continues to give lower CV above transition

[Full report: sprints/sprint_041.md]

### Sprint 043 — First-Order Search: q=10, q=20 Potts — All Second-Order
**Status:** Complete (4 experiments).

**New technique: direct MPS contraction for ρ_ij.** Replaces Gell-Mann correlation reconstruction for d≥10. Computes all-pairs MI in 1.4s at d=10, n=8 (vs estimated minutes for 9801 Gell-Mann correlations). Enables MI-CV at arbitrarily large local dimension.

**q=10 Potts: CROSSING confirmed at g_c≈0.246.** Clear MI-CV crossing between n=8 and n=12: CV increases with n for g≤0.20 (disordered), decreases for g≥0.25 (ordered). NOT first-order.

**q=20 Potts: universal large-q regime.** At χ=10, q=20 ground states are identical to q=10 (energies match to 4+ significant figures). Half-chain entropy rises continuously (0→0.95→1.04), no discontinuity. The relevant excitations live in {|0⟩, |1⟩, |q-1⟩} subspace for all q≥10.

**Conclusion: 1D quantum Potts with transverse field is second-order for ALL tested q (2-20).** g_c decreases monotonically: 1.0 (q=2) → 1.0 (q=3) → 0.89 (q=4) → 0.41 (q=5) → 0.25 (q=10). The 2D classical "q>4 → first-order" rule does not apply in 1D quantum.

**Surprises:**
- q=20 ground states IDENTICAL to q=10 at χ=10 — universal large-q regime
- g_c shift rate slows: Δg=0.48 (q=4→5) vs Δg=0.16 (q=5→10) — diminishing shift per q
- Direct MPS contraction removes the d²-bottleneck entirely

[Full report: sprints/sprint_043.md]

### Sprint 044 — g_c Scaling Law: Predictive Formula for Potts Critical Points
**Status:** Complete (4 experiments).

**g_c scaling law extracted and verified.** Fitted g_c(q) for 1D quantum Potts to 6 candidate forms. Best fit: g_c ≈ 0.87*(q-3)^(-0.85) for q≥4, with g_c=1 for q=2,3 (self-duality protected). The pole at q=3 reflects the breaking of self-duality. Exponent 0.85 is tantalizingly close to q=3 Potts ν=5/6.

**Blind prediction validated at q=7.** Power law (fitted without q=7) predicted g_c(7)=0.263 vs measured 0.259 — only 1.6% error. MI-CV crossing confirmed at g_c≈0.259 (n=8,12, χ=20). Second-order transition, consistent with all other q.

**χ convergence matters for large d.** χ=10 gives 25% CV inflation at d=7. Must use χ≥20 for d≥7.

**Surprises:**
- Blind prediction accurate to 1.6% — scaling law is genuinely predictive
- Best exponent 0.85 ≈ 5/6 — possible deep connection to Potts universality
- g_c(7)=0.259 vs g_c(10)=0.246 — nearly flat, large-q regime starts at q≈7
- The "cliff" is between q=3 and q=5, driven by self-duality breaking at q=3

[Full report: sprints/sprint_044.md]

### Sprint 045 — ν(q) Extraction: ν(q=5)≈2, Universality Class Changes with q
**Status:** Complete (6 experiments).

**ν(q=5) ≈ 2.0 from data collapse.** Joint (ν, g_c) optimization across n=8,12,16 (all direct MPS contraction, χ=20). Free fit: ν=2.24, g_c=0.446. ν=2 nearly as good (ratio 1.04). Independently confirmed: slope~n^0.49 → ν≈2.05. Ising (ν=1) is 3.7x worse, Potts q=3 (ν=5/6) is 5.8x worse.

**ν(q) is NON-MONOTONIC.** q=2: ν=1, q=3: ν=5/6, q=5: ν≈2. Decreases q=2→3, then increases sharply. Large ν at q=5 means transition is "almost first-order" — wide crossover, enormous finite-size effects.

**Sprint 042d Gell-Mann MI was WRONG for q=5.** Errors up to 11x (g=0.50: CV=1.36 → 0.12). Gell-Mann correlation reconstruction unreliable for d≥5. Direct MPS contraction (Sprint 043) is essential.

**Surprises:**
- ν=2 at q=5 — opposite to mean-field prediction of ν→1/2
- Non-monotonic ν(q): decreases then increases
- Gell-Mann method fails catastrophically at q=5
- Finite-size crossing points 30% below thermodynamic g_c (large ν amplifies finite-size effects)

[Full report: sprints/sprint_045.md]

### Sprint 046 — ν(q=4) Extraction: Marginal Point Shows BKT-like Behavior
**Status:** Complete (6 experiments).

**ν(q=4) ≥ 2.2 and trending upward with system size.** Slope ratio near g_c≈0.89: ν=1.92 (n=8,12) → 2.71 (n=12,16). This INCREASES with n — signature of BKT or logarithmic corrections at the marginal q=4 point.

**Standard FSS collapse FAILS at q=4.** Constrained at g_c=0.89, optimizer hits ν→∞ with poor quality. Free collapse finds g_c≈0.53 (fitting CV minimum, not phase transition). MI-CV curves don't cross near g_c — monotonically ordered n=16>n=12>n=8 at all g≥0.50.

**Sprint 040 Gell-Mann data at d=4 was misleading.** Gell-Mann showed crossings at g_c≈0.893 that don't exist with direct MPS. Gell-Mann now unreliable for d≥4 (previously claimed d≤4 safe, but only validated at d=2,3).

**ν(q) updated: q=4 is the marginal/BKT peak.**
| q | ν | Nature |
|---|---|--------|
| 2 | 1.0 | Standard 2nd-order |
| 3 | 5/6 | Standard 2nd-order (minimum) |
| 4 | ≥2.2 (diverging?) | Marginal/BKT-like |
| 5 | ~2.0 | Large-ν 2nd-order |

**Surprises:**
- ν estimate INCREASES with n (1.92→2.71) — first BKT signature in the q series
- NO MI-CV crossings near g_c — qualitatively different from all other tested q
- Gell-Mann misleading at d=4, not just d≥5
- CV minimum at g≈0.4-0.5 far below g_c — ordered phase has complex MI structure

[Full report: sprints/sprint_046.md]

### Sprint 047 — ν(q=7) Extraction: Crossings Return, q=4 Confirmed as BKT Peak
**Status:** Complete (4 experiments).

**q=7 MI-CV shows CROSSING CURVES.** n=8 and n=12 MI-CV curves cross at g_c≈0.244 (6% below scaling law prediction of 0.259). Ordered phase: n=12 CV < n=8 CV (e.g. 0.421 vs 0.597 at g=0.15). Disordered phase: n=12 CV > n=8 CV (e.g. 0.541 vs 0.438 at g=0.30). This is qualitatively identical to q=2,3,5 and qualitatively DIFFERENT from q=4 (no crossings).

**ν(q=7) ≈ 0.5 (near mean-field).** Disordered-side slope ratio [0.26,0.30] gives ν≈0.49. Uncertain from only 2 sizes at coarse g-resolution, but clearly much smaller than q=5 (2.0) or q=4 (≥2.2).

**q=4 confirmed as unique BKT peak.** It is the ONLY tested q (out of 2,3,4,5,7) that lacks MI-CV crossings. ν(q) has a sharp peak at q=4: 1.0 → 5/6 → ≥2.2 → 2.0 → 0.5.

**Surprises:**
- Crossings RETURN at q=7 — q=4 is uniquely crossing-less
- ν drops 4x from q=5 (2.0) to q=7 (0.5) — approaches mean-field
- n=8 has sharp CV kink at g_c that smooths at n=12
- MI-CV crossing g_c=0.244, 6% below predicted g_c=0.259

[Full report: sprints/sprint_047.md]

### Sprint 048 — ν(q=10): Dead-Pair Bias Invalidates Raw MI-CV at Large d
**Status:** Complete (5 experiments).

**q=10 raw MI-CV shows NO crossings at χ=20.** n=12 CV is below n=8 CV at ALL tested g values (gap -0.08 to -0.18). Sprint 043's χ=10 "crossing confirmed" is INVALIDATED — bad ground state at χ=10, not a real crossing.

**Dead-pair bias mechanism discovered.** At d=10, n=8 has 25% dead MI pairs (< 0.01) vs n=12 has 17%. This asymmetric fraction inflates n=8 CV relative to n=12, creating spurious crossing absence. When dead pairs are filtered: n=12 CV > n=8 CV at ALL g values — same BKT-like pattern as q=4 (no crossings, monotonically ordered).

**Energy convergence ≠ MI convergence (χ test).** At n=8 g=0.20: χ=20 CV=0.601, χ=30 CV=0.588, χ=40 CV=0.370. Energy converges (5th decimal) but mean MI jumps 44% (1.26→1.82). Dead pairs are partially χ artifacts. MI-CV at d=10 needs χ >> 40 (probably χ > d²=100).

**Raw MI-CV unreliable for d ≥ 10.** Dead-pair fraction varies systematically with n, corrupting size-dependent analyses. Filtering removes the artifact but doesn't restore crossings — q=10 filtered behavior matches q=4 (BKT-like).

**Surprises:**
- Raw MI-CV gives OPPOSITE size ordering from filtered (n12<n8 raw vs n12>n8 filtered)
- n=8 CV remarkably flat (3% variation) — d=10 smears the transition
- q=10 matches q=4 BKT pattern, NOT q=7 mean-field — BKT may recur at large q
- Sprint 043's landmark result (q=10 crossings) was doubly wrong (χ + dead pairs)

[Full report: sprints/sprint_048.md]

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
