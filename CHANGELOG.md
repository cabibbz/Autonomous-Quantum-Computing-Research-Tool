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

## Compressed Sprint History (Sprints 035-070)

- **035-048** — BW scaling, MI-CV FSS/classification, Potts phase transitions. Sprints 044-047 INVALIDATED (wrong g_c).
- **049-053** — g_c correction (true g_c from self-duality). g_c(q) formula. ν(q) ≈ 0.82-0.86.
- **054-060** — c(q) ≈ 0.40·ln(q-1)+0.55. CFT operator content. x₁ peaks at q=3. OPE coefficients.
- **061-064** — Large-q limit. c·x₁ ≈ 0.112. BW locality degrades with q.
- **065-067** — Hybrid ≠ clock universality (3 probes). No first-order signal. No floating phase in hybrid.
- **068-070** — 2D hybrid model. g_c(2D) measured. q=2 continuous confirmed. q=5 inconclusive.
- **071** — Ly=2 cylinder: g_c(q=2)=0.451, g_c(q=5)=0.714. DMRG impractical for q≥5 cylinder.
- **072** — q=3 cylinder g_c=0.565. Cyl/1D ratio monotonically decreasing. Ly=3 q=2: g_c=0.655.
- **073** — q=3 Ly=3: g_c=0.797 (49.7% to 2D). q=2 Ly=4: g_c=0.688. q=3 converges 1.9× slower.
- **074** — q=3 Ly=4 infeasible. q=5 Ly=3: g_c=0.974. Convergence saturates for q≥3.
- **075** — q=2 Ly=5: g_c=0.701 (86.6% to 2D). Power-law 1/Ly² convergence, not exponential.

Full details for all compressed sprints are in sprints/sprint_NNN.md.

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

### Sprint 077 — S_q Potts at q=7,10: g_c Found, NO First-Order Signal at n≤8
**Status:** Complete (3 experiments).

**S_q Potts g_c(q=7) = 0.144 (077a).** From (4,6) gap×N crossing. Targeted n=8 test: gap×N=0.649, matches (4,6) crossing to 0.5%. Three sizes converge within 2.4% band. Previous "69% gap collapse" (Sprint 076c) was artifact of wrong g_c estimate.

**Conformal tower at q=7 (077b).** S_q has 6-fold degenerate 1st excited (S₇ symmetry). Hybrid has 2-fold (Z₇ pairs). R_ε = 3.48 (S_q) vs 8.07 (hybrid) — 2.3x ratio. Completely different operator content.

**NO first-order signal at q=5,7,10 up to n≤8 (077c).** Δ·N drift: S_q q=5 (-4.8%), q=7 (+1.3%), q=10 (+0.6%). All <5%. Hybrid similarly stable. S_q Δ·N ≈ 0.6 nearly q-independent; hybrid Δ·N drops with q (0.48→0.25).

**S_q g_c(q=10) = 0.101** from (4,5) crossing. g_c ratio (hybrid/S_q) grows: 2.2 (q=5), 3.7 (q=7), 6.8 (q=10).

**Surprises:**
- S_q q=7 n=8 gap×N convergence BETTER than q=5 (1.3% vs 4.8% drift)
- Complex CFT walking regime extends to n≥8 even at q=10
- S_q Δ·N nearly q-independent (~0.6) while hybrid varies 2x across q=5-10
- Sprint 076c "69% gap collapse" completely explained by wrong g_c

[Full report: sprints/sprint_077.md]

### Sprint 078 — S_q Potts Self-Duality: g_c = 1/q Exact, DMRG Walking, Complex CFT c Confirmed
**Status:** Complete (3 experiments).

**g_c(S_q Potts) = 1/q exactly from self-duality (078a).** Kramers-Wannier duality verified numerically for q=2,5,7,10. All crossings approach 1/q from above with ~1/n² corrections. q=2 best: 0.50025 (dev 0.05%), q=5: 0.20076 (dev 0.4%).

**q=2 S_q ≠ hybrid (078a).** S_q field = X, hybrid = X+X† = 2X at q=2. Factor-of-2 gives S_q g_c = 0.500 ≠ hybrid 0.250. S_q = hybrid ONLY at q=3.

**No walking breakdown at n≤12 for q=5 (078b).** DMRG gap×N INCREASES: 1.88 (n=8) → 2.00 (n=12). System becomes more CFT-like, not less. Walking length ξ* > 12.

**S_q c(q=5) = 1.15, matches complex CFT Re(c) ≈ 1.138 to ~1% (078c).** Calabrese-Cardy fits with R² > 0.999. q=3 c_eff = 0.89-0.91, converging toward exact c = 4/5 with ~12% FSS overshoot. S_q c > hybrid c (1.15 vs 1.11) confirms different universality classes.

**Surprises:**
- q=2 S_q ≠ hybrid — field operator differs by factor 2
- DMRG gap×N INCREASES with n — opposite to walking breakdown
- Complex CFT quantitatively confirmed: Re(c) within 1% of DMRG measurement
- q=3 FSS overshoot persists at 12% even at n=24

[Full report: sprints/sprint_078.md]

### Sprint 079 — c(q=7,10): Walking Regime Breaks Down at Higher q
**Status:** Complete (3 experiments).

**Complex CFT c_eff ≈ Re(c) works at q=5 but FAILS at q≥7 (079a-c).** Exact diag and DMRG entropy profiles at g_c=1/q for q=7,10. Correct Coulomb gas formula: √Q=2cos(π/p), c=1-6/[p(p-1)], Re(c)=1+6α²/(π²+α²) for q>4. Verified q=2→0.5, q=3→0.8, q=5→1.138.

**c_eff/Re(c) ratio degrades monotonically with q.** At n=8: q=3 (+12%), q=5 (+0.3%), q=7 (-18%), q=10 (-40%). q=5 is the unique "sweet spot" where walking regime extends beyond all accessible sizes.

**DMRG = exact at q=7 n=8 (079b).** chi=56 DMRG and exact diag give identical c=1.1108. Zero truncation error — the low c is genuine physics, not numerical artifact.

**q=7 c_eff DECREASING with n (079a).** n=8: 1.111, n=12: 1.059. Opposite to q=3,5 (where c_eff overshoots and converges from above). The walking breakdown shows up as c_eff trending away from Re(c) at larger n.

**c_eff at n=8 ≈ 1.1 for ALL q≥5.** Nearly q-independent at moderate sizes (1.14, 1.11, 0.95 for q=5,7,10). The complex CFT q-dependence only manifests at n << ξ*(q).

**Surprises:**
- DMRG has zero chi truncation error at q=7 n=8 (ground state has low entanglement)
- c_eff(q=7) DECREASING with n — opposite to q=3,5 behavior
- q=5 unique sweet spot: only q where c_eff ≈ Re(c) to 1%
- c_eff nearly q-independent (~1.1) at moderate sizes for all q≥5

**POTENTIALLY NOVEL:** Systematic c_eff/Re(c) mapping across q=3-10 at matched sizes.

[Full report: sprints/sprint_079.md]

### Sprint 080 — c_eff at q=6,8,9: Walking Regime Boundary Mapped
**Status:** Complete (3 experiments).

**q=6 c_eff/Re(c) = 0.916 — marginal walking (080a).** c_eff = 1.147 at n=8, stable with n (+0.12% drift n=6→8). Walking extends beyond accessible exact diag sizes at q=6. Re(c) = 1.253.

**q=8 c_eff/Re(c) = 0.738, q=9 = 0.668 (080b).** Walking clearly broken. q=8 c_eff DECREASING with n (−1.0% n=6→7). All (q-1)-fold degeneracies confirmed (3-fold gaps for all q).

**Walking boundary mapped completely (080c).** c_eff/Re(c) = 1.004 × exp(−0.105(q−5)). Ratio=1 crossing at q = 5.0 ± 0.1. q=5 is the exact walking threshold. Size dependence is the key discriminator: c_eff stable for q≤6, degrading for q≥7.

**gap×N INCREASES with q** (2.01→2.17 for q=6→8) even as walking breaks down — breakdown is in entropy, not the gap.

**Surprises:**
- q=6 c_eff stable with n (+0.12%) — walking may extend to q=6 at larger n
- gap×N increases with q despite walking breakdown
- Both linear and exponential fits agree: ratio=1 at q≈5.0
- c_eff converges to ~1.0-1.15 for ALL q≥5 at moderate sizes

**POTENTIALLY NOVEL:** Complete walking boundary curve q=3-10 with exponential decay law. First measurements at q=6,8,9.

[Full report: sprints/sprint_080.md]

### Sprint 081 — ξ*(q=6) via DMRG: Walking is MARGINAL BREAKING
**Status:** Complete (3 experiments).

**q=6 c_eff drops 2.9% from n=8→12 (081a).** DMRG at g_c=1/6: n=10 (chi=40, 48s) c_eff=1.130, n=12 (chi=30, 205s) c_eff=1.115. Drift rate dc/d(ln n) = -0.048. The n=6-8 stability (+0.12%) from Sprint 080 was misleading — too short a range.

**gap×N healthy at 2.10 (081b).** DMRG excited state at n=10, chi=40 (313s). Gap still increasing with n, consistent with walking breakdown being in entropy not spectrum.

**Three-way comparison reveals crossover (081c).** dc/d(ln n): q=5 (+0.014), q=6 (-0.048), q=7 (-0.091). q=6 drift is 2× closer to q=7 than q=5. Walking boundary is a smooth crossover, not a phase boundary. Extrapolated q=6 breakdown at n≈38.

**Surprises:**
- q=6 c_eff stability at n=6-8 was MISLEADING — drop only visible at n≥10
- q=6 drift rate 3.4× faster than q=5 (opposite sign!)
- Extrapolated breakdown n≈38 — tantalizingly close to DMRG reach
- d=6 DMRG dramatically slower than d=5 (172s vs ~10s at n=8)

**POTENTIALLY NOVEL:** First systematic drift-rate comparison dc/d(ln n) across q=5-7 resolving walking crossover. q=6 breakdown length ξ*≈38.

[Full report: sprints/sprint_081.md]

### Sprint 082 — Spin-Spin Correlator: x_σ(q) Nearly Universal, Walking is Entropy-Only
**Status:** Complete (3 experiments).

**Open-BC DMRG correlator gives WRONG exponents (082a).** Raw power-law fit on open chain (q=5 n=8-24) gives η≈1.25 — inflated 5× by boundary conformal factors. R² degrades with n (0.996→0.760). Lesson: NEVER use raw open-BC power law for x extraction.

**Periodic chain correlator perfectly conformal (082b).** Exact diag on periodic chain with chord distance fit: G(r) = C × [(N/π)sin(πr/N)]^{-2x_σ}. R² ≥ 0.99999 for all q=2-8. q=2: x_σ=0.123 (exact 0.125, 1.5%), q=3: 0.132 (exact 0.133, 1.3%).

**x_σ ≈ 0.13 nearly universal across q=2-8 (082c).** Peaks at q=5 (0.136), varies ±5%. Walking vs non-walking is INVISIBLE in x_σ at accessible sizes. No oscillatory corrections detected (Im(x_σ) too small → oscillation period ~80× longer than accessible range).

**Sound velocity v(q) extracted: v = gap×N/(2π·x_σ).** Decreases monotonically: 1.02 (q=2) → 0.75 (q=5) → 0.66 (q=8). The decreasing gap×N with q is primarily velocity reduction, NOT x_σ change.

**Walking is an ENTROPY phenomenon.** Correlators (x_σ) and energy gaps (gap×N) both appear perfectly conformal for ALL q=2-8. Only c_eff (entanglement entropy) deviates for q>5. The gap-entropy decoupling extends to correlators.

**Surprises:**
- x_σ nearly constant across q — no walking signature in spin dimension
- Open-BC raw power law inflates η by 5× — critical methodological finding
- Conformal form exact to 0.04% even for q=7-8 where c_eff is breaking
- v(q) monotonically decreasing — explains gap×N trend completely
- No oscillatory corrections detectable (Im(x_σ) ≈ 0.02, period ~80 in ln(r))

**POTENTIALLY NOVEL:** First x_σ(q) measurement for S_q Potts chain via correlators. First v(q) extraction. Discovery that walking breakdown manifests in entropy/velocity but NOT in x_σ or correlator form.

[Full report: sprints/sprint_082.md]

### Sprint 083 — Casimir Energy: c_implied Matches Re(c) for ALL q, Walking is Entropy-Only
**Status:** Complete (3 experiments).

**Casimir energy formula E₀/N = ε_∞ - πvc/(6N²) fits with R² > 0.9999 for all q=2-8 (083a,b).** Periodic chain at g_c=1/q, multiple sizes per q. Extracted vc (product of velocity and central charge) independently of entropy or correlators.

**c_implied = vc/v_gap matches Re(c) to ±3% for ALL q (083c).** Defines c_implied as the central charge that makes Casimir velocity agree with gap-derived velocity. Result: c_implied/Re(c) = 0.993 (q=2), 1.005 (q=3), 1.008 (q=4), 1.027 (q=5), 1.024 (q=6), 1.017 (q=7), 0.999 (q=8). The ground state energy sees Re(c) from complex CFT even at q=8, where entropy c_eff is 40% off.

**Walking breakdown is EXCLUSIVELY an entropy phenomenon.** Energy (Casimir), spectrum (gap×N), correlators (x_σ) — all governed by complex CFT. Only entanglement entropy deviates for q>5. The reduced density matrix is uniquely sensitive to walking breakdown.

**vc(q) monotonically increasing toward ~1.** vc: 0.50 (q=2) → 0.94 (q=8). Best fit: vc = 0.31·ln(q) + 0.34. Pairwise vc decreases with N for all q (-1% to -3.6% drift), converging from above.

**Surprises:**
- c_implied/Re(c) = 0.999 at q=8 — most precise match at the q where entropy is most wrong
- Casimir energy formula works perfectly even in walking-broken regime
- Walking breakdown has zero signature in ground state energy
- vc approaching 1 at large q — may saturate

**POTENTIALLY NOVEL:** First demonstration that Casimir energy obeys complex CFT Re(c) across the walking boundary. Establishes energy-entropy hierarchy: energy observables track Re(c), entropy does not.

[Full report: sprints/sprint_083.md]

### Sprint 084 — Entanglement Spectrum: Walking Breakdown = Entropy Concentration in (q-1) Multiplet
**Status:** Complete (3 experiments).

**Entanglement spectrum has (q-1)-fold degenerate first excited level for ALL q=2-8 (084a,b).** Exact diag on periodic chains at g_c=1/q. Half-chain bipartition. S_q permutation symmetry imprints directly on entanglement spectrum: q=5 has 4-fold, q=7 has 6-fold, q=8 has 7-fold degeneracy.

**Entanglement gap Δξ ≈ 0.47 + 0.78·ln(q), INCREASES with q (084c).** Walking breakdown is NOT about closing the entanglement gap. Δξ grows from 1.03 (q=2) to 2.13 (q=8). The spectrum is more gapped, not less, at broken walking.

**Level 1 absorbs progressively more entropy: 48% (q=2) → 69% (q=8) (084c).** The (q-1)-fold degenerate multiplet captures increasing entropy weight as q grows, at the expense of both ground state (33%→19%) and tail (19%→12%). This redistribution is the microscopic mechanism for c_eff deviation.

**Tail entropy correlates with c_eff/Re(c) (Pearson r = 0.80).** Higher entanglement levels carry the walking signature. Energy/gap/correlator observables only see ground + first level (perfectly conformal). Entropy sums all levels → sensitive to tail depletion.

**Normalized R₂₁ = (ξ₂ - ξ₀)/Δξ decreases: 4.31 (q=2) → 2.42 (q=8).** Spectrum compresses above first excited level. Clean walking discriminator.

**Surprises:**
- Entanglement gap INCREASES with q — not a gap-closing phenomenon
- (q-1) degeneracy in entanglement spectrum matches energy spectrum exactly
- Level 1 entropy fraction is a monotonic function of q (48%→69%)
- Participation ratio increases with q (1.74→3.27) but this is just the (q-1) multiplicity
- Walking breakdown fully explained by a single spectral redistribution mechanism

**POTENTIALLY NOVEL:** First entanglement spectrum decomposition across walking boundary for S_q Potts chain. First identification of entropy concentration in (q-1)-fold multiplet as the microscopic mechanism for walking breakdown. No prior literature found on entanglement spectrum at q>4.

[Full report: sprints/sprint_084.md]

### Sprint 085 — Rényi Entropies: α=3 Recovers Re(c) in Walking-Broken Regime
**Status:** Complete (3 experiments).

**Rényi c_α from size pairs across q=2-8 (085a-c).** Computed S_α for α=0.5,1,2,3,5,10,∞ on periodic S_q Potts chains. Extracted c_α using CFT formula for periodic chain (two entanglement cuts): c_α = 6·ΔS_α/((1+1/α)·ln(N₂/N₁)). Initially used wrong prefactor (c/12, one cut) — corrected to c/6 (two cuts).

**α=1 (von Neumann) is MOST accurate for walking q≤6 (085c).** c₁/Re(c): q=2: 0.996, q=5: 0.999 (0.06% accuracy!), q=6: 0.995. von Neumann entropy is uniquely accurate in the walking regime — counterintuitive given it sums all eigenvalues.

**α=3 uniquely recovers Re(c) in walking-broken regime (085c).** c₃/Re(c): q=7: 1.027 (2.7%), q=8: 0.991 (0.9%). The "magic" Rényi index shifts from α≈1 (walking) to α≈3 (broken) — a crossover in optimal probe.

**c_∞ does NOT recover Re(c) — original hypothesis WRONG (085c).** c_∞ deviates comparably to c₁ for q≥7 (both ~7% off). Spectral redistribution affects ALL eigenvalues including λ_max.

**Rényi spread (c₂-c_∞)/Re(c) is a new monotonic walking discriminator (085c).** Values: -0.046 (q=2), +0.008 (q=3), +0.086 (q=5), +0.106 (q=6), +0.118 (q=7), +0.129 (q=8). Changes sign near q=3, monotonically increases. Cleaner than c_eff/Re(c).

**α=2 always overshoots Re(c) (085c).** c₂/Re(c) = 1.02-1.13 for ALL q. Most stable (smallest variation across q) but most biased.

**Surprises:**
- von Neumann (α=1) is literally the most accurate Rényi entropy for walking regime
- α=3 (not α=∞) is the "corrective" index in broken regime — deep non-triviality
- Rényi spread is monotonic in q — better discriminator than c_eff/Re(c)
- Factor-of-2 periodic/open BC distinction nearly led to wrong physical conclusions
- c_α profile shape changes at walking boundary: peak shifts from α=3-5 to α=2

**POTENTIALLY NOVEL:** First systematic Rényi c_α(q,α) mapping across walking boundary for S_q Potts chain. ~~Discovery of α=3 as optimal probe for walking-broken regime~~ (RETRACTED Sprint 086 — single-size extraction artifact). Rényi spread as new walking discriminator (needs recheck with size pairs).

[Full report: sprints/sprint_085.md]

### Sprint 086 — Rényi Entropy Scaling via DMRG: α=3 Was an Extraction Artifact
**Status:** Complete (3 experiments).

**DMRG Rényi entropies at q=5 n=8,10,12 and q=7 n=6,8,10,12 (086a-b).** Open-BC DMRG Schmidt spectra at g_c=1/q. Profile fits to CC gave poor R² (0.79-0.88) due to boundary corrections. Midchain S_α extracted for size-pair analysis.

**Sprint 085's α=3 finding was a SINGLE-SIZE EXTRACTION ARTIFACT (086c).** Single-size formula c_α = 12·S/((1+1/α)·ln(N/π)) includes non-universal constant c'_α which varies with α. At n=6-7 periodic BC, c'₃ coincidentally compensates walking deviation. Size-pair extraction (cancels c'_α) shows α=1-2 is best for ALL q.

**Corrected optimal α progression: 0.5 (real CFT) → 1 (walking) → 2 (broken).** Periodic size-pair: q=2,3 best at α=0.5; q=5,6 at α=1; q=7,8 at α=2. Smoother progression than Sprint 085 suggested.

**No Rényi index recovers Re(c) for broken walking.** q=7 periodic pair: best c_α/Re(c) = 0.90 (α=2). q=8: 0.88. Walking breakdown is genuine for ALL Rényi indices.

**DMRG walking breakdown in c₁ size pairs (086b).** q=7 c₁/Re(c): 1.05 (6,8 pair) → 0.83 (10,12 pair). 22% drop. q=5: only 5% drop. Clear walking vs stable walking separation.

**Entanglement tail grows 8× at q=7 n=6→12.** λ_max: 0.891→0.852. (q-1) multiplet: 0.109→0.147. Tail: 0.014%→0.11%. Direct microscopic evidence of progressive walking breakdown.

**Surprises:**
- Single-size vs size-pair extraction gives qualitatively different "optimal α"
- Open-BC profile R² < 0.9 even at n=12 — boundary corrections dominate
- Periodic and open BC give different optimal α (periodic: α=2, open: α=1 for q=7)
- q=7 entanglement tail grows 8× while λ_max drops only 4%

[Full report: sprints/sprint_086.md]

### Sprint 087 — Entanglement Spectrum Scaling via DMRG: Tail Weight is UNBOUNDED
**Status:** Complete (3 experiments).

**q=5 entanglement spectrum n=8-24 (087a).** DMRG midchain Schmidt spectrum at g_c=1/5. λ_max decreases 0.858→0.787, (q-1) multiplet grows 14.2%→20.9%, tail weight grows 0.05%→0.44%. All trends monotonic, no saturation.

**q=7 entanglement spectrum n=6-12 (087b).** g_c=1/7, limited to n=12 (293s). λ_max: 0.891→0.852, multiplet: 10.9%→14.7%, tail: 0.014%→0.11%.

**Tail weight grows as UNBOUNDED power law (087c).** q=5: w_tail ~ n^1.72 (R²=0.993). q=7: w_tail ~ n^2.52 (R²=0.989). q=7 exponent larger but prefactor 9× smaller. No saturation. Entanglement gap closes logarithmically: Δξ = -0.43·ln(n) (q=5), -0.50·ln(n) (q=7).

**Energy-entropy decoupling is TEMPORARY.** Tail weight reaches 10% at n≈147 (q=5), n≈72 (q=7). Walking breakdown creates a hierarchy of crossover scales. %S(lev0) saturates (~0.225), weight flows from multiplet to tail.

**c_eff size pairs from DMRG.** q=5: converges 1.165→1.010 ×Re(c) from (8,12)→(20,24). q=7: diverges 1.055→0.841. Walking vs breakdown cleanly separated.

**Surprises:**
- Tail weight is UNBOUNDED — no saturation even at n=24
- q=7 power law exponent 47% larger (2.52 vs 1.72) — breakdown accelerates
- Entanglement gap closes logarithmically — will vanish at n ~ 10⁴
- Weight flows lev1 → tail, NOT lev0 → tail
- %S(lev0) approaches fixed point (~0.225 for q=5)

**POTENTIALLY NOVEL:** First power-law scaling of entanglement spectrum tail weight across walking boundary. First demonstration that energy-entropy decoupling is finite-size with quantitative crossover scales.

**NOTE (Sprint 088 correction):** Exponents 1.72 (q=5) and 2.52 (q=7) were from real-space curve_fit (biased by large-n). Log-log fit gives b≈2.0 for q=2,3,5 — universal, not walking-specific.

[Full report: sprints/sprint_087.md]

### Sprint 088 — Tail Weight at q=2,3: UNIVERSAL Power Law, Not Walking-Specific
**Status:** Complete (3 experiments).

**q=3 entanglement spectrum n=8-24 (088a).** DMRG at g_c=1/3 (c=4/5 real CFT). λ_max: 0.860→0.791, w_mult: 0.139→0.205, w_tail: 0.00047→0.00429. (q-1)=2 multiplet perfectly degenerate. %S(tail) reaches 4.2% at n=24.

**q=2 entanglement spectrum n=8-24 (088b).** DMRG at g_c=1/2 (c=1/2 Ising). λ_max: 0.890→0.839, w_mult: 0.110→0.158, w_tail: 0.00037→0.00324. %S(tail) reaches 4.4% at n=24.

**Tail weight exponent is UNIVERSAL: b ≈ 2.0 for q=2,3,5 (088c).** Log-log power-law fit: q=2 b=1.98, q=3 b=2.01, q=5 b=2.01 — identical within fit error. q=7 shows b=2.97 but only 4 data points (n=6-12), uncertain. Sprint 087's exponents (1.72, 2.52) were real-space fitting artifacts.

**Walking breakdown is NOT in tail growth.** %S(tail) ≈ 4% at n=24 for ALL q=2,3,5. Entanglement gap closure Δξ ~ -0.43·ln(n) also universal. Walking-specific effects must reside in level REDISTRIBUTION (how multiplet vs ground state share weight), not in the tail.

**Surprises:**
- q=2,3,5 tail exponents identical to <1% — far more universal than expected
- Real-space vs log-log fitting gives qualitatively different conclusions (1.72 vs 2.01)
- %S(tail) at n=24 is nearly q-independent (3.8-4.4%)
- Sprint 087 narrative partially corrected: tail growth is universal, only redistribution is q-dependent

[Full report: sprints/sprint_088.md]
