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

- **076** — S_q Potts vs hybrid head-to-head at q=5. g_c(S_q)=0.200. 4-fold vs 2+2 degeneracy. S_q FSS 10x larger.
- **077** — S_q at q=7,10: g_c=0.144/0.101. No first-order signal at n≤8. Δ·N≈0.6 q-independent.
- **078** — Self-duality: g_c=1/q exact. DMRG walking holds at n≤12. c(q=5)=1.15 matches Re(c)=1.138.
- **079** — c_eff at q=7,10: walking breaks down. c_eff/Re(c) degrades monotonically. q=5 unique sweet spot.

Full details for sprints 076-079 are in sprints/sprint_NNN.md.

- **080** — c_eff at q=6,8,9: walking boundary mapped. c_eff/Re(c) = 1.004·exp(−0.105(q−5)). q=5 exact walking threshold. gap×N increases with q even as walking breaks.
- **081** — q=6 marginal breaking: c_eff drops 2.9% n=8→12, drift dc/d(ln n)=-0.048. Walking boundary is smooth crossover. ξ*≈38 for q=6.

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

### Sprint 089 — q=7 Tail Exponent Resolution & Level Redistribution as Walking Discriminator
**Status:** Complete (3 experiments).

**q=7 tail exponent b≈3.0 was PRE-ASYMPTOTIC (089a).** Extended DMRG to n=14 (chi=30, 237s) and n=16 (chi=25, 138s). With n≥8 data (5 points): b = 2.07 ± 0.14, consistent with universal b≈2.0. Sprint 088's b≈3.0 from 4 points (n=6-12) resolved — n=6 was in a crossover regime.

**Static entropy partition is the walking discriminator (089b).** %S(lev0) and %S(lev1) are MONOTONIC in q. Weight always flows multiplet→tail. d(%S)/d(ln n) rates are NOT monotonic — dynamic rates don't cleanly discriminate. %S(lev0) saturation: 0.338 (q=2) → 0.203 (q=7), R² > 0.999.

**Multiplet dominance M/[(q-1)/q] crosses 1.0 at walking boundary (089c).** M = %S(lev1)/[%S(lev0)+%S(lev1)]. At n=16: M = 0.676 (q=2), 0.724 (q=3), 0.771 (q=5), 0.797 (q=7). Ratio M/[(q-1)/q]: 1.33 (q=2), 1.07 (q=3), 0.96 (q=5), 0.93 (q=7). Crosses 1.0 between q=3 and q=5 — the real-to-complex CFT boundary. Democracy index flips: ground-depleted for real CFT, multiplet-depleted for complex CFT.

**Mean tail exponent across all q: b = 2.017 ± 0.032.** Universal.

**Surprises:**
- q=7 b≈3.0 completely resolved — was pre-asymptotic crossover at small n
- M/[(q-1)/q] = 1.0 crossing coincides with real→complex CFT boundary
- Democracy index flips sign at q≈4 — unexpected structural connection
- c_eff/Re(c) anticorrelates with M (r = −0.82) — multiplet trapping explains c_eff deviation
- %S(lev1) → 0 as n→∞ for all q — multiplet is transient entropy reservoir

**POTENTIALLY NOVEL:** First identification of multiplet dominance M as walking discriminator. M/[(q-1)/q] crossover at q=4 connects real-to-complex CFT transition to entanglement spectrum entropy partition. Universal tail exponent b=2.0 confirmed for all q=2-7.

[Full report: sprints/sprint_089.md]

### Sprint 090 — q=4 Entanglement Spectrum: M/[(q-1)/q] Crossover CONFIRMED at q≈4
**Status:** Complete (3 experiments).

**q=4 DMRG entanglement spectrum n=8-24 (090a).** g_c=1/4, chi=50-130. M/[(q-1)/q] = 1.021 (n=8), 1.010 (n=12), 1.003 (n=16), 0.998 (n=20), 0.994 (n=24). **Crosses 1.0 between n=16 and n=20.** (q-1)=3-fold multiplet perfectly degenerate. %S(lev0) saturates at 0.251 (≈1/4). c_eff converges from above: 1.25→1.09.

**q=4 tail exponent b = 2.024 — universal b≈2.0 confirmed (090b).** Updated mean across q=2,3,4,5,7: b = 2.018 ± 0.029. q=4 slots perfectly into universal pattern. %S(tail) at n=24 = 4.0%, between q=3 (4.2%) and q=5 (3.8%).

**M/[(q-1)/q] = 1.0 crossover at q≈4 CONFIRMED (090c).** Full curve at n=16: 1.352 (q=2), 1.086 (q=3), **1.003 (q=4)**, 0.964 (q=5), 0.930 (q=7). Linear interpolation: q_cross = 4.07. **q_cross converges to q≈4.0 as n→∞:** 4.49 (n=8) → 3.92 (n=24). The democracy index quantitatively pins the real-to-complex CFT boundary.

**Surprises:**
- M/[(q-1)/q] at q=4 n=16 is 1.003 — within 0.3% of exactly 1.0
- q_cross converges monotonically toward q=4 with increasing n
- Periodic vs open BC give different M/[(q-1)/q] at same finite n (0.95 vs 1.02 at n=8)
- %S(lev0) → 0.251 ≈ 1/4 for q=4 — may saturate at exactly 1/q

**POTENTIALLY NOVEL:** First entanglement spectrum measurement at q=4 S_q Potts. First confirmation that democracy index crossover occurs at q=4, pinpointing the real-to-complex CFT boundary in the entanglement spectrum. q_cross → 4.0 in the thermodynamic limit.

[Full report: sprints/sprint_090.md]

### Sprint 091 — Entanglement Hamiltonian H_E: BW Fidelity Across Walking Boundary
**Status:** Complete (3 experiments).

**BW fidelity confounded by nA in naive comparison (091a).** Periodic chain at g_c=1/q, q=2-7 with different n (and thus nA). q=7 (nA=3) showed 99.9% BW fidelity, q=2 (nA=6) showed 82.3% — but this was an nA artifact, not a q effect. Lesson: always compare at fixed nA.

**At fixed nA=4, non-Potts fraction grows EXPONENTIALLY with q (091b).** q=2: 0.029%, q=3: 0.081%, q=4: 0.955%, q=5: 2.729%. Fit: ~exp(1.6·q). NNN and 3-body Potts operators contribute <0.2% — the residual is genuinely non-Potts operators absent from the physical Hamiltonian.

**Biggest BW jump at q=3→4 (11.7×), NOT q=4→5 (2.9×) (091c).** The real-to-complex CFT boundary at q=4 is ALSO a BW locality boundary. d(log non-Potts)/dq spikes to 2.46 at q=3→4, vs 1.05 elsewhere. BW alpha also jumps: 2.4 (q=2,3) → 3.2 (q≥4).

**q=5 non-Potts grows 1.7× faster with nA than q=2.** BW locality scaling: q=2 drops from 99.98% (nA=3) to 62.6% (nA=7), q=5 from 99.93% to 97.1%. Walking-enhanced corrections amplified at larger subsystems.

**Surprises:**
- BW fidelity dominated by nA, not q — methodological caution for future comparisons
- Biggest non-Potts jump at q=3→4 (11.7×) coincides with real-to-complex CFT boundary
- Non-Potts operators dominate BW residual — NNN/3-body Potts operators negligible
- BW alpha discontinuity at q=4: jumps from 2.4 to 3.2
- q=5 non-Potts slope 1.7× steeper than q=2 — walking amplifies BW corrections

**POTENTIALLY NOVEL:** First BW fidelity across walking boundary for S_q Potts. First identification of q=4 as BW locality boundary. First demonstration that BW corrections are dominated by non-Potts operators.

[Full report: sprints/sprint_091.md]
