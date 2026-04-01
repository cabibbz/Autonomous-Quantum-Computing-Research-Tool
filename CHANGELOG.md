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

- **035** — BW size scaling: Pauli fraction fails (8% at n=20), spectrum R²=1.0. BW form correct but envelope wrong.
- **036** — MI-CV FSS: genuine divergence ~n^1.1 at TFIM. Correlation-function MI reconstruction validated.
- **037** — MI-CV classification complete: crossing (Ising), step (1st-order), dome (BKT).
- **038** — Data collapse confirms Ising ν=1. q=3 Potts MI-CV shows crossings (at wrong g, see Sprint 049).
- **039** — Potts data collapse: ν=5/6 at wrong g. Slope exponent as second discriminator.
- **040** — q=4 Potts MI-CV: crossings at g≈0.893. Slope lower than q=3 (log corrections at marginal q=4).
- **041** — q=5 clock MI-CV: crossings persist at g≈0.673. Clock ≠ Potts for q≥4.
- **042** — True q=5 Potts MI-CV: second-order at g≈0.41, NOT first-order. Custom PottsChain model built.
- **043** — q=10,20 Potts: all second-order. Direct MPS contraction for ρ_ij. Universal large-q regime.
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
