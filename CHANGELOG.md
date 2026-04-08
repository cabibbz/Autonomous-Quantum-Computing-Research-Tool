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

- **082** — x_σ(q) ≈ 0.13 nearly universal for q=2-8. v(q) monotonically decreasing. Walking invisible in correlators. Open-BC power law inflates η by 5×.
- **083** — Casimir c_implied/Re(c) ≈ 1.00±0.03 for ALL q=2-8. Energy tracks Re(c) even where entropy deviates 40%. Walking is exclusively an entropy phenomenon.
- **084** — Entanglement spectrum: (q-1)-fold degenerate multiplet. Level 1 absorbs 48%→69% of entropy (q=2→8). Entropy concentration as microscopic walking breakdown mechanism.

- **085** — Rényi c_α(q,α) mapping. α=1 best for walking (q≤6), optimal α shifts to 2 for broken (q≥7). ~~α=3 recovery~~ RETRACTED in Sprint 086. Rényi spread is monotonic walking discriminator.

- **086** — Rényi DMRG: α=3 recovery was extraction artifact. Corrected: optimal α = 0.5 (real) → 1 (walking) → 2 (broken). q=7 tail grows 8×. No Rényi index recovers Re(c) for broken walking.
- **087** — Entanglement spectrum DMRG: tail weight unbounded power law ~n^b. c_eff size pairs confirm walking vs breakdown. NOTE: exponents corrected in Sprint 088.
- **088** — Tail exponent UNIVERSAL b≈2.0 for q=2,3,5. Walking breakdown NOT in tail growth — must be in level redistribution.
- **089** — q=7 b≈3.0 resolved (pre-asymptotic). M/[(q-1)/q] crosses 1.0 at walking boundary. Democracy index flips at q≈4.
- **090** — q=4 M/[(q-1)/q]=1.003 at n=16. q_cross→4.0 in thermodynamic limit. b=2.024 confirms universal.
- **091** — BW fidelity across walking: non-Potts grows ~exp(1.6q). Biggest jump at q=3→4. BW alpha discontinuity at q=4.
- **092** — Non-Potts operators = mixed clock-shift (XZ-type). Walking amplifies but doesn't change type. Flat at nA=3.
- **093** — BW temperature profile: sin-envelope correct, residual bulk-concentrated (real CFT) vs uniform (walking). 34× amplification at q=5.
- **094** — BW threshold nA*(q): two-regime, B(q)=0.48q+1.09. Walking shifts nA* down by 1 site. NOTE: nA*=5 was artifact (Sprint 095).

Full details for sprints 086-094 are in sprints/sprint_NNN.md.

### Sprint 095 — BW R² vs n at Fixed nA: Equal-Bipartition Anomaly Discovered
**Status:** Complete (3 experiments).

**q=2 BW R² with varying n at fixed nA (095a).** Periodic exact diag, nA=3-6, n up to 18. Key: 1-R² is flat in n for nA≤4 (varies <2× as n doubles). nA=6 collapses at ALL n (1-R²≈0.17 at n=18, nA/n=0.33). nA=5 at n≥16: 1-R²=5.6e-4 (20× better than Sprint 094's n=10).

**q=5 periodic nA=3 (095b).** n=7: 1-R²=1.66e-3, n=8: 2.11e-3. Walking amplification 10-15× at matched nA/n. DMRG attempt at chi≥125 too expensive (251s at n=10).

**Equal-bipartition penalty discovered (095c).** At nA/n=0.5, BW accuracy is systematically worse: 2× (nA=3), 2.1× (nA=4), **18× (nA=5)**, 1.4× (nA=6). Sprint 094's "threshold at nA=5" was an ARTIFACT. True threshold: nA*=6 for q=2.

**Surprises:**
- nA=5 is NOT at the BW threshold — 18× equal-bipartition artifact inflated Sprint 094's result
- BW corrections are UV (lattice), not IR — increasing n doesn't help
- nA=6 threshold is genuine: persists at nA/n=0.33
- Equal-bipartition penalty peaks at threshold nA — most dangerous where it matters most
- Walking amplification (q=5/q=2) is ~10-15× regardless of nA/n ratio

**POTENTIALLY NOVEL:** First BW fidelity vs subsystem-to-system ratio mapping. Discovery of equal-bipartition anomaly in BW accuracy. Revision of BW threshold from nA*=5 to nA*=6 for q=2.

[Full report: sprints/sprint_095.md]

### Sprint 096 — BW Threshold Mechanism: Non-Potts Operators, Not Spectrum
**Status:** Complete (4 experiments).

**Entanglement spectrum smooth across BW threshold (096a-b).** q=2 n=14 nA=3-7: tail weight changes 8% from nA=5→6, but 1-R² jumps 212×. Confirmed for q=3,5. Spectrum does NOT predict BW breakdown.

**BW residual dominated by max-range operators (096c).** At nA=6, 96% of residual in operators spanning ≥4 sites. Top operators: ZIIIXI, ZIIXXI — non-Potts mixed types.

**Adding longer-range Potts operators barely helps (096d).** Free-fit 21 Potts ops: 1-R² improves from 0.171 to 0.163 (5%). Non-Potts fraction itself has 244× threshold jump at nA=6. BW envelope is optimal for Potts subspace.

**Mechanism:** BW breaks when H_E develops significant non-Potts operators (mixed clock×shift). Threshold is in EIGENVECTORS of ρ_A, not eigenvalues. UV lattice effect.

[Full report: sprints/sprint_096.md]

### Sprint 097 — H_E Operator Compactness: BW Breakdown is Fundamental
**Status:** Complete (3 experiments).

**Pre-threshold H_E is maximally compact (097a).** q=2 n=14, full Pauli decomposition. nA=3-5: exactly 2nA-1 BW operators capture >99.9%. PR=4-7. At nA=6 (threshold): need 206 ops for 90%, 1052 for 99%. PR jumps to 12.

**Non-BW content is DIFFUSE at threshold (097b).** nA=6: non-BW PR=357, need 743 operators for 90% of non-BW weight. XZ-mixed 44%, XYZ 21%. High-body (4-6) long-range (3-5) operators dominate. No compact "BW + corrections" ansatz exists.

**Universal across q (097c).** q=3 n=10 clock-shift basis: nA=3,4 compact (10-14 BW ops, >99.9%). nA=5 threshold: non-Potts 22.4% (even larger than q=2 at 16.3%). Free-fit all Potts ops only 3% improvement.

**BW is the optimal compact H_E approximation.** Beyond threshold, H_E requires exponentially many parameters — fundamentally beyond any local Hamiltonian ansatz.

[Full report: sprints/sprint_097.md]

### Sprint 098 — Harden Casimir Finding: CONFIRMED NOVEL
**Status:** Complete (2 experiments).

**GPU-extended Casimir fit (098a).** Vectorized Hamiltonian builder for large matrices. New sizes: q=3 n=12 (531k), q=5 n=10 (9.8M, 159s), q=7 n=8 (5.7M, 148s), q=8 n=7 (2.1M). All q now have 4-5 data points.

**Pairwise convergence analysis (098b).** Pairwise c_implied/Re(c) from consecutive (N₁,N₂) pairs: q=7 (7,8)→1.000 (exact!), q=8 (6,7)→0.992. Extrapolated ∞: mean=0.998±0.012. 1/N⁴ corrections systematic (|d|/q≈0.044, q-independent).

**Casimir 16× more consistent with Re(c) than entropy.** Pairwise-last spread 0.009 vs entropy spread 0.145. At q=7: Casimir off by 0.0%, entropy off by 21.6%.

**Upgraded from POTENTIALLY NOVEL to CONFIRMED NOVEL.** Five checks passed: 5+ points, pairwise convergence, 16× consistency, systematic corrections, independence from entropy.

[Full report: sprints/sprint_098.md]

### Sprint 099 — Complex CFT: Im(c) Oscillation Detection (NEGATIVE)
**Status:** Complete (3 experiments).

**Dense Casimir scan (099a).** E₀(N) at all integer N: q=2 N=4-14 (11 pts), q=5 N=4-10 (7 pts, GPU for N=9,10), q=7 N=4-8 (5 pts, GPU for N=8). Odd sizes (N=5,7,9) newly measured for q=5.

**Oscillatory vs monotonic fit (099b).** Compared 3 models at same parameter count. Monotonic 1/N⁶ correction beats oscillatory cos(ω·ln N)/N⁴ for ALL q (M3/M2 = 5051 for q=2, 1.54 for q=5, 1.61 for q=7). Complex CFT oscillation period (6.53 in ln N for q=5) requires N ratio ~684 — our N=4-10 range covers only 14% of one cycle.

**Pairwise separation (099c).** q=5 pairwise c/Re(c) is non-monotonic (first decreases, then increases from N≈6). Traced to VELOCITY effect: pairwise vc itself is monotonically decreasing for all q. Detrended vc RMS: q=2 4.75e-5, q=5 3.97e-4 (8×), q=7 9.96e-4 (21×). Richardson extrapolation drifts for q>4.

**Im(c) oscillations undetectable at exact diag sizes.** Would need DMRG at N=10-100 or non-Hermitian formulation.

[Full report: sprints/sprint_099.md]

### Sprint 100 — DMRG Casimir: Open vs Periodic BC
**Status:** Complete (3 experiments).

**DMRG Casimir q=2 open BC (100a).** 7 points N=8-24. DMRG matches exact diag to 10⁻¹⁴ at N≤14. Open-BC 4-param fit: c=0.485 (c/0.5=0.971). Boundary corrections ~3% at N~24. Pairwise converges slowly: 0.88→0.93.

**Open-BC Casimir for q=5,7 — FAILS (100b).** q=5 (6 pts N=6-12): c/Re(c)=0.64. q=7 (5 pts N=4-8): c/Re(c)=0.36. Boundary corrections scale with q: 3% (q=2), 36% (q=5), 64% (q=7). Open-BC extraction is impractical for q≥5 at accessible DMRG sizes.

**Periodic-BC reanalysis — confirms Re(c) tracking (100c).** Sprint 099a dense periodic data reanalyzed. Pairwise c/Re(c): q=2 (5,6)→1.000, q=5 (7,8)→1.002, q=7 (5,6)→1.000. Periodic BC 10-100× more accurate than open BC at same N (1/N² vs 1/N corrections).

**Key finding:** Casimir-Re(c) result CANNOT be extended to DMRG sizes via open BC. Periodic exact diag (Sprint 098) remains definitive. Literature confirms Im(c) detection impossible from real Hamiltonian Casimir energy.

[Full report: sprints/sprint_100.md]

### Sprint 101 — Symmetry-Resolved Entanglement Entropy (SREE)
**Status:** Complete (3 experiments, 101a partially retracted).

**SREE for S_q Potts q=2-10 at g_c=1/q (101b).** Charge-sector decomposition via Z_q projectors. Charge-0 enrichment increases with q: p(0)*q from 1.55 (q=2) to 5.15 (q=10). p(0)→1/2 as q→∞ (not 1/q). S_n/S_t ≈ 0.908 universal across q at n=6.

**Size scaling (101c).** S_n/S_t is q-INDEPENDENT at fixed geometry: 0.953 (n=4), 0.908 (n=6), 0.875 (n=8) for all q tested. Implies S_number ∝ c(q). p(0)*q decreases slowly with n but stays >>1 at accessible sizes.

**101a RETRACTED:** Used hybrid field (X+X†) with S_q critical coupling g_c=1/q. For q≥4 this put the system far from criticality. Corrected in 101b.

**SREE is NOT a walking discriminator** at accessible sizes — all q-dependence is smooth and monotonic. No feature at q=4 or q=5.

**POTENTIALLY NOVEL:** Universal q-independent S_n/S_t ratio at fixed geometry for critical S_q Potts. First SREE measurement for q=4-10 Potts model.

[Full report: sprints/sprint_101.md]

### Sprint 102 — Fidelity Susceptibility Across the Walking Boundary
**Status:** Complete (3 experiments).

**χ_F scan near g_c for q=2,3,5,7 (102a).** S_q Potts periodic chain, δg=10⁻⁴. Peak heights at n=6: 0.73 (q=2), 4.09 (q=3), 36.29 (q=5), 166.88 (q=7) — 228× range. Peak at exact g_c for q≥5 (no FSS shift). Peak converges from below for q=2,3.

**FSS analysis + q=2 n=14 (102b).** χ_F_max ~ N^α. q=2: α=0.98 → ν=1.009 (exact 1.0). q=3: α=1.38 → ν=0.841 (exact 0.833). q=5: α=2.09 → ν_eff=0.648 (gap-derived 0.83). q=5 exponent matches first-order prediction (α=2 in 1D).

**q=4 crossover (102c).** n=6,8: α=1.69, ν=0.743 (exact 2/3). χ_F at n=6 grows exponentially with q: ~exp(1.06q), growth 2.88× per unit q.

**Key finding:** Walking regime makes χ_F scale as ~N² (first-order) despite gap/correlators showing continuous behavior. Fourth observable with walking-specific anomaly (after entropy, Casimir, correlators).

**POTENTIALLY NOVEL:** First fidelity susceptibility measurement across walking boundary. Discovery of first-order-like χ_F scaling in walking regime. Need more sizes at q=5 to harden.

[Full report: sprints/sprint_102.md]

### Sprint 103 — Harden χ_F Scaling: α(q) Confirmed Novel
**Status:** Complete (3 experiments).

**q=5 GPU extension (103a).** n=6,8,9,10 (4 sizes, GPU for n≥9). α=2.091 (4-point fit, err=0.002). Pairwise: (6,8)→2.088, (8,9)→2.094, (9,10)→2.100 — converging UPWARD. n=10 (dim=9.8M) at 156s on GPU.

**q=6,7 mapping (103b).** q=6 n=6,8: α=2.37. q=7 n=6,7,8: α=2.65 (pairwise 2.63, 2.67). Super-first-order scaling (α>2) for ALL q≥5, increasing with q.

**Full α(q) analysis (103c).** α(q) = 0.315q + 0.469 linear for q≥4 (residuals ±0.05). α=2 crossing at q≈4.9 (walking boundary). χ_F(n=6) ~ exp(1.06q). q=4 α=1.69 (15% below BKT prediction due to log corrections).

**Upgraded to CONFIRMED NOVEL.** Five checks: 4+ sizes at q=5 (0.3% stability), three walking q values, q=2,3 match exact ν to <2%, independent reproduction of q=7 n=6, simple linear functional form.

[Full report: sprints/sprint_103.md]

### Sprint 104 — Energy-Entropy Hierarchy Universality Test (NEGATIVE)
**Status:** Complete (3 experiments).

**J1-J2 chain at three regimes (104a).** Spin-1/2 chain, periodic BC, S_z=0 sector exact diag, N=8-20. Compared c_eff (entropy) vs c_Cas (Casimir) convergence to c=1. XX (Δ=0): tied (1.0×). Heisenberg (Δ=1): entropy wins 21× (c_eff 0.07% off, c_Cas 1.4% off). BKT (J2=0.2412): Casimir wins 14× (c_Cas 0.03% off, c_eff 0.37% off).

**Velocity extraction (104b).** Independent v from S_z=0→S_z=1 gap. XX: v_gap=1.002 (exact 1.0). Heisenberg: v_gap=1.387 (exact π/2=1.571, 12% log correction). BKT: v_gap=1.177 (smoothly converging, minimal log correction). Multi-size Casimir fit confirms: c_Cas(fit) = 0.998-1.005 for all models.

**Cross-model comparison (104c).** Potts max |Δc_eff| = 26% (q=7-8). J1-J2 max |Δc_eff| = 0.6%. Potts effect is 41× LARGER. Heisenberg c_Cas deviation follows 1/ln(N) (fit: 0.283/ln(N) - 0.084).

**Key finding: Energy-entropy hierarchy is NOT universal.** Direction depends on correction type: marginal operator logs → entropy wins (Heisenberg); BKT log³ corrections → Casimir wins; walking spectrum reorganization → Casimir wins by O(1) (Potts only). The O(1) entropy deviation is unique to walking.

[Full report: sprints/sprint_104.md]

### Sprint 105 — χ_F Scaling at J1-J2 BKT: Walking Super-Scaling is Unique
**Status:** Complete (3 experiments).

**Wide χ_F scan (105a).** J1-J2 Heisenberg, S_z=0 periodic BC, N=8-20. No peak near J2_c=0.2412 at any size — χ_F monotonically increasing through BKT region. Peak instead at J2≈0.487 (Majumdar-Ghosh point). N=16 has level-crossing spike.

**BKT vs MG resolution (105b).** BKT region (J2=0.15-0.35): MONOTONIC for all N. MG region (J2=0.40-0.55): peak saturates at N≈16 (χ_F=1.88) then DECREASES. Pairwise α: 4.24→0.28→-0.16. MG is a finite-size level crossing, not divergent.

**Exact-point scaling (105c).** χ_F at J2_c=0.2412: grows as N^0.58 (global), but pairwise α decreases 1.01→0.29 — trending to logarithmic. Three-way comparison:
- BKT: α→0 (invisible, exponential gap closing too slow)
- MG first-order: saturates then decreases (finite discontinuity)
- Potts walking: α=2.09-2.65, increasing with N (genuine super-first-order)

**Walking χ_F super-scaling confirmed as walking-specific.** Combined with Sprint 104 (Casimir-Re(c) walking-specific), both confirmed novel results now properly scoped.

[Full report: sprints/sprint_105.md]

### Sprint 106 — χ_F Spectral Decomposition: Walking Super-Scaling Mechanism
**Status:** Complete (3 experiments).

**Spectral decomposition at g_c (106a).** S_q Potts periodic chain, 20 eigenstates per config. χ_F is 100% dominated by a SINGLE excited state: the (q-1)-fold degenerate S_q multiplet (level q-1). The spectral gap (level 1) has ZERO matrix element with H_field — symmetry-forbidden. For q≥3, one multiplet captures 100%; for q=2, 91-93%.

**Size scaling of components (106b).** Decomposition: α = β_me + 2z_m - 1, where gap_m ~ N^{-z_m} and |me|² ~ N^{β_me}. Extended to n=14 (q=2), n=10 (q=3), n=9 (q=5, GPU), n=7 (q=7, GPU). Decomposition is exact (predicted = measured to all digits). z_m: 0.989 (q=2) → 1.310 (q=7). β_me: 0.066 (q=2) → 1.010 (q=7). Both increase linearly with q.

**Mechanism analysis (106c).** Added q=4 (5 sizes) and q=6 (3 sizes). Linear fits: z_m(q) = 0.065q + 0.843, β_me(q) = 0.182q - 0.201. Reconstructed α(q) = 0.312q + 0.485 matches known 0.315q + 0.469 within 3%.

**Two sources of walking super-scaling:** (1) Multiplet gap closes faster than 1/N (z_m > 1 for q≥3), (2) Field matrix element grows with system size (β_me > 0 for q≥3). For q=2 (real CFT): gap only. For q=7 (broken walking): both contribute (72%/28%).

**POTENTIALLY NOVEL:** First spectral decomposition of χ_F at walking transition. Spectral gap symmetry-forbidden. Exact α = β_me + 2z_m - 1 with linear q-dependence.

[Full report: sprints/sprint_106.md]

### Sprint 107 — Harden χ_F Spectral Decomposition: CONFIRMED NOVEL
**Status:** Complete (4 experiments).

**GPU extension (107a/107a2).** q=5 extended to 5 sizes (n=6-10, GPU for n=10 dim=9.8M, 279s). q=7 extended to 3 sizes (n=6-8, GPU for n=8 dim=5.8M, 185s). All sizes: dominant fraction = 1.000. Updated fits: z_m(q) = 0.065q + 0.845, β_me(q) = 0.188q − 0.238. α = β_me + 2z_m − 1 exact to machine precision.

**Cross-check (107b).** Spectral vs finite-difference χ_F compared at q=2,3,5. Spectral/fd ratio = captured fraction: q=5 ~98% (single multiplet), q=2 ~87-93% (multi-state). Confirms mechanism: walking regime is truly single-multiplet dominated.

**Entanglement connection (107c).** Energy multiplet gap (z_energy ≈ 1.0-1.3) and entanglement gap (z_ent ≈ 0.2) share S_q symmetry origin but scale independently. Δξ/gap_m grows with N. The two (q-1)-fold multiplets are independent physical phenomena.

**Upgraded to CONFIRMED NOVEL.** Five checks: 5+ sizes at q=5 (pairwise α converges upward: 2.076→2.099), cross-validation, falsification survived, exact decomposition formula, independence from entanglement spectrum.

[Full report: sprints/sprint_107.md]

### Sprint 108 — q=4 BKT Crossover: Log Corrections Map Three Regimes
**Status:** Complete (3 experiments).

**q=4 GPU extension (108a).** 7 sizes n=4-10 (GPU for n=9,10). Dominant fraction=1.000 at all sizes (single multiplet). Pairwise α converges downward: 1.825→1.771. Global α=1.788.

**Log correction analysis (108b).** Pairwise exponent drift vs 1/ln(N) for q=3,4,5 (7 sizes each). Walking (q=5): α_log = -0.003/ln(N), R²=0.001 — zero correction. BKT (q=4): α_log = +0.243/ln(N), R²=0.92 — moderate. Continuous (q=3): α_log = +0.460/ln(N), R²=0.98 — strong. Extrapolated α_∞: q=3→1.22, q=4→1.66, q=5→2.08.

**Multiplet structure (108c).** Gap ratio (multiplet/spectral) decreases with q: 6.2 (q=3) → 5.2 (q=4) → 4.6 (q=5). Spectral gap is (q-1)-fold degenerate and symmetry-forbidden. Dominant state is a singlet above the multiplet.

**Key finding:** Linear formula α=0.315q+0.469 valid ONLY for walking regime (q≥5) where log corrections vanish. For q≤4, finite-size α is inflated. Walking regime is a cleaner measurement regime than continuous or BKT. **NOTE:** 1/ln(N) extrapolation for q=2,3 RETRACTED in Sprint 109 — power-law 1/N² is correct.

[Full report: sprints/sprint_108.md]

### Sprint 109 — α(q) Correction Form: Power-Law Recovers Exact ν
**Status:** Complete (3 experiments).

**q=3 GPU extension (109a).** 9 sizes n=4-14 (GPU for n=12,14). Pairwise α converges to 1.416 at (n=12,14) — heading straight to exact 1.4 (ν=5/6). Single multiplet dominance (frac=1.000) confirmed at all sizes.

**q=2 and q=4 extension (109b).** q=2: 11 sizes n=4-18. Pairwise α converges to 1.012 at (n=16,18) — approaching exact 1.0 (ν=1). q=4: 8 sizes n=4-11 (GPU for n=11). Pairwise α flat at 1.771 for last 3 pairs — BKT slow convergence.

**Power-law vs log correction analysis (109c).** Power-law 1/N^ω fits vs 1/ln(N) fits. Power-law overwhelmingly wins for q=2 (AIC: -140 vs -87) and q=3 (AIC: -109 vs -71). Fitted ω≈2 — standard FSS correction, not leading irrelevant operator. Power-law α_∞: q=2→1.001 (exact 1.0, 0.1% error), q=3→1.405 (exact 1.4, 0.4% error). Log α_∞: q=2→0.874 (12.6% error), q=3→1.301 (7.1% error). Neither model works for q=4 BKT.

**Partial retraction of Sprint 108:** 1/ln(N) form was wrong for q=2,3 (continuous transitions have power-law, not log corrections). Extrapolated α_∞(q=3)=1.22 retracted → correct value 1.40. Walking (q≥5) zero-correction finding CONFIRMED.

[Full report: sprints/sprint_109.md]

### Sprint 110 — Extend α(q) to q=5-9: Formula Revision
**Status:** Complete (3 experiments).

**q=6 GPU extension (110a).** 4 sizes n=6-9 (GPU for n=9, dim=10M, 305s). Vectorized builder. Pairwise α converges upward: 2.358→2.380→2.402. Global α=2.377. Single-multiplet dominance (frac=1.000) confirmed.

**q=8,9 new measurements (110b).** q=8 n=6,7: α=2.897. q=9 n=6,7: α=3.159. Both below Sprint 103 formula prediction (2.99, 3.30). Single-multiplet dominance universal through q=9.

**α(q) formula refit (110c).** Sprint 103's α=0.315q+0.469 **REVISED** — was biased by including q=4 (BKT). Pure walking data (q=5-9): α ≈ 0.262q + 0.815 (linear, RMS=0.019) or α ≈ -0.009q²+0.389q+0.386 (quadratic, RMS=0.011, AIC-preferred ΔAIC=3.3). Bias-corrected: α ≈ 0.272q + 0.754. Component fits updated: z_m=0.082q+0.741, β_me=0.098q+0.333. Reconstructed α exactly matches direct fit.

**Partial revision of Sprint 103:** Formula coefficients updated but confirmed novel phenomenon (linear α(q) in walking regime) unchanged.

[Full report: sprints/sprint_110.md]

### Sprint 111 — q=10 χ_F Measurement: Power-Law α(q) Preferred
**Status:** Complete (2 experiments).

**q=10 spectral decomposition (111a).** 3 sizes n=5,6,7 (GPU for n=7, dim=10M, 700s). Pairwise α: (5,6)→3.295, (6,7)→3.417 — converging upward. Global α=3.349. Single-multiplet dominance (frac=1.000) at all sizes. z_m=1.566, β_me=1.285.

**α(q) refit with q=10 (111b).** 6 data points (q=5-10). AIC model comparison: power-law 0.69·q^0.69 preferred (ΔAIC=4.1 over linear). Linear 0.260q+0.827 (RMS=0.018) vs quadratic (RMS=0.013) vs power-law (RMS=0.013). Component fits stable: z_m=0.083q+0.734, β_me=0.094q+0.360. Reconstructed α matches direct linear fit exactly.

**Key findings:** (1) Single-multiplet dominance universal through q=10. (2) Power-law α(q) slightly preferred, but discrimination weak — all models fit comparably. (3) χ_F super-scaling (α>2) persists where walking has broken down (c_eff/Re(c)=0.60), confirming single-multiplet mechanism is independent of entropy-based walking.

[Full report: sprints/sprint_111.md]

### Sprint 112 — q=12 χ_F Measurement: Quadratic α(q) Strongly Preferred
**Status:** Complete (2 experiments).

**q=12 spectral decomposition (112a).** 3 sizes n=4,5,6 (GPU for n=6, dim=3.0M, 160s). Pairwise α: (4,5)→3.553, (5,6)→3.734 — converging upward. Global α=3.631. Single-multiplet dominance (frac=1.000) at all sizes. z_m=1.662, β_me=1.307. All three prior fits (linear, power-law, quadratic from Sprint 111) overshoot.

**α(q) refit with q=12 (112b).** 7 data points (q=5-12). AIC model comparison: **quadratic −0.010q²+0.41q+0.30 strongly preferred** (ΔAIC=12.2 over linear). Ranking: quadratic (0) > √q (+2.2) > log (+5.3) > power-law (+5.5) > linear (+12.2). Linear model RULED OUT. α(q) is sublinear — growth slows with q, approaching saturation near α≈4-5.

**Key findings:** (1) Linear α(q) from Sprint 110 definitively wrong at q≥12. (2) Sublinear growth — quadratic or √q form. (3) Single-multiplet dominance universal through q=12. (4) χ_F super-scaling persists at q=12 where walking is fully broken (c_eff/Re(c)≈0.51).

[Full report: sprints/sprint_112.md]

### Sprint 113 — DMRG Fidelity Susceptibility: Boundary Conditions Matter
**Status:** Complete (3 experiments).

**DMRG overlap χ_F validated (113a).** q=2 open BC, n=8-20 (DMRG extends to n=20). DMRG matches exact diag to 6 decimal places at all overlapping sizes (n=8-14). Open-BC pairwise α converges: 1.097→1.037. Global α=1.065, heading toward exact 1.0.

**q=4 and q=5 DMRG χ_F (113b).** Open-BC α dramatically different from periodic: q=4 α_open=1.51 (vs periodic 1.77), q=5 α_open=1.60 (vs periodic 2.09). Boundary effects dominate.

**Periodic vs open BC comparison (113c).** Overlap method matches spectral method exactly at periodic BC (q=4: 1.79 vs 1.77, q=5: 2.09 vs 2.09). BC effect on α: Δα=0.03 (q=2), 0.29 (q=4), 0.49 (q=5). **Periodic/open χ_F ratio diverges with N** (~N^{Δα}). Boundary fidelity susceptibility is enhanced at walking transitions.

**Key finding: Open-BC DMRG CANNOT extend α(q) to larger sizes.** Boundary effects for q≥4 are too large and growing. The q=4 open-BC α≈1.5 is coincidental, not evidence for BKT ν=2/3. Periodic exact diag remains the definitive measurement.

**POTENTIALLY NOVEL:** First systematic mapping of χ_F boundary condition dependence across walking boundary. Discovery that boundary fidelity susceptibility is enhanced (divergent ratio) at walking transitions.

[Full report: sprints/sprint_113.md]

### Sprint 114 — q=15 χ_F: Logarithmic α(q) Growth Emerging
**Status:** Complete (2 experiments).

**q=15 spectral decomposition (114a).** 3 sizes n=3,4,5 (GPU for n=5, dim=759k, 56s). Pairwise α: (3,4)→3.815, (4,5)→4.070 — converging upward. Global α=3.921. Single-multiplet dominance (frac=1.000) at all sizes. z_m=1.784, β_me=1.353. ALL prior models overpredict: quadratic by 0.28, √q by 0.40, log by 0.25, linear by 0.81.

**α(q) refit with 8 points q=5-15 (114b).** Quadratic overwhelmingly AIC-best (ΔAIC=14.9 over logarithmic). BUT quadratic −0.0134q²+0.451q+0.159 peaks at q≈17 then DECREASES — physically unreasonable. Logarithmic α(q) ≈ 1.73 ln(q) − 0.68 is the best physically unbounded model (RMS=0.041). Components also logarithmic: z_m ≈ 0.75 ln(q) − 0.21, β_me ≈ 1.31 ln(q) − 2.05.

**Key findings:** (1) α(q) sublinearity is stronger than any prior model predicted. (2) Logarithmic growth is the most physically plausible form — both components are individually logarithmic. (3) Single-multiplet dominance universal through q=15 where walking is fully broken.

[Full report: sprints/sprint_114.md]

### Sprint 115 — q=20 χ_F: Logarithmic α(q) Confirmed as AIC-Best
**Status:** Complete (2 experiments).

**q=20 spectral decomposition (115a).** 3 sizes n=3,4,5 (GPU for n=5, dim=3.2M, 333s). Pairwise α: (3,4)→4.416, (4,5)→4.835 — converging upward. Global α=4.590. Single-multiplet dominance (frac=1.000) at all sizes. z_m=2.047 (first time >2.0), β_me=1.496. Quadratic predicted 3.82 — off by 0.77 (DECISIVELY RULED OUT). Logarithmic predicted 4.49 — closest (off by 0.10).

**α(q) refit with 9 points q=5-20 (115b).** **Logarithmic NOW AIC-BEST** (ΔAIC=0, was +14.9 behind quadratic in Sprint 114). Updated: α(q) ≈ 1.78 ln(q) − 0.79 (RMS=0.046). Quadratic falls to ΔAIC=+9.4 and peaks at q≈24 (unphysical). All other models ΔAIC≥11.3. Components also logarithmic: z_m ≈ 0.76 ln(q) − 0.23, β_me ≈ 1.13 ln(q) − 1.68 (both AIC-best).

**Key findings:** (1) q=20 measurement flips AIC ranking — logarithmic goes from second-best to best. (2) Quadratic definitively ruled out (17% error at q=20). (3) z_m crosses 2.0 — multiplet gap closes faster than 1/N² (super-quartic χ_F). (4) Single-multiplet dominance universal through q=20 where walking is fully broken.

[Full report: sprints/sprint_115.md]

### Sprint 116 -- q=25 chi_F: Log+Loglog Subleading Correction Emerges
**Status:** Complete (2 experiments).

**q=25 spectral decomposition (116a).** 3 sizes n=3,4,5 (GPU for n=5, dim=9.8M, 1616s -- pushing GPU limit). Pairwise alpha: (3,4)->4.932, (4,5)->5.504 -- converging upward. Global alpha=5.170. Single-multiplet dominance (frac=1.000) at all sizes. z_m=2.286, beta_me=1.598. Pure log predicted 4.94 -- underpredicts by 4.7% (+0.23).

**alpha(q) refit with 10 points q=5-25 (116b).** **Log+loglog marginally AIC-best** (dAIC=1.4 over pure log): alpha(q) = 2.62 ln(q) - 1.77 ln(ln(q)) - 1.26. Pure log updated: 1.86 ln(q) - 0.96 (RMS=0.072). Quadratic peaks at q=32.6 (dAIC=+7.7, unphysical). Sqrt/power-law dAIC>=4.2. z_m(q) = 0.786 ln(q) - 0.290 (log AIC-best for components).

**Key findings:** (1) alpha(q=25)=5.17 confirms continued growth. (2) Subleading ln(ln(q)) correction improves fit but dAIC=1.4 is not decisive (need >4). (3) z_m=2.29 continues past 2.0. (4) n=5 at q=25 (9.8M) took 1616s -- practical GPU limit reached. (5) q=30 predictions: log+loglog 5.47, pure log 5.37.

[Full report: sprints/sprint_116.md]

### Sprint 117 -- q=30 chi_F: Logarithmic Model Stabilizes
**Status:** Complete (2 experiments).

**q=30 spectral decomposition (117a).** n=3,4 only (n=5=24.3M exceeds GPU). Pair(3,4) alpha=5.384. z_m=2.384. Single-multiplet dominance at q=30. Pure log predicted 5.37, log+loglog predicted 5.47.

**alpha(q) refit with 11 points (117b).** dAIC gap narrowed from 1.4 to 0.8. LOO cross-validation: log+loglog 0.081 vs pure log 0.085 -- nearly identical. **Pure logarithmic alpha(q) = 1.87 ln(q) - 0.97 is the preferred model** (Occam's razor). alpha(q) mapping q=5-30 now complete with diminishing returns.

[Full report: sprints/sprint_117.md]

### Sprint 118 -- q=4 chi_F Extended + Model Identity Audit
**Status:** Complete (1 experiment + audit).

**q=4 chi_F with 8 sizes n=4-11 (118a).** Pairwise alpha converges: 1.825->1.794->1.780->1.774->1.772->1.771->1.771. **alpha(q=4) = 1.771 +/- 0.001** (stable last 3 pairs). Differs from exact q=4 Potts alpha=2.0 by 11.5% — likely logarithmic corrections at the marginal Ashkin-Teller point.

**⚠ MODEL IDENTITY AUDIT:** Code audit confirmed ALL experiments from Sprint 076 onward use the **standard S_q Potts model** (Σ X^k field), NOT the "Potts-clock hybrid" (X+X† field) described in KNOWLEDGE.md. All six claimed novel findings are on the known S_q Potts chain. KNOWLEDGE.md, STATE.md, and sprint report corrected. The novel contributions are the probes (chi_F spectral decomposition, systematic scaling data), not the model itself.

~~Original Sprint 118 interpretation "confirms hybrid != S_q Potts" was WRONG~~ — the code IS S_q Potts.

[Full report: sprints/sprint_118.md]

### Sprint 119 — Hybrid Potts-Clock chi_F Spectral Decomposition
**Status:** Complete (2 experiments).

**Hybrid q=5 chi_F spectral (119a).** 7 sizes n=4-10 (GPU for n≥8). g_c=0.438. Single-multiplet dominance frac=1.000 at ALL sizes — identical to S_q. Spectral gap symmetry-forbidden. Global alpha=1.408, z_m=0.909, beta_me=0.590, nu_eff=0.831. **Prediction confirmed: alpha≈1.41, nu≈0.83.** Pairwise alpha DECREASING (1.56→1.26), opposite to S_q (INCREASING 2.08→2.10).

**Hybrid q=3,7,10 chi_F spectral (119b).** q=3 (sanity: hybrid=S_q): alpha=1.46→1.40, matches exact. q=7 (5 sizes n=4-8): alpha=0.952, z_m=0.770. q=10 (4 sizes n=4-7): alpha=0.552, z_m=0.688. Pairwise alpha STRONGLY decreasing at q≥7.

**Key findings:**
1. Single-multiplet dominance UNIVERSAL (both models, all q)
2. Decomposition alpha = beta_me + 2*z_m - 1 EXACT in both models
3. Hybrid alpha DECREASES with q (1.41→0.95→0.55) while S_q INCREASES (2.09→2.65→3.35)
4. z_m < 1 in hybrid (gap closes slower than 1/N) vs z_m > 1 in S_q (faster)
5. Walking super-scaling (alpha>2) is S_q-specific, NOT universal

**POTENTIALLY NOVEL:** First chi_F spectral decomposition on hybrid Potts-clock model. Universal mechanism with model-specific exponents. z_m crossing 1 as microscopic walking/continuous discriminator.

[Full report: sprints/sprint_119.md]

### Sprint 120 — Hybrid q=4 chi_F: z_m=1 Transition Point
**Status:** Complete (2 experiments).

**Hybrid q=4 g_c scan (120a).** Gap×N crossings at n=6,8,10 (GPU). Coarse scan g∈[0.33,0.43], 30 points. **g_c(hybrid,q=4) = 0.3933 ± 0.0002.** Excellent convergence: (6,8)→0.3934, (8,10)→0.3933.

**Hybrid q=4 chi_F spectral (120b).** 8 sizes n=4-11 (up to dim=4.2M, GPU). Plus q=3 sanity check. Single-multiplet dominance frac=1.000 at ALL sizes. **Global alpha=1.548, z_m=1.004, beta_me=0.541, nu_eff=0.785.** Pairwise alpha DECREASING (1.64→1.49), opposite to S_q q=4 (expected upward toward 2.0).

**Key findings:**
1. **z_m=1.004 at q=4** — exactly marginal. Hybrid z_m crosses 1 between q=3 (1.03) and q=5 (0.91), with q=4 at the boundary
2. **Hybrid alpha peaks at q=4** (1.55), then declines: alpha(q)=[1.40, 1.55, 1.41, 0.95, 0.55] for q=[3,4,5,7,10]. Non-monotonic
3. **Model divergence starts at q=4** — 12% gap (hybrid 1.55 vs S_q 1.77), growing to 64% at q=7
4. Pairwise alpha drift is OPPOSITE between models at q=4: hybrid downward, S_q upward. Different universality classes
5. Literature confirms S_q q=4 has log corrections from marginal operator (Salas-Sokal 1997): alpha→2.0 asymptotically

[Full report: sprints/sprint_120.md]

### Sprint 121 — z_m(q) Continuous Fit & S_q q=4 Extended Sizes
**Status:** Complete (2 experiments).

**z_m(q) fitting (121a).** 5 functional forms (linear, log, 1/√q, rational, power-law) fitted to hybrid z_m data at q=[3,4,5,7,10]. All converge: **q_cross = 3.58 ± 0.04.** Best fit: rational z_m = 1.368/(1+0.102q), R²=0.973. Alpha peak also at q≈3.58 from quadratic-in-ln(q) fit.

**S_q q=4 extended (121b).** 8 sizes n=4-11 (up to dim=4.2M, GPU). Single-multiplet dominance frac=1.000. **Global alpha=1.777±0.003, z_m=1.092±0.001, beta_me=0.604±0.002.** Pairwise alpha converges from above: 1.825→1.771. S_q z_m slowly drifting down: 1.097→1.079.

**Key findings:**
1. **q_cross = 3.58 ± 0.04** — walking→continuous boundary at non-integer q. Hybrid q=4 is barely walking (z_m=1.004), not exactly marginal
2. **S_q q=4: alpha=1.777, z_m=1.092** — firmly walking, consistent with asymptotic alpha=2.0 plus log corrections
3. **Alpha drift comparison:** S_q converges to ~1.77 (stable), hybrid decreasing from 1.64→1.49 (still drifting). 19% divergence at (10,11)
4. S_q z_m drift (1.097→1.079) suggests logarithmic approach to z_m=1 at larger sizes — consistent with BKT marginal operator

[Full report: sprints/sprint_121.md]
