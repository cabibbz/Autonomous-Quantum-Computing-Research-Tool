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

**Key finding:** Linear formula α=0.315q+0.469 valid ONLY for walking regime (q≥5) where log corrections vanish. For q≤4, finite-size α is inflated. Walking regime is a cleaner measurement regime than continuous or BKT.

[Full report: sprints/sprint_108.md]
