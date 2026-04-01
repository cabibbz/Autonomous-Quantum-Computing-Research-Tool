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
