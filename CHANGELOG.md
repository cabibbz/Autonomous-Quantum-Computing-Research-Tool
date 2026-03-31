# Changelog — Quantum Explorer

## QPU Budget
- Monthly allocation: 600 seconds
- Used this period: 0s
- Period start: 2026-03-31

## Constraints (learned the hard way)
- **Max qubits for partial_trace: 10.** Beyond that, CPU time explodes.
- **Max bash command runtime: 60 seconds.** Design experiments to fit.
- **Save results after EACH experiment.** Never batch.
- **Separate scripts per experiment.** One timeout shouldn't kill everything.

---

## Sprint Log

### Sprint 001 — 2026-03-31 — Bell States & CHSH
**Status:** Partial (2/4 experiments; GHZ scaling timed out)

**Completed:**
- Verified all 4 Bell states: perfect correlations, entropy = 1.0 bit
- CHSH inequality: **S = 2.834** (violates classical bound of 2.0, matches Tsirelson bound 2.828)

**Failed:**
- GHZ entropy scaling (n>15 qubits): partial_trace too slow on CPU
- Random vs GHZ entropy: same issue

**Surprises:**
- CHSH angles are subtle: wrong Bob angles give S=0 even though individual correlations look correct
- Simulator is very accurate for small circuits (S matches theory to 3 decimals at 100k shots)

**Next:** Map CHSH violation landscape, explore noise effects, small-scale entanglement comparison

[Full report: sprints/sprint_001.md]

### Sprint 002 — 2026-03-31 — CHSH Landscape, Noise, & Entanglement Zoo
**Status:** Complete (3/3 experiments)

**Completed:**
- **2a: CHSH angle landscape** — only 19.9% of measurement angles violate classical bound. Max S=2.828 at known optimal.
- **2b: Noise sensitivity** — depolarizing noise kills CHSH violation at **p≈9.5%**. S degrades exponentially.
- **2c: Entanglement zoo** — GHZ, W, Cluster, Product states at n=4,6,8. Half-cut entropy is 1.0 for ALL entangled states (can't distinguish them). But single-qubit and qubit-loss metrics reveal real differences.

**Surprises:**
- **Cluster states get MORE entangled when you lose a qubit** (entropy 1.0 -> 2.0). Removing a qubit projects neighbors into entangled pairs. No classical analog.
- W state entanglement "dilutes" with system size — single-qubit entropy drops from 0.81 (n=4) to 0.54 (n=8)
- Half-cut entropy is too coarse to distinguish qualitatively different entanglement structures

**Next:** Deeper investigation of cluster state qubit-loss phenomenon. Entanglement robustness under noise by state type. Try real IBM hardware for CHSH.

[Full report: sprints/sprint_002.md]

### Sprint 003 — 2026-03-31 — Cluster State Entanglement Under Qubit Loss
**Status:** Complete (2/2 experiments)

**Completed:**
- **3a: Progressive loss** — GHZ stays rigid (1.0 always), W decays monotonically, Cluster jumps to 2.0 and plateaus before dropping
- **3b: Position-dependent loss** — Sharp boundary at the half-cut: losing qubits 0-2 gives entropy 2.0, qubits 3-7 gives 1.0. Two-qubit loss can reach **entropy 3.0** (triple the original!)

**Surprises:**
- Entropy increase from qubit loss depends on the qubit's position relative to the measurement bipartition, not on the qubit's role in the state
- Two-qubit loss can unlock **3 bits** of entanglement that were "hidden" in the original state
- GHZ is astonishingly rigid: lose 6 of 8 qubits, still exactly 1 bit of entanglement

**Key insight:** Cluster states contain "trapped" entanglement that is released by qubit loss. This connects to measurement-based quantum computing — computation by entanglement destruction.

**Next:** 2D cluster states, projective measurement vs trace, mutual information, IIT/Phi

[Full report: sprints/sprint_003.md]

### Sprint 004 — 2026-03-31 — Mutual Information & Tripartite Information
**Status:** Complete (2/2 experiments)

**Completed:**
- **4a: Pairwise mutual information** — Finally distinguishes all three states. GHZ: uniform maximal MI. W: uniform weak MI. Cluster: sparse, only nearest-neighbor MI.
- **4b: Tripartite information** — GHZ has I3=+1.0 (redundant/classical-like). Cluster has I3 down to **-1.0** for consecutive triples — genuinely irreducible 3-body quantum correlations.

**Surprises:**
- GHZ's correlations are structurally classical (positive I3) despite being "maximally entangled" — it's a broadcast channel in superposition
- Cluster state geometry is directly visible in I3: only consecutive-on-chain triples show negative I3
- Hierarchy of measures: entropy < single-qubit entropy < MI < I3 in discriminating power

**Key insight:** Negative tripartite information = information spread across 3 subsystems that no pair can recover. This is the operational meaning of genuine multipartite entanglement and connects to quantum error correction.

**Next:** Integrated information (Phi), quantum discord, real hardware test

[Full report: sprints/sprint_004.md]

### Sprint 005 — 2026-03-31 — Quantum Discord: Classical vs Quantum Correlations
**Status:** Complete (2/2 experiments)

**Completed:**
- **5a: Pairwise quantum discord** — GHZ has ZERO discord (100% classical pairwise). W has 74% quantum fraction — most quantum state. Cluster has zero discord for all pairs.
- **5b: Discord at 1-vs-2 scale** — Cluster state has zero discord even measuring one qubit against a pair. W discord increases to 0.37. GHZ stays zero.

**Surprises:**
- GHZ is "maximally entangled" but has ZERO quantum correlations in any pair — its 2-qubit RDM is a classical mixture
- Cluster states have zero discord at ALL pairwise scales despite having I3=-1.0 (irreducibly multipartite correlations)
- W state is the ONLY one with nonzero discord — it's the most "quantum" at the pairwise level despite being less entangled than GHZ

**Key insight:** Discord and I3 are orthogonal measures. I3 measures correlation *structure* (decomposable vs irreducible). Discord measures correlation *nature* (classical vs quantum). Cluster states are "structural quantum" — quantumness is in the pattern of correlations, not in any individual correlation. This explains why MBQC works with local classical measurements.

**Three archetypes confirmed:**
- GHZ = classical broadcast in superposition (MI=1, discord=0, I3=+1)
- W = quantum sharing (MI=0.38, discord=0.28, 74% quantum)
- Cluster = structural quantum (MI sparse, discord=0, I3=-1)

**Next:** I3 for W state (complete archetype table), discord under noise, concurrence/negativity, 2D cluster states, Integrated Information (Phi)

[Full report: sprints/sprint_005.md]

### Sprint 006 — 2026-03-31 — Entanglement Monotones: Concurrence, Negativity & Monogamy
**Status:** Complete (3/3 experiments)

**Completed:**
- **6a: Pairwise concurrence & negativity** — GHZ: zero (all multipartite). W: uniform C=2/n across all pairs. Cluster: zero (all multipartite, like GHZ).
- **6b: Noise sensitivity** — Half-cut negativity decays identically for all three under depolarizing noise (non-discriminating). W pairwise concurrence dies at p≈14%.
- **6c: Entanglement spectrum** — Negativity across ALL bipartitions. GHZ: flat (0.5 everywhere). W: grows with cut size (0.37→0.50). Cluster: explosive (0.5 to 3.5, geometry-dependent).

**Surprises:**
- Cluster state can have **7x more entanglement** than GHZ at certain bipartitions — "entanglement hotspots"
- W state has **zero residual tangle** (CKW monogamy saturated) — ALL its entanglement is pairwise, zero genuine tripartite
- GHZ and W are opposite extremes of monogamy: GHZ = all multipartite, W = all pairwise
- Entanglement spectrum is the most powerful single discriminator found so far

**Key insight:** Three distinct entanglement *topologies*: Democratic (GHZ, flat spectrum), Distributed (W, growing spectrum), Geometric (Cluster, explosive spectrum). The topology determines computational utility — Cluster's geometry-dependent entanglement is why MBQC works.

**Next:** I3 for W state, 2D cluster states, entanglement spectrum under qubit loss, local noise models, GME witnesses

[Full report: sprints/sprint_006.md]

### Sprint 007 — 2026-03-31 — 2D Cluster States: Does Geometry Shape Entanglement?
**Status:** Complete (3/3 experiments)

**Completed:**
- **7a: W-state I3 + 2D cluster MI & I3** — W-state I3 = +0.195 (weak redundant, like GHZ). 2D cluster (2x3) has 8 negative-I3 triples vs 1D's 6 — geometry amplifies irreducible correlations.
- **7b: Entanglement spectrum** — 2D eliminates weak bipartitions. ALL |A|=2 cuts have negativity 1.5 (1D: 0.5–1.5). 2D raises the entanglement floor — "entanglement democracy."
- **7c: 2D qubit loss** — 2D baseline entropy = 3.0 (vs 1D = 1.0). Loss is completely position-independent in 2D (always 2.0). Second qubit loss causes no additional degradation. Opposite behavior to 1D (entropy decreases vs increases).

**Surprises:**
- 2D cluster starts with 3x the half-cut entanglement of 1D
- Qubit loss in 2D is **position-independent** — corner and edge behave identically. In 1D, position matters enormously.
- **Second qubit loss is free** in 2D — no additional entropy degradation
- 1D and 2D have **opposite** loss behavior: 1D releases trapped entanglement (entropy up), 2D gracefully degrades (entropy down)

**Key insight:** 1D→2D is a qualitative transition from geometric to **topological** entanglement. 2D achieves uniform, position-independent robustness — the foundation of topological quantum error correction (surface codes).

**Next:** Integrated Information Theory (Phi), local noise models, squarer 2D geometries (3x3)

[Full report: sprints/sprint_007.md]
