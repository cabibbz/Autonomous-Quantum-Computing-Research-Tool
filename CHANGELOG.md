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
