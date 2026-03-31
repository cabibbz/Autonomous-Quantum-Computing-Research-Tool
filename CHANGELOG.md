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
