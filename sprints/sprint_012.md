# Sprint 012 — Quantum Scrambling: How Fast Does Information Spread?

**Date:** 2026-03-31
**Status:** Complete (3/3 experiments)

## Motivation

Sprint 011 showed that entangling gates restructure correlations (pairwise → multipartite). Quantum scrambling is the extreme limit of this process: local information becomes delocalized across the entire system, unrecoverable from any local measurement. Scrambling connects to black hole physics (Hayden-Preskill protocol), quantum computational advantage (random circuit sampling), and quantum chaos (OTOCs as quantum Lyapunov exponents).

## Experiments & Results

### 12a: Random Circuit Scrambling Dynamics

**Setup:** n=6 qubits, random 2-qubit gates in brick-wall layers. Two geometries: 1D nearest-neighbor vs all-to-all. 10 random instances, 0–15 layers. Tracked half-cut entropy, total pairwise MI, and average I3.

**Results:**
| Metric | 1D (layer 5) | All-to-all (layer 5) | Both (layer 15) |
|--------|-------------|---------------------|-----------------|
| Half-cut entropy | 1.72 | 2.12 | ~2.3 (79% Page) |
| Total MI | 2.31 | 2.03 | ~1.5 |
| Average I3 | -0.227 | -0.253 | ~-0.30 |

**Key findings:**
- All-to-all scrambles ~2x faster initially but both converge by layer ~10
- Half-cut entropy plateaus at ~80% of Page entropy (2.91) — not reaching the random state limit even at 15 layers
- Total pairwise MI *decreases* as scrambling proceeds — same correlation restructuring seen in Sprint 011 but now continuous
- I3 converges to ≈ -0.30 for both geometries — strong, geometry-independent multipartite correlations in the scrambled state

### 12b: OTOCs — Out-of-Time-Ordered Correlators

**Setup:** OTOC F(t) = |⟨ψ|W†(t)V†W(t)V|ψ⟩|² with W=Z₀ (perturbation on qubit 0), V=Zⱼ (probe at distance j). Measures how fast information from qubit 0 reaches qubit j.

**Results — 1D brick-wall (scrambling time t_scr = first layer with F < 0.5):**
| Probe qubit | Distance | t_scr | F(12) |
|-------------|----------|-------|-------|
| Z₁ | 1 | 1 | 0.012 |
| Z₂ | 2 | 3 | 0.011 |
| Z₃ | 3 | 3 | 0.018 |
| Z₄ | 4 | 5 | 0.004 |
| Z₅ | 5 | 5 | 0.012 |

**Results — All-to-all:**
| Probe qubit | Distance | t_scr | F(12) |
|-------------|----------|-------|-------|
| All qubits | any | 2–3 | ~0.02 |

**Key findings:**
- **1D has a clear light cone**: t_scr scales linearly with distance (~1 layer per qubit). This is a Lieb-Robinson bound in action — information propagates at a finite "speed of light" in 1D.
- **All-to-all has no light cone**: uniform t_scr ≈ 2–3 regardless of distance. Long-range connectivity allows instant delocalization — this is "fast scrambling" (saturates the Maldacena-Shenker-Stanford bound).
- All OTOCs decay to ~0.01–0.03 by layer 12 — complete scrambling in both cases, just at different rates.

### 12c: Scrambling Efficiency — Structured vs Random

**Setup:** Compare GHZ, W, Cluster preparation circuits to random circuits. Then add random layers post-preparation and track convergence.

**Preparation circuits alone (5 two-qubit gates each):**
| State | Half-cut S | S/gate | Total MI | OTOC |
|-------|-----------|--------|----------|------|
| GHZ | 1.00 (33%) | 0.200 | 15.0 | **1.000** |
| W | 0.97 (32%) | 0.194 | 2.38 | **1.000** |
| Cluster | 1.00 (33%) | 0.200 | 2.00 | **1.000** |

**Post-preparation (adding 10 random layers):**
| Starting state | S(+0) | S(+10) | OTOC(+10) |
|---------------|-------|--------|-----------|
| GHZ | 1.00 | 2.26 | 0.017 |
| W | 0.97 | 2.25 | 0.022 |
| Cluster | 1.00 | 2.29 | 0.021 |
| Random (no prep) | 0.00 | 2.13 | 0.023 |

**Key findings:**
- **All structured states have OTOC = 1.0** — they create entanglement but do NOT scramble! This is the central surprise of the sprint.
- GHZ has enormous MI (15.0 — every pair is maximally correlated) but zero scrambling. Entanglement ≠ scrambling.
- Post-preparation scrambling converges identically regardless of starting state (~S=2.25 at +10 layers)
- Structured states give a ~1 bit head-start in entropy over random circuits but NO advantage in scrambling speed
- None reach 90% of Page entropy — scrambling to the random-state limit is hard and slow

## Analysis & Key Insights

### 1. Entanglement ≠ Scrambling
This is the headline result. GHZ creates maximal half-cut entanglement and maximal pairwise MI, yet its OTOC is exactly 1.0 — a local Z perturbation on qubit 0 remains perfectly localized. Why? Because GHZ correlations are *structured*: every qubit carries the same bit of information (|00...0⟩ vs |11...1⟩). Scrambling requires information to be *delocalized* — spread across many qubits such that no local measurement can recover it. Structured entanglement concentrates information; scrambling distributes it.

### 2. The Information Light Cone
1D random circuits show a clear Lieb-Robinson-like bound: information propagates at ~1 qubit per layer. This is the quantum analog of the speed of light — a fundamental limit on how fast quantum information can spread in geometrically local systems. All-to-all connectivity bypasses this, achieving "fast scrambling" in O(log n) depth — the theoretical minimum (Maldacena-Shenker-Stanford bound).

### 3. The Scrambling Hierarchy
Combining Sprint 011 and Sprint 012:
- **Layer 1–3**: Entanglement creation (S increases, MI increases)
- **Layer 3–8**: Correlation restructuring (S increases, MI *decreases* as pairwise → multipartite)
- **Layer 8–15**: Approach to scrambled state (S, MI, I3 all plateau near universal values)

The "scrambled state" is characterized by: S ≈ 80% Page entropy, MI ≈ 1.5 (mostly multipartite), I3 ≈ -0.30 (strong irreducible correlations).

### 4. Universality
Perhaps the deepest finding: scrambled states are *universal*. Regardless of starting point (GHZ, W, Cluster, or |0⟩), random circuits converge to the same statistical ensemble. The preparation circuit's structure is completely erased. This is the quantum circuit analog of thermalization — the second law of quantum information.

## Surprises
1. **GHZ has OTOC = 1.0 despite MI = 15.0** — maximal correlations, zero scrambling. Entanglement and scrambling are fundamentally different.
2. **Clear light cone in 1D**: t_scr = {1, 3, 3, 5, 5} exactly mirrors the brick-wall geometry (2 layers per "brick width")
3. **All-to-all achieves fast scrambling** (t_scr ≈ 3, distance-independent) — saturating the MSS bound for n=6
4. **Structured states DON'T help scrambling** — post-preparation convergence is identical to starting from |0⟩
5. **80% Page entropy barrier** — 15 layers of random gates can't reach the random state limit for n=6. Full scrambling requires exponential depth?

## Connections
- **Black holes**: Hayden-Preskill says you can recover information thrown into a black hole after the scrambling time. Our OTOC measurements directly give this time.
- **Quantum advantage**: Google's random circuit sampling works because circuits quickly reach a hard-to-simulate distribution. Our 80% Page barrier suggests there's a gap between "hard to simulate" and "fully random."
- **Sprint 011**: The MI-destruction during circuit construction is the *beginning* of scrambling. Sprint 012 shows the full trajectory: create → restructure → scramble.
- **Sprint 010**: The body-order hierarchy (W=1-body, GHZ=2-body, Cluster=3-body) maps onto scrambling: higher body-order doesn't mean more scrambled.

## Next Ideas
- **Scrambling vs circuit complexity**: minimum circuit depth to prepare a state vs time to scramble it — are they related?
- **Hayden-Preskill protocol**: simulate throwing a qubit into a "black hole" (scrambling circuit) and recovering it from the radiation
- **Real hardware**: run OTOCs on IBM hardware — noise should HELP scrambling (decoherence is a free scrambler)
- **Quantum error correction**: how does scrambling interact with error correction codes? Does the stabilizer structure resist scrambling?
