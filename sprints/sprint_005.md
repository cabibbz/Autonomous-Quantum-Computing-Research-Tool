# Sprint 005 — Quantum Discord: Separating Quantum from Classical Correlations

**Date:** 2026-03-31
**Status:** Complete

## Motivation

We've built an information-theoretic hierarchy: entropy < single-qubit entropy < MI < I3. Each level reveals structure the previous misses. The next natural measure is **quantum discord** — it decomposes mutual information into quantum and classical parts.

Key question: GHZ has positive I3 (redundant, "classical-like"). Does it have low discord too? If yes, GHZ correlations are genuinely classical in character. If discord is high, then I3 is misleading about classicality.

## Experiments

### 5a: Pairwise Quantum Discord (n=6)

Discord D(A|B) = MI(A:B) - J(A:B), where J is classical mutual information maximized over local measurements on B.

**Results:**

| State | Mean MI | Mean Discord | Quantum Fraction |
|-------|---------|-------------|-----------------|
| GHZ | 1.000 | 0.000 | 0% |
| W | 0.382 | 0.282 | 74% |
| Cluster | 0.133 | 0.000 | 0% |

**GHZ: 100% classical pairwise.** The 2-qubit RDM is (|00><00| + |11><11|)/2 — a classical mixture. Zero entanglement, zero quantum correlations in any pair. All correlations are perfectly captured by measuring in the Z basis.

**W: 74% quantum.** The most quantum state at the pairwise level. The W state's 2-qubit RDM has coherences (superposition of |01> and |10>) that no classical measurement can capture.

**Cluster: zero discord, but only edge pairs have any MI at all.** Interior pairs (1,2), (2,3), etc. have MI = 0. Edge pairs (0,1) and (4,5) have MI = 1.0 but discord = 0.

### 5b: Discord at Different Scales — D(A|BC)

Compute discord between one qubit and a pair, using local measurements on the pair. Tests whether multipartite quantum correlations (visible in I3) appear as discord at a larger scale.

**Results for consecutive triples:**

| State | MI(A:BC) | Discord D(A|BC) | Classical Corr | Quantum Fraction |
|-------|----------|-----------------|---------------|-----------------|
| GHZ | 1.000 | 0.000 | 1.000 | 0% |
| W | 0.568 | 0.374 | 0.194 | 66% |
| Cluster | 1.000 | 0.000 | 1.000 | 0% |

**Cluster state: zero discord even at 1-vs-2 level.** Despite having I3 = -1.0 (genuinely irreducible 3-body correlations), all correlations can be extracted by local classical measurements.

## Key Insights

### 1. The Information Hierarchy is Now Complete

entropy < single-qubit entropy < MI < I3 < discord (different axis)

Discord and I3 measure orthogonal things:
- I3 measures whether correlations are pairwise-decomposable (structure)
- Discord measures whether correlations are classically accessible (nature)

### 2. Three Archetypes of Multipartite Entanglement

| Property | GHZ | W | Cluster |
|----------|-----|---|---------|
| Pairwise MI | 1.0 (uniform) | 0.38 (uniform) | 0/1 (edges) |
| Discord | 0 | 0.28 (74%) | 0 |
| I3 | +1.0 | ? | -1.0 |
| Character | Classical broadcast | Quantum sharing | Structural quantum |

- **GHZ** = "classical broadcast in superposition" — all pairs maximally correlated, all correlations classical, globally redundant
- **W** = "quantum sharing" — weak but genuinely quantum pairwise correlations, the only state with nonzero discord
- **Cluster** = "structural quantum" — no pairwise quantum correlations at all, quantumness is in the pattern not the parts

### 3. Connection to Measurement-Based Quantum Computing

The cluster state's zero discord at all pairwise levels, combined with its utility for universal quantum computation via local measurements, reveals something deep: **the quantum computational resource in cluster states is not in any individual correlation but in the global structure of correlations.** You compute by making classical measurements — the quantum power comes from the topology of the correlation network, not from any particular quantum correlation.

## What's Next

- **Compute I3 for W state** to complete the archetype table
- **Quantum discord under noise** — how quickly do W state's quantum correlations degrade?
- **Entanglement measures** (concurrence, negativity) to compare with discord
- **2D cluster states** — does geometry change the discord structure?
- **Integrated Information (Phi)** — a measure of system integration, complementary to discord
