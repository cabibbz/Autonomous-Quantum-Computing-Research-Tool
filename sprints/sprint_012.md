# Sprint 012 — Quantum Scrambling: How Fast Does Information Spread?

**Date:** 2026-03-31
**Status:** In Progress

## Motivation

Sprint 011 showed that entangling gates don't just create correlations — they restructure them. Random quantum circuits are the canonical model for quantum scrambling, where local information becomes delocalized across the entire system. Scrambling connects to:

- **Black hole information paradox** — Hayden-Preskill showed that scrambling determines when information can be recovered from a black hole
- **Quantum computational advantage** — random circuit sampling (Google's quantum supremacy) relies on circuits being good scramblers
- **Quantum chaos** — OTOCs (out-of-time-ordered correlators) are the quantum analog of Lyapunov exponents

Key questions:
1. How fast do random circuits scramble information compared to our structured states?
2. Does circuit geometry (1D nearest-neighbor vs all-to-all) affect scrambling speed?
3. How do our archetypes (GHZ, W, Cluster) compare to random states as scramblers?

## Experiments

### 12a: Random Circuit Scrambling Dynamics
Track entropy, MI, and I3 as layers of random 2-qubit gates are applied. Compare 1D (nearest-neighbor, brick-wall) vs all-to-all connectivity.

### 12b: OTOCs — Out-of-Time-Ordered Correlators
Compute OTOCs F(t) = |⟨ψ|W†(t)V†W(t)V|ψ⟩|² to measure how a local perturbation spreads. Track decay as a function of circuit depth.

### 12c: Scrambling Efficiency of Structured vs Random Circuits
Compare how quickly GHZ, W, Cluster preparation circuits scramble vs random circuits of the same depth. Are structured circuits more or less efficient scramblers?

---

## Results

### 12a: Random Circuit Scrambling Dynamics

*(pending)*

### 12b: OTOCs

*(pending)*

### 12c: Scrambling Efficiency

*(pending)*
