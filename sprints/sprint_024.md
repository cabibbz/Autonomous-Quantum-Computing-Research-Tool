# Sprint 024 — Surface Code vs Shor: Same [[9,1,3]], Different Architecture

**Date:** 2026-03-31
**Status:** In Progress

## Motivation

Sprint 023 showed that the figure of merit determines the winner: avg Holevo → Shor wins, min-basis Holevo → [[5,1,3]] wins. But we haven't tested the most practically important code family: the **surface code**.

The rotated surface code at distance d=3 is [[9,1,3]] — the **exact same parameters** as the Shor code. Same number of physical qubits (9), same number of logical qubits (1), same code distance (3). But completely different architecture:
- **Shor:** Concatenated (phase-flip code inside bit-flip code)
- **Surface:** Topological (2D lattice with boundary stabilizers)

This is the cleanest possible architectural comparison. No qubit-count or distance confounds.

**Literature search:** arXiv/Google Scholar — no prior work directly compares surface [[9,1,3]] vs Shor [[9,1,3]] using Holevo information or basis-averaged metrics. Papers on surface code noise thresholds exist (arXiv:2403.08706, Nature 2024) but focus on large-scale logical error rates, not small-scale information-theoretic comparison.

## Experiments

### 24a: Surface Code Entanglement Structure
- Build rotated d=3 surface code stabilizers
- Compute MI, I3, negativity spectrum
- Compare to Shor [[9,1,3]] and [[5,1,3]]

### 24b: Surface vs Shor Under Combined T1+T2 Noise
- Head-to-head basis-averaged Holevo across (γ,λ) grid
- Which architecture wins under realistic noise?

### 24c: All-Code Min-Basis Holevo Ranking
- Comprehensive ranking: surface, Shor, [[5,1,3]], Steane, 3-qubit, uncoded
- The definitive small-code comparison

---

## Results

*(to be filled as experiments complete)*
