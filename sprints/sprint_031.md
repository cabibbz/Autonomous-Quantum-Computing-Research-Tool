# Sprint 031 — Entanglement Spectrum: What the Reduced Density Matrix Hides Beyond Entropy

**Date:** 2026-04-01
**Status:** In Progress

## Motivation

We've characterized entanglement through scalar measures (entropy, MI, I3, discord, negativity) and discovered 5 archetypes. But the reduced density matrix contains far more information than any single number can capture — its full eigenvalue spectrum (the "entanglement spectrum") encodes edge modes, topological order, and universality classes.

The Li-Haldane conjecture (2008) posits that the entanglement spectrum mirrors the physical edge spectrum of topological phases. A 2025 paper showed this correspondence can *break down* in XXZ chains — entanglement phase transitions can be disconnected from bulk transitions. This directly connects to our Sprint 030 finding that archetype boundaries ≠ phase boundaries (I3 sign change at Δ≈0.7 inside XY phase).

**Key question:** Does the entanglement spectrum reveal structure invisible to scalar entanglement measures? Can it explain the archetype-phase boundary mismatch?

**Literature searched:** arXiv 2509.03588 (entanglement phase transitions in Haldane phase, 2025), Nature Comms 2023 (energy-entanglement spectrum wormhole connection), Quantum 2022 (Li-Haldane conjecture probing).

## Experiments

### 31a: Entanglement Spectrum of 5 Archetypes
**Goal:** Compute full eigenvalue spectrum of half-cut RDM for GHZ, W, Cluster 1D, Cluster 2D, and TFIM critical (Scale-Free). Characterize spectrum shape: flat, peaked, power-law, gapped.

### 31b: Entanglement Spectrum Across TFIM Phase Diagram
**Goal:** Track spectrum evolution through ordered → critical → disordered. Does the entanglement gap close at criticality? Does spectrum shape encode universality?

### 31c: Entanglement Spectrum Across XXZ Phase Diagram
**Goal:** Track spectrum across FM → XY → Néel. Does I3 sign change at Δ≈0.7 correspond to a spectrum transition? Test the "dissociation" of entanglement and bulk transitions.
