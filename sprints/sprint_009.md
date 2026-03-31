# Sprint 009 — Structured Noise: How Different Error Channels Shape Entanglement

**Date:** 2026-03-31
**Status:** In Progress

## Motivation

All prior noise experiments used depolarizing noise — symmetric, uniform, unrealistic. Real quantum hardware has *structured* noise: amplitude damping (T1 energy relaxation, |1⟩→|0⟩), phase damping (T2 dephasing, loss of coherence without energy loss). These channels break different symmetries and should discriminate our state archetypes differently.

Phase damping is the dominant error in superconducting qubits (T2 < T1 always). Understanding which states survive which noise channels connects directly to quantum error correction design.

## Experiments

### 9a: Amplitude Damping — Energy Relaxation
**Question:** How does T1-like noise (each qubit independently decays |1⟩→|0⟩) degrade entanglement across GHZ, W, Cluster 1D, Cluster 2D?

### 9b: Phase Damping — Dephasing
**Question:** How does T2-like noise (each qubit independently loses phase coherence) compare? Phase damping preserves populations but kills superpositions.

### 9c: Noise Channel Fingerprint
**Question:** At equal noise strength, which states are robust to which noise type? Build a noise-state discrimination matrix.

---

## Results

*(updated incrementally)*
