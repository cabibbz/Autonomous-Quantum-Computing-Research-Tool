# Sprint 022 — Noise-Adapted Codes: When Does Bias-Awareness Beat Isotropy?

**Date:** 2026-03-31
**Status:** In Progress

## Motivation

Sprint 021 established that basis isotropy is the defining property of good QEC codes — [[5,1,3]] dominates universally because it protects all Bloch sphere directions equally. But real superconducting hardware has strongly biased noise: T2 << T1 means dephasing (Z errors) dominate over relaxation (X/Y errors). The XZZX surface code literature shows that bias-tailored codes can achieve dramatically higher thresholds at large scale.

**The gap:** At small scale (≤10 qubits), when does a bias-aware code choice beat the isotropic [[5,1,3]]? Is there a bias ratio where specialization wins?

## Literature Search

- XZZX surface code: rotated stabilizers exploit Z-bias, achieving higher thresholds under biased noise
- Bias-tailored LDPC codes (arXiv:2202.01702): framework for bias-tailoring beyond 2D topological codes
- Google's below-threshold surface code (Nature 2024): real hardware demonstration
- Key insight from literature: bias-tailoring is well-studied at large scale but small-scale crossover analysis is sparse

## Codes Under Test

| Code | Qubits | Stabilizers | Protects | Bias |
|------|--------|-------------|----------|------|
| 3-qubit bit-flip | 3 | XXI, IXX | X errors | X-specialized |
| 3-qubit phase-flip | 3 | ZZI, IZZ | Z errors | Z-specialized |
| [[5,1,3]] | 5 | XZZXI, IXZZX, XIXZZ, ZXIXZ | All single-qubit | Isotropic |
| Uncoded | 1 | — | Nothing | — |

## Experiments

### 22a: Specialized codes under biased noise
**Question:** Under Z-biased noise (T2 << T1), does the phase-flip code outperform [[5,1,3]]?

### 22b: Optimal code map across (γ, λ) landscape
**Question:** For each noise point, which code wins on basis-averaged Holevo?

### 22c: Cost of isotropy — efficiency under known bias
**Question:** How much performance does [[5,1,3]] leave on the table when the bias is known?

---

## Results

*(To be filled as experiments complete)*
