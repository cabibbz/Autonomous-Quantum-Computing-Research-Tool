# Sprint 022 — Noise-Adapted Codes: When Does Bias-Awareness Beat Isotropy?

**Date:** 2026-03-31
**Status:** Complete (3/3 experiments)

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

### 22a: Specialized codes under biased noise

**Setup:** Fixed γ=0.05 (mild T1), swept λ from 0 to 0.50 (Z-bias ratio up to 10:1). Also reverse: fixed λ=0.05, swept γ (X-bias). Tested Z, X, Y logical bases. Basis-averaged Holevo as metric.

**Result: [[5,1,3]] wins at ALL bias ratios, for BOTH bias directions.**

| Z-bias ratio (λ/γ) | 3q-bitflip avg | 3q-phaseflip avg | [[5,1,3]] avg | Winner |
|---------------------|---------------|-----------------|--------------|--------|
| 0 | 0.949 | 0.868 | 0.990 | [[5,1,3]] |
| 2 | 0.717 | 0.804 | 0.937 | [[5,1,3]] |
| 5 | 0.537 | 0.754 | 0.853 | [[5,1,3]] |
| 10 | 0.390 | 0.643 | 0.702 | [[5,1,3]] |

Key findings:
- **Phase-flip code beats bit-flip under Z-bias** (as expected) but never catches [[5,1,3]]
- **Specialized code asymmetry is catastrophic**: bit-flip under Z-bias has asymmetry 0.9 (nearly useless for X/Y info)
- **[[5,1,3]] asymmetry stays below 0.06** even at extreme bias
- The **gap narrows** at high noise but [[5,1,3]] never loses
- Under X-biased noise, the same pattern holds in reverse (bit-flip > phase-flip, but both lose to [[5,1,3]])

**Surprise:** The phase-flip code (Z-specialized) does NOT dominate under Z-bias. At λ/γ=10, it achieves avg=0.643 vs [[5,1,3]]'s 0.702. Specialization cannot compensate for leaving 2/3 of the Bloch sphere unprotected.

### 22b: Isotropy advantage map across (γ, λ) landscape

**Setup:** 10×10 grid over γ∈[0.01,0.40], λ∈[0.01,0.40]. Basis-averaged Holevo. Compared [[5,1,3]], both 3-qubit codes, uncoded. "Bias-aware oracle" picks best 3-qubit code at each point.

**Result: No specialized code ever wins ANY grid point.**

Winner map (5=[[5,1,3]], U=uncoded):
- [[5,1,3]] wins 86/100 grid points (all of physically relevant low-to-moderate noise)
- Uncoded wins 14/100 (high noise corner: γ,λ > 0.27)
- Specialized codes win 0/100

Key findings:
- Peak [[5,1,3]] advantage: **+0.178** at balanced moderate noise (γ=λ≈0.10-0.14)
- [[5,1,3]] vs uncoded break-even: diagonal around γ=λ≈0.33
- The real competition is **isotropic vs uncoded**, not isotropic vs specialized
- Bias-aware oracle (picking best 3-qubit code for each noise point) still never beats [[5,1,3]]
- At extreme noise, code overhead (5 qubits = 5 noise targets) becomes the liability

**Surprise:** Specialization is **never** the optimal strategy at any noise point. The landscape has exactly two regimes: (1) isotropy wins (low-moderate noise), (2) no code helps (high noise). There is no "specialization zone."

### 22c: Why isotropy always wins — distance decomposition and quality vs quantity

**Four analyses:**

**1. Per-error-type protection (Pauli channel):**
[[5,1,3]] has **perfect balance** (1.000) against X, Y, Z errors at all error rates. 3-qubit codes have balance ~0.5 — catastrophically vulnerable to their unprotected error type. Even uncoded has perfect balance (it's just bad everywhere equally).

**2. Per-qubit efficiency:**
| Noise (γ,λ) | Uncoded H/n | [[5,1,3]] H/n | Ratio |
|---|---|---|---|
| (0.05, 0.05) | 0.840 | 0.193 | 4.4x worse |
| (0.10, 0.10) | 0.732 | 0.176 | 4.2x worse |
| (0.20, 0.20) | 0.565 | 0.133 | 4.2x worse |

Per-physical-qubit efficiency: uncoded wins overwhelmingly (4x better).

**3. Hybrid strategy (5-qubit budget):**
- [[5,1,3]]: 0.963 Holevo, 1 logical qubit
- 3q-phaseflip + 1 uncoded: 1.665 Holevo, 2 logical qubits (4 physical)
- 5× uncoded: 4.202 Holevo, 5 logical qubits

For **total information throughput**, uncoded wins 4.4x. QEC sacrifices quantity for quality.

**4. Distance decomposition:**
| Code | d_X | d_Y | d_Z | Code distance |
|---|---|---|---|---|
| 3q-bitflip | 3 | 1 | 1 | **1** |
| 3q-phaseflip | 1 | 1 | 3 | **1** |
| [[5,1,3]] | 3 | 3 | 3 | **3** |

Code distance = min(d_X, d_Y, d_Z). Specialized codes are effectively distance-1 overall.

**Key insight:** Isotropy wins because the weakest error direction determines the code distance, and unknown quantum states have components in ALL directions. But there's a deeper lesson: QEC's overhead means it's a quality-over-quantity trade. 5 uncoded qubits carry 4.4x more total information than 1 [[5,1,3]]-encoded logical qubit. **QEC only makes sense when you need high-fidelity individual logical qubits** (i.e., for computation where gate errors compound), not for maximum information throughput.

## Summary

Sprint 022 definitively answers: **bias-awareness never beats isotropy at small scale.** Three key results:

1. **No crossover exists.** Even at 10:1 noise bias, [[5,1,3]] beats every specialized code on basis-averaged Holevo. The noise landscape has exactly two regimes: isotropy wins (86%) or uncoded wins (14%). There is no "specialization zone."

2. **Distance decomposition explains why.** Specialized codes are effectively distance-1 against their unprotected errors. The code distance is the minimum over all error directions, and unknown quantum states probe all directions.

3. **Quality vs quantity tradeoff.** Per-physical-qubit, uncoded is 4x more efficient. QEC trades quantity for quality — it only makes sense when individual qubit fidelity matters more than throughput (i.e., for computation, not communication).

**Connection to literature:** The XZZX surface code's success at large scale relies on having enough qubits to maintain isotropy against common errors while gaining extra distance against the biased error. At our small scale (3-5 qubits), there aren't enough qubits to do both — you can have isotropy OR specialization, not both. The crossover to bias-tailoring advantage likely requires O(d²) qubits where d is the target distance.

## Next Sprint Ideas

- **Concatenated bias-tailoring:** Can concatenating a phase-flip code inside a bit-flip code (like Shor's 9-qubit) achieve effective isotropy from specialized components?
- **Quantum channel capacity under bias:** Is the capacity ordering different under biased noise? Does the capacity gap widen or narrow?
- **Real hardware comparison:** The noise fingerprint framework from Sprint 016 can now be tested — predict which code wins on real IBM hardware given its measured T1/T2 ratio
- **Entanglement-assisted codes:** Sprint 019 showed 2-4x capacity boost from pre-shared entanglement. How does this interact with bias?
- **Surface code at distance 2:** Can we fit a planar d=2 surface code in our 10-qubit budget and test bias-tailoring?
