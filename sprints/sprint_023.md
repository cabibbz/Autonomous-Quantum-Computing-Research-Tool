# Sprint 023 — Concatenated Bias-Tailoring: Does Two-Level Structure Unlock Specialization?

**Date:** 2026-03-31
**Status:** In Progress

## Motivation

Sprint 022 showed that single-level specialized codes (3-qubit bit-flip, 3-qubit phase-flip) NEVER beat the isotropic [[5,1,3]] code under any noise bias ratio. But the Shor [[9,1,3]] code has a unique two-level concatenated structure: inner phase-flip protection (within each 3-qubit block) and outer bit-flip protection (across blocks). This mirrors what large-scale bias-tailored codes (XZZX surface codes) do — use concatenation to get BOTH isotropy and specialization.

**Key question:** Does the Shor code's concatenated structure create a specialization advantage that single-level codes can't achieve?

## Literature Search

- "Elevator Codes" (arXiv:2601.10786): Concatenated bias-tailored codes for biased noise. Concatenated repetition codes outperform thin surface codes at bias >7×10⁴.
- "Hardware-efficient QEC via concatenated bosonic qubits" (Nature 2025): Cat qubits with native Z-bias, corrected by outer repetition code. Hardware-level bias preservation.
- Key gap: No systematic basis-averaged Holevo comparison of Shor vs [[5,1,3]] under combined T1+T2 at small scale.

## Experiments

### 23a: Shor [[9,1,3]] vs All Codes Across Noise Landscape

[Results pending]

### 23b: Decomposing Shor's Two-Level Protection

[Results pending]

### 23c: Concatenation Depth vs Bias — Predicting the Crossover

[Results pending]
