# Sprint 021 — Combined T1+T2 Noise: The Realistic Error Landscape for QEC Codes

**Date:** 2026-03-31
**Status:** In Progress

## Motivation

Every previous noise analysis (Sprints 009, 016, 018, 019) examined amplitude damping and phase damping *separately*. Real superconducting hardware has both simultaneously, with T2 ≤ 2T1 (Lindblad bound). The T2/T1 ratio is the key hardware parameter — dephasing-dominated (T2 << T1) vs relaxation-dominated (T2 ≈ 2T1) regimes may favor different codes.

**Literature search:** No paper systematically maps QEC code performance across the T2/T1 ratio for small stabilizer codes ([[5,1,3]], Steane, 3-qubit, toric). The gap is clean and novel. Key references:
- arXiv:2512.09189 — Pauli twirling misdescribes thermal noise by 2-10x for surface codes
- arXiv:2603.04564 — First break-even under combined noise, but dephasing suppressed by DD
- arXiv:2402.16937 — Exact toric code coherent info under Pauli decoherence (not combined channel)

**Prediction:** Code performance ordering will depend on T2/T1 ratio, with crossover points where one code overtakes another. [[5,1,3]] should do relatively better when T2 << T1 (dephasing dominated) since it corrects phase errors. 3-qubit code should do best when T2 ≈ 2T1 (relaxation dominated).

## Experiments

### 21a: Combined T1+T2 Fidelity Landscape
**Goal:** Sweep (γ, λ) parameter space for 3-qubit, [[5,1,3]], and toric codes. Map fidelity/Holevo across the 2D noise plane.

### 21b: Code Selection Phase Diagram
**Goal:** At each point in (total_noise, T2/T1_ratio) space, determine which code has best logical fidelity. Find crossover boundaries.

### 21c: Holevo Information Under Combined Noise
**Goal:** Track logical information preservation as a function of T2/T1 ratio at fixed total noise. Compare with channel capacity.

---

## Results

*(To be filled as experiments complete)*
