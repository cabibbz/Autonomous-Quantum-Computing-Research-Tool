# Sprint 023 — Concatenated Bias-Tailoring: Does Two-Level Structure Unlock Specialization?

**Date:** 2026-03-31
**Status:** Complete (3/3 experiments)

## Motivation

Sprint 022 showed that single-level specialized codes (3-qubit bit-flip, 3-qubit phase-flip) NEVER beat the isotropic [[5,1,3]] code under any noise bias ratio. But the Shor [[9,1,3]] code has a unique two-level concatenated structure: inner phase-flip protection (within each 3-qubit block) and outer bit-flip protection (across blocks). This mirrors what large-scale bias-tailored codes (XZZX surface codes) do — use concatenation to get BOTH isotropy and specialization.

**Key question:** Does the Shor code's concatenated structure create a specialization advantage that single-level codes can't achieve?

## Literature Search

- "Elevator Codes" (arXiv:2601.10786): Concatenated bias-tailored codes for biased noise. Concatenated repetition codes outperform thin surface codes at bias >7×10⁴.
- "Hardware-efficient QEC via concatenated bosonic qubits" (Nature 2025): Cat qubits with native Z-bias, corrected by outer repetition code. Hardware-level bias preservation.
- Key gap: No systematic basis-averaged Holevo comparison of Shor vs [[5,1,3]] under combined T1+T2 at small scale.

## Experiments

### 23a: Shor [[9,1,3]] vs All Codes Across Noise Landscape

6×6 (γ,λ) grid, basis-averaged Holevo. Compared uncoded, 3q-bitflip, 3q-phaseflip, [[5,1,3]], Shor [[9,1,3]].

**Results:**
- Shor wins **67%** of noise landscape on basis-averaged Holevo!
- [[5,1,3]] wins 28%, uncoded wins 6% (high noise corner)
- Single-level specialized codes (3q-bf, 3q-pf) win **zero** points
- Shor advantage largest at moderate-high γ (amplitude damping dominant)
- [[5,1,3]] advantage at high λ (dephasing dominant)
- Shor much less isotropic: avg asymmetry 0.34 vs [[5,1,3]]'s 0.02

**Surprise:** Shor wins on average Holevo despite Sprint 022 showing isotropy always wins! The key: Shor has 9 qubits (more redundancy) AND two-level structure.

### 23b: Decomposing Shor's Two-Level Protection

Compared Shor vs 9-qubit repetition codes (single-level, same qubit count) to separate structure from redundancy.

**Results:**
- **Structure ALWAYS beats redundancy** — Shor beats best 9-qubit repetition code at ALL 6 test points
- Advantage ranges from +0.09 (dephasing-dominated) to +0.44 (amplitude-damping dominated)
- Shor's X-basis Holevo is nearly perfect (0.998) even under heavy dephasing — inner phase-flip code specifically protects X
- Z and Y bases collapse under dephasing (Z/X ratio drops to 0.41 at high noise)
- **Shor is self-dual** under inner/outer concatenation swap: phase-flip-inside-bit-flip = bit-flip-inside-phase-flip for 3-qubit codes

**Key insight:** The concatenated structure creates EXTREME anisotropy — near-perfect protection for one Bloch sphere direction, poor for others. This is the structural signature of bias-tailored codes.

### 23c: Isotropy-Adjusted Comparison — The Figure of Merit Matters

Tested multiple figures of merit: average Holevo, min-basis Holevo (worst case), per-qubit efficiency.

**Results:**

| Metric | [[5,1,3]] wins | Shor wins | Uncoded wins |
|--------|----------------|-----------|--------------|
| Avg Holevo | 28% | 67% | 6% |
| **Min-basis Holevo** | **81%** | **8%** | **11%** |
| Per-qubit avg | 0% | 0% | 100% |
| Per-qubit min | 0% | 0% | 100% |

- Min-basis Holevo (worst-case basis = actual quantum computation ability): [[5,1,3]] dominates 81%
- Shor only wins 3/36 min-basis points (narrow strip at γ≈0.12-0.26, λ=0.05)
- Per-qubit efficiency: uncoded ALWAYS wins — QEC never beats uncoded for information density
- Head-to-head at balanced noise: Shor avg=0.569 > [[5,1,3]] avg=0.528, but Shor min=0.415 < [[5,1,3]] min=0.524

## Key Insights

1. **The figure of merit IS the answer.** "Which code is best?" is meaningless without specifying the metric. Average-case vs worst-case completely reverses the ranking.

2. **For quantum computation, min-basis Holevo is the correct metric** — you need to protect ALL Bloch sphere directions for universal quantum computing. [[5,1,3]]'s isotropy makes it the clear winner (81% of landscape).

3. **For classical-like storage of a known basis**, Shor's extra redundancy helps (67% on average Holevo). But this is not quantum error correction's purpose.

4. **Concatenation genuinely creates structural advantage** — Shor always beats 9-qubit repetition codes. The two-level structure handles different error types at different levels. But at 9 qubits, it can't achieve both specialization AND isotropy.

5. **The Shor code is self-dual** under inner/outer swap — the "phase-flip inside bit-flip" vs "bit-flip inside phase-flip" distinction doesn't exist for 3-qubit inner codes.

6. **QEC codes never beat uncoded per qubit** — error correction trades information density for information quality. It's investment, not amplification.

7. **This resolves the Sprint 022 story:** Isotropy always wins *for quantum computation* (min-basis). Shor's average-Holevo advantage is real but illusory for QC — it comes from extreme X-axis specialization that sacrifices Z/Y protection.

## Connection to Large-Scale Bias-Tailoring

The literature shows bias-tailored codes (XZZX surface code, elevator codes) succeed at large scale because they have enough qubits for BOTH high distance against biased errors AND adequate distance against all errors. At 9 qubits (Shor), you get high distance in one direction but d_min = 3 (same as [[5,1,3]]). The crossover to bias-tailoring advantage requires O(d²) qubits — consistent with our finding that Shor (9 qubits, d=3) can't beat [[5,1,3]] (5 qubits, d=3) on min-basis metrics.

## Sources

- [Elevator Codes (arXiv:2601.10786)](https://arxiv.org/html/2601.10786)
- [Hardware-efficient QEC via concatenated bosonic qubits (Nature 2025)](https://www.nature.com/articles/s41586-025-08642-7)
- [Bias-tailored quantum LDPC codes (arXiv:2202.01702)](https://arxiv.org/abs/2202.01702)
- [Tailored XZZX codes for biased noise](https://arxiv.org/pdf/1812.08186)
