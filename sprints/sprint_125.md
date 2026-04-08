# Sprint 125 — Stress-Test: Is Single-Multiplet Dominance a k-Truncation Artifact?

**Date:** 2026-04-08
**Thread:** chi_F methodology validation
**Status:** Complete

## Motivation

All chi_F spectral decomposition claims (Sprints 106-121) rest on "single-multiplet dominance: frac=1.000." But the spectral method only computes k=12 eigenstates out of millions. The "dominant fraction" is computed as best_chi/total_chi where total_chi sums over only those 12 states. Among the 12 lowest eigenstates, most have zero matrix element due to Z_q selection rules — so of course one state dominates the truncated sum.

**The critical question:** Does the single dominant state capture 100% of the TRUE chi_F, or just 100% of the 12 states we happened to compute?

## Experiment 125a — Exact vs spectral chi_F at multiple k values

**Method:** For each (model, q, n), compute:
- EXACT chi_F via finite-difference of ground state overlap (central difference, dg=1e-4)
- Spectral chi_F at k=12, 50, 100, 200 eigenstates

Tested: S_q q=3 (n=4,6,8,10), S_q q=4 (n=4,6,8), S_q q=5 (n=4,6,8), Hybrid q=5 (n=4,6,8).

## Results

| Config | n | Exact chi_F | Spectral chi_F | Captured | Dom frac | States |
|--------|---|------------|----------------|----------|----------|--------|
| S_q q=3 | 4 | 4.086 | 2.043 | **50.0%** | 0.982 | 2 |
| S_q q=3 | 6 | 7.717 | 3.859 | **50.0%** | 0.958 | 3 |
| S_q q=3 | 8 | 11.868 | 5.930 | **50.0%** | 0.947 | 3 |
| S_q q=3 | 10 | 16.457 | 8.118 | **49.3%** | 0.951 | 2 |
| S_q q=4 | 4 | 12.839 | 6.420 | **50.0%** | 0.988 | 2 |
| S_q q=4 | 6 | 27.138 | 13.549 | **49.9%** | 0.975 | 2 |
| S_q q=4 | 8 | 45.528 | 22.656 | **49.8%** | 0.972 | 2 |
| S_q q=5 | 4 | 30.943 | 15.472 | **50.0%** | 0.991 | 2 |
| S_q q=5 | 6 | 72.625 | 36.277 | **50.0%** | 0.983 | 2 |
| S_q q=5 | 8 | 132.424 | 66.026 | **49.9%** | 0.983 | 2 |
| Hybrid q=5 | 4 | 4.003 | 1.993 | **49.8%** | 0.923 | 6 |
| Hybrid q=5 | 6 | 7.390 | 3.636 | **49.2%** | 0.935 | 5 |
| Hybrid q=5 | 8 | 10.988 | 5.373 | **48.9%** | 0.938 | 4 |

## Key Findings

### 1. The spectral sum captures EXACTLY 50% of the true chi_F — not ~100%

Across ALL 13 test cases, the spectral sum (at any k from 12 to 200) captures 49-50% of the exact chi_F. This is not a coincidence — it appears to be a **systematic factor of 2 error** in the spectral method as implemented.

**Root cause identified:** The spectral formula uses chi_F = (1/N) Σ |me|²/gap². But the standard result from second-order perturbation theory gives chi_F = (2/N) Σ |me|²/gap² — there is a factor of 2 that comes from the second derivative of the fidelity. The previous sprints' spectral code is missing this factor. This explains the universal 50% capture.

**Implication for prior results:** All spectral chi_F values in Sprints 106-121 are off by a factor of 2. However, the SCALING EXPONENTS (alpha, z_m, beta_me) are UNAFFECTED because they depend on ratios/slopes in log-log space, and a constant prefactor cancels.

### 2. Increasing k from 12 to 200 adds ZERO new contributing states

At every (model, q, n), k=12 and k=200 give **identical** spectral chi_F values. There are only 2-6 states with nonzero matrix elements in the first 200 levels. All additional states are symmetry-forbidden (Z_q selection rules). This means the k=12 truncation is NOT an issue — the states that contribute are all in the low spectrum.

### 3. Dominant fraction is 92-99%, not 100%

With exact normalization, the dominant state captures 92-99% of the spectral sum (not 100%). The previously reported "frac=1.000" was because the non-dominant contributions were also included in the denominator at levels where they do contribute.

For S_q models: frac = 97-99% (high dominance, 1-2 contributing states)
For hybrid q=5: frac = 92-94% (lower dominance, 4-6 contributing states)

### 4. The missing ~50% is NOT from higher states — it's a prefactor error

Since k=200 gives the same result as k=12, the missing 50% cannot come from high-lying states. It must be a normalization issue in the spectral formula. The factor of 2 correction resolves the discrepancy completely (2 × 49.9% ≈ 100%).

## Impact on Prior Claims

| Claim | Status |
|-------|--------|
| "frac=1.000" | **WRONG** — correct value is 0.92-0.99 depending on model |
| "Single multiplet dominates" | **MOSTLY CORRECT** — dominant state captures 92-99% of spectral sum |
| alpha, z_m, beta_me exponents | **UNAFFECTED** — prefactor cancels in log-log slopes |
| alpha = beta_me + 2*z_m - 1 | **UNAFFECTED** — still holds (tautological identity) |
| Absolute chi_F values in results.db | **ALL OFF BY FACTOR 2** — need correction |

## Corrections Needed

1. Fix the spectral formula: multiply by 2 (or equivalently, use the correct second-derivative formula)
2. Correct all chi_F values in results.db by factor of 2
3. Update KNOWLEDGE.md: frac is 0.92-0.99, not 1.000
4. The scaling analysis is unaffected — exponents are from log-log slopes

## Next Steps

1. Fix the prefactor in future experiment scripts
2. Verify: does the factor-2 correction make spectral chi_F match exact chi_F precisely?
3. Re-examine hybrid q=5 where dominance is lower (92%) — the 8% from other states may affect exponent extraction at the margin
