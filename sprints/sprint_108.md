# Sprint 108 — q=4 BKT Crossover: Spectral Decomposition Deep Dive

**Date:** 2026-04-02
**Status:** In progress

## Motivation

q=4 S_q Potts sits exactly at the BKT boundary — it's the marginal case between continuous (q≤4) and first-order (q>4) in the thermodynamic limit. Sprint 106c measured α=1.69 with sizes n=4-8, and the linear formula predicts α=1.729. The question: does q=4 have BKT-type log corrections that modify the spectral decomposition, or is it already fully walking-type?

Key questions:
1. Does α(q=4) converge to the linear prediction 1.729 with more sizes, or deviate?
2. Are z_m and β_me individually stable, or do they have log corrections (BKT signature)?
3. How does the multiplet structure at q=4 compare to q=3 (continuous) and q=5 (walking)?

## Literature context

q=4 S_q Potts: c=1 exactly (Baxter). BKT-type transition with multiplicative log corrections. Known log corrections to energy gap and correlation length. Fidelity susceptibility at BKT points is expected to scale as N^{2/ν} with log corrections (Schwandt et al. 2009). For q=4 Potts, ν=2/3 → α_BKT = 2·(3/2) - 2 = 1, but with strong log corrections.

## Experiments

### 108a — GPU-extended q=4 spectral decomposition (n=4-10)

**7 sizes, n=4-10.** GPU for n=9 (dim=262k, 23s) and n=10 (dim=1M, 96s). Total 125s.

Key findings:
- **Dominant fraction = 1.000 for ALL sizes.** Single (q-1)=3-fold multiplet at level 4 captures 100% of χ_F. Same as q=5 (walking), unlike q=2 (87-93%).
- **Pairwise α converges DOWNWARD:** 1.825 → 1.794 → 1.780 → 1.774 → 1.772 → 1.771
- **Global α = 1.788** (7 sizes), above Sprint 106c's 1.69 (5 sizes n=4-8)
- **z_m slowly decreasing:** 1.097 → 1.082 (pairwise). Consistent with log correction.
- **β_me non-monotonic:** 0.630 → 0.597 → 0.601 → 0.607 (pairwise). Stabilizing near 0.60.
- **Decomposition exact:** α_pred = β_me + 2z_m - 1 matches α_meas to all digits at every pair.

| n | dim | gap_m | |me|² | χ_F | α_pw |
|---|-----|-------|-------|------|------|
| 4 | 256 | 0.9515 | 22.96 | 6.34 | — |
| 5 | 1024 | 0.7449 | 26.43 | 9.53 | 1.825 |
| 6 | 4096 | 0.6100 | 29.50 | 13.21 | 1.794 |
| 7 | 16384 | 0.5155 | 32.33 | 17.38 | 1.780 |
| 8 | 65536 | 0.4457 | 35.02 | 22.03 | 1.774 |
| 9 | 262144 | 0.3923 | 37.59 | 27.14 | 1.772 |
| 10 | 1048576 | 0.3500 | 40.07 | 32.71 | 1.771 |

**Surprise:** α is converging ABOVE the linear formula (1.729), not below. Sprint 106c's 1.69 was too low because it only had 5 sizes (n=4-8) and the small-n points pull it down. The asymptotic α appears to be ~1.77, roughly 2.5% above the linear formula.

### 108b — BKT log correction test: q=3 vs q=4 vs q=5

**7 sizes each (q=3,4: n=4-10; q=5: n=4-9).** Fit pairwise exponents vs 1/ln(N) to detect log corrections.

Key findings:
- **q=5 (walking): α has ZERO log correction.** α_log = -0.003/ln(N), R²=0.001. Pairwise α flat: 2.089→2.090. The decomposition formula α = β_me + 2z_m - 1 is exact, and while z_m and β_me individually drift (~0.07/ln(N)), they COMPENSATE exactly.
- **q=4 (BKT): α has moderate log correction.** α_log = +0.243/ln(N), R²=0.92. Extrapolated α_∞ = 1.657. This means current α(q=4) ≈ 1.77 is INFLATED by log corrections — the asymptotic value is lower.
- **q=3 (continuous): α has strong log correction.** α_log = +0.460/ln(N), R²=0.98. Extrapolated α_∞ = 1.217. The continuous regime has the strongest finite-size inflation.

| q | α_measured | α_∞ (extrapolated) | α_log coeff | z_m_∞ | β_me_∞ |
|---|-----------|--------------------|-----------:|------:|-------:|
| 3 | 1.43 (last pw) | 1.217 | +0.460 | 1.033 | 0.151 |
| 4 | 1.77 (last pw) | 1.657 | +0.243 | 1.054 | 0.549 |
| 5 | 2.09 (last pw) | 2.084 | -0.003 | 1.120 | 0.844 |

**Three-regime structure confirmed:**
1. **Continuous (q=3):** Large log corrections inflate α. True α_∞ ≈ 1.22. z_m flat, β_me drives drift.
2. **BKT boundary (q=4):** Moderate log corrections. True α_∞ ≈ 1.66, below current measured 1.77. z_m has log correction too.
3. **Walking (q=5):** No log corrections. α is the true asymptotic value. z_m and β_me drift cancel.

**The walking regime is distinguished not just by larger α, but by the ABSENCE of log corrections.** This makes walking exponents more reliable at finite size than continuous or BKT ones.

### 108c — Multiplet structure comparison (q=3, 4, 5)

**Detailed low-energy spectrum for q=3 (n=6-12), q=4 (n=6-10), q=5 (n=6-8).**

Key findings:
- **Spectral gap degeneracy = (q-1) for ALL q.** q=3: 2-fold, q=4: 3-fold, q=5: 4-fold. All symmetry-forbidden (|me|²=0 exactly).
- **Dominant multiplet is SINGLET for all q.** The single state at level q (i.e., level 3 for q=3, level 4 for q=4, level 5 for q=5) carries 100% of χ_F. Reported as "degeneracy=1" in all cases.
- **Gap ratio (multiplet/spectral) decreases with q:** q=3: 6.2, q=4: 5.2, q=5: 4.6. The multiplet is getting CLOSER to the spectral gap as q increases.
- **Gap ratio slowly decreasing with N** at all q (converging from above). q=3: 6.24→6.17. q=4: 5.28→5.13. q=5: 4.64→4.51.
- **Second excited multiplet** has degeneracy q: q=3: 3-fold, q=4: 4-fold, q=5: 5-fold (at n=6) or 4-fold (at n=8, splitting). All have |me|²=0.

| q | n | gap_spec | gap_mult | ratio | spec_deg | mult_deg |
|---|---|----------|----------|-------|----------|----------|
| 3 | 6 | 0.1232 | 0.7684 | 6.24 | 2 | 1 |
| 3 | 12 | 0.0610 | 0.3764 | 6.17 | 2 | 1 |
| 4 | 6 | 0.1155 | 0.6100 | 5.28 | 3 | 1 |
| 4 | 10 | 0.0682 | 0.3500 | 5.13 | 3 | 1 |
| 5 | 6 | 0.1084 | 0.5029 | 4.64 | 4 | 1 |
| 5 | 8 | 0.0799 | 0.3605 | 4.51 | 4 | 1 |

**Interpretation:** The spectral gap has (q-1) degenerate states but is symmetry-forbidden from coupling to the field. The χ_F-dominant state is a SINGLET above this multiplet. As q increases, this singlet gets closer to the spectral gap, explaining why the field response (and hence χ_F) grows faster with N for larger q.

## Summary

Three key results from Sprint 108:

1. **q=4 α ≈ 1.77 at finite size, but α_∞ ≈ 1.66** after log correction extrapolation. The linear formula α(q) = 0.315q + 0.469 = 1.729 overshoots the true asymptotic value.

2. **Walking (q=5) has zero log correction to α** (α_log = -0.003, R²=0.001). Continuous (q=3) and BKT (q=4) both have large positive log corrections. This means walking α is trustworthy at finite size, while q≤4 values are inflated.

3. **Gap ratio (multiplet/spectral) decreases monotonically with q:** 6.2→5.2→4.6. The χ_F-dominant singlet approaches the symmetry-forbidden spectral gap as q increases, concentrating the field response.

**Revised α(q) formula accounting for log corrections:**
- q=3: α_∞ = 1.22 (measured 1.43 at n=10)
- q=4: α_∞ = 1.66 (measured 1.77 at n=10)
- q=5: α_∞ = 2.08 (measured 2.09 at n=9, no correction needed)
- The linear formula α = 0.315q + 0.469 is only valid for q≥5 (walking regime) where log corrections vanish.

**NOT upgrading to novel** — the log correction finding is interesting but expected for BKT (q=4) and continuous (q=3) universality classes. The new insight is that walking ELIMINATES log corrections, making it a cleaner measurement regime.
