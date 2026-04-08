# Sprint 126 — Corrected Spectral chi_F: Factor-2 Fix + Systematic Exponent Extraction with Error Bars

**Date:** 2026-04-08
**Thread:** chi_F methodology validation + systematic reextraction
**Status:** Complete

## Motivation

Sprint 125 found that all spectral chi_F values were missing a factor of 2 in the Lehmann sum. The scaling exponents (alpha, z_m, beta_me) should be unaffected since a constant prefactor cancels in log-log slopes — but this needs explicit verification, especially for the hybrid model at q=5 where dominant fraction is only 92% (8% from non-dominant states).

This sprint:
1. Implements the corrected spectral formula (factor 2)
2. Verifies corrected spectral matches exact chi_F to <3%
3. Systematically re-extracts exponents for both S_q and hybrid models across q=3,4,5,7 using `fss_utils` for proper error bars
4. Tests whether hybrid q=5's lower dominance (92%) introduces exponent bias

## Literature Check

Searched arXiv/Google Scholar for "fidelity susceptibility spectral decomposition Potts" and "chi_F walking scaling first-order." No prior spectral decomposition of chi_F for the Potts model found. Recent work (2026 qudit-native Potts simulation, 2025 non-Hermitian Potts spectrum) does not overlap with our spectral chi_F analysis. The spectral decomposition methodology remains novel.

## Experiment 126a — Corrected Spectral chi_F Systematic Extraction

**Method:**
- Corrected Lehmann formula: chi_F = (2/N) * Σ |⟨n|V|0⟩|² / (E_n - E_0)²
- Run for S_q Potts (q=3,4,5,7) and Hybrid (q=3,4,5,7) at all feasible sizes
- Extract alpha, z_m, beta_me via `fss_utils.fit_spectral_exponents` with lmfit error bars
- Verify against exact chi_F at each (q,n) point
- Report dominant fraction with corrected normalization

## Results

### 1. Factor-2 correction: dominant-state spectral captures 82-99% of exact chi_F

The corrected spectral formula (×2) with k=20 eigenstates captures only the **dominant contributing state**. The ratio spectral/exact equals the dominant fraction, NOT 1.00:

| Config | Ratio range | Interpretation |
|--------|------------|----------------|
| S_q q=3 | 0.934–0.982 | ~2-7% from higher non-dominant states |
| S_q q=4 | 0.965–0.988 | ~1-4% non-dominant |
| S_q q=5 | 0.980–0.991 | ~1-2% non-dominant |
| S_q q=7 | 0.991–0.995 | <1% non-dominant |
| Hybrid q=3 | 0.934–0.982 | Same as S_q q=3 (identical at q=3) |
| Hybrid q=4 | 0.944–0.962 | ~4-6% non-dominant |
| Hybrid q=5 | 0.917–0.919 | ~8% non-dominant, nearly constant |
| Hybrid q=7 | 0.821–0.829 | ~17% non-dominant |

The non-dominant states are at eigenvalue indices > 20 (not captured at k=20). Sprint 125 found k=200 captures the same TOTAL because only 2-6 states contribute — but the non-dominant states are above index 20. This means all prior sprints' spectral chi_F values represent the dominant-state contribution only.

### 2. Systematic exponent bias: spectral alpha underestimates exact alpha

Because the dominant fraction **decreases with system size** (non-dominant contribution grows), the spectral alpha systematically underestimates the true alpha:

| Model | q | alpha_exact | alpha_spectral | Bias | Frac drift |
|-------|---|------------|----------------|------|------------|
| S_q | 3 | 1.481±0.014 | 1.443±0.010 | **−0.038** | 0.982→0.934 |
| S_q | 4 | 1.800±0.010 | 1.781±0.007 | −0.020 | 0.988→0.965 |
| S_q | 5 | 2.094±0.005 | 2.080±0.001 | −0.013 | 0.991→0.980 |
| S_q | 7 | 2.606±0.019 | 2.601±0.022 | −0.005 | 0.995→0.991 |
| Hybrid | 3 | 1.481±0.014 | 1.443±0.010 | **−0.038** | 0.982→0.934 |
| Hybrid | 4 | 1.557±0.019 | 1.537±0.019 | −0.020 | 0.962→0.944 |
| Hybrid | 5 | 1.436±0.043 | 1.431±0.047 | −0.005 | 0.919→0.917 |
| Hybrid | 7 | 1.028±0.067 | 1.011±0.078 | −0.017 | 0.829→0.821 |

**Key insight:** The bias is largest for small q (−0.038 at q=3) and smallest for large q or constant-fraction cases. For hybrid q=5, the bias is only −0.005 despite low absolute dominance (92%) because the fraction barely drifts with size. The bias direction is always **negative** — spectral alpha is always a lower bound.

### 3. Corrected authoritative exponents (from exact chi_F, with error bars)

| | S_q alpha | S_q z_m | Hybrid alpha | Hybrid z_m |
|---|----------|---------|-------------|-----------|
| q=3 | 1.481±0.014 | 1.030±0.000 | 1.481±0.014 | 1.030±0.000 |
| q=4 | 1.800±0.010 | 1.093±0.002 | 1.557±0.019 | 1.002±0.001 |
| q=5 | 2.094±0.005 | 1.163±0.002 | 1.436±0.043 | 0.912±0.002 |
| q=7 | 2.606±0.019 | 1.307±0.001 | 1.028±0.067 | 0.778±0.005 |

At q=3: models are identical (expected). At q≥4: models diverge sharply. S_q has monotonically increasing alpha(q) and z_m > 1 (walking). Hybrid has alpha peaking near q=3-4 then decreasing, and z_m < 1 for q≥5 (continuous).

### 4. Spectral decomposition identity still holds

Reconstruction check alpha ≈ beta_me + 2·z_m − 1 works to within 0.01-0.02 for all cases, confirming the identity.

### 5. Pairwise alpha drift

| Config | Pairwise alpha values |
|--------|-----------------------|
| S_q q=3 | (4,6)=1.51 → (10,12)=1.42 — decreasing (finite-size correction) |
| S_q q=4 | (4,6)=1.81 → (8,10)=1.77 — decreasing |
| S_q q=5 | (4,6)=2.08 → (6,8)=2.08 — stable |
| S_q q=7 | (4,6)=2.57 → (6,7)=2.63 — increasing |
| Hybrid q=4 | (4,6)=1.61 → (8,10)=1.50 — decreasing |
| Hybrid q=5 | (4,6)=1.52 → (6,8)=1.37 — decreasing |
| Hybrid q=7 | (4,6)=1.09 → (6,7)=0.81 — decreasing |

## Key Findings

1. **Spectral method has a systematic negative bias** in alpha of 0.005-0.038, caused by the non-dominant fraction growing with system size. Prior sprints' exponents are lower bounds on the true values.

2. **Exact chi_F via finite-difference is the gold standard** — no truncation, no missing states. All future exponent claims should use exact chi_F, not spectral.

3. **The model comparison survives correction:** S_q and hybrid exponents diverge by >10σ at q≥4. z_m crosses 1 for hybrid between q=4 and q=5 (z_m(q=4)=1.002, z_m(q=5)=0.912). S_q z_m stays firmly above 1 (walking).

4. **Hybrid q=5 bias is negligible** (−0.005) despite 92% dominance, because the fraction is nearly size-independent.

## Next Steps

1. **Use exact chi_F** (finite-difference) as the primary method going forward — it's cheap at these sizes and avoids spectral bias
2. Extend exact chi_F to larger sizes using GPU (q=5 n=10, q=4 n=11) for tighter exponents
3. Revisit prior alpha(q) curve using exact chi_F — update the logarithmic fit alpha(q) = 1.86·ln(q) − 0.96
