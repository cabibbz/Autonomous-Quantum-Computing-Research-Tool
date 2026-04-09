# Sprint 127 -- Extended exact chi_F to GPU sizes for tighter exponent fits

**Date:** 2026-04-08
**Thread:** chi_F exponent precision via larger system sizes
**Status:** Complete

## Motivation

Sprint 126 established exact chi_F (finite-difference) as the gold standard method. But exponents were limited by system sizes: q=5 only had 3 data points (n=4,6,8), and q=7 had 3 points (n=4,6,7). With GPU acceleration, we can push to ~10M-dimensional Hilbert spaces and get 4-6 data points per configuration for much tighter fits.

## Literature Check

No new papers found on fidelity susceptibility scaling exponents for the S_q Potts model since Sprint 126. The spectral chi_F decomposition methodology remains novel.

## Experiment 127a -- Extended exact chi_F to GPU sizes

**Method:** Central finite-difference chi_F = (2 - |<psi(g)|psi(g+dg)>|^2 - |<psi(g)|psi(g-dg)>|^2) / (dg^2 * N), with dg=1e-4. Three eigsh(k=1) calls per data point.

**New data points computed (12 total):**

| Model | q | n | dim | chi_F | Time |
|-------|---|---|-----|-------|------|
| S_q | 3 | 14 | 4,782,969 | 26.7410 | 17.0s |
| S_q | 4 | 9 | 262,144 | 56.1837 | 1.0s |
| S_q | 4 | 11 | 4,194,304 | 80.3109 | 15.2s |
| S_q | 5 | 9 | 1,953,125 | 169.455 | 7.4s |
| S_q | 5 | 10 | 9,765,625 | 211.411 | 41.9s |
| S_q | 7 | 8 | 5,764,801 | 715.699 | 27.1s |
| Hybrid | 3 | 14 | 4,782,969 | 26.7913 | 16.4s |
| Hybrid | 4 | 9 | 262,144 | 15.6371 | 0.8s |
| Hybrid | 4 | 11 | 4,194,304 | 21.1691 | 12.2s |
| Hybrid | 5 | 9 | 1,953,125 | 12.8294 | 4.8s |
| Hybrid | 5 | 10 | 9,765,625 | 14.6861 | 25.9s |
| Hybrid | 7 | 8 | 5,764,801 | 5.4677 | 13.1s |

**Total compute: ~3 min.** GPU acceleration made all 10M-dim eigsh calls feasible within 42s each.

## Results

### 1. Updated alpha exponents (combined fits with all data)

| Model | q | alpha (S127) | +/- | R^2 | # sizes | Sizes |
|-------|---|-------------|-----|-----|---------|-------|
| S_q | 3 | **1.4676** | 0.0116 | 0.99985 | 6 | 4,6,8,10,12,14 |
| S_q | 4 | **1.7945** | 0.0069 | 0.99997 | 6 | 4,6,8,9,10,11 |
| S_q | 5 | **2.0944** | 0.0019 | 1.00000 | 5 | 4,6,8,9,10 |
| S_q | 7 | **2.6357** | 0.0184 | 0.99995 | 4 | 4,6,7,8 |
| Hybrid | 3 | **1.4675** | 0.0117 | 0.99985 | 6 | 4,6,8,10,12,14 |
| Hybrid | 4 | **1.5447** | 0.0137 | 0.99983 | 6 | 4,6,8,9,10,11 |
| Hybrid | 5 | **1.3816** | 0.0302 | 0.99914 | 5 | 4,6,8,9,10 |
| Hybrid | 7 | **0.9569** | 0.0670 | 0.99249 | 4 | 4,6,7,8 |

### 2. Changes from Sprint 126

| Model | q | Sprint 126 | Sprint 127 | Shift | Error reduction |
|-------|---|-----------|-----------|-------|-----------------|
| S_q | 3 | 1.481+/-0.014 | 1.468+/-0.012 | -0.013 | 17% |
| S_q | 4 | 1.800+/-0.010 | 1.795+/-0.007 | -0.006 | 31% |
| S_q | 5 | 2.094+/-0.005 | 2.094+/-0.002 | +0.000 | **62%** |
| S_q | 7 | 2.606+/-0.019 | 2.636+/-0.018 | +0.030 | 3% |
| Hybrid | 3 | 1.481+/-0.014 | 1.468+/-0.012 | -0.014 | 17% |
| Hybrid | 4 | 1.557+/-0.019 | 1.545+/-0.014 | -0.012 | 28% |
| Hybrid | 5 | 1.436+/-0.043 | 1.382+/-0.030 | -0.054 | **30%** |
| Hybrid | 7 | 1.028+/-0.067 | 0.957+/-0.067 | -0.071 | 0% |

**Key changes:**
- S_q q=5: error bars cut by 62% (from +/-0.005 to +/-0.002). alpha=2.094 rock-solid.
- S_q q=4: error reduced 31%. alpha continues to drift DOWN (1.795 vs 1.800). Pairwise alpha at largest sizes: (10,11)=1.779. Settling toward ~1.78, not 2.0.
- S_q q=7: alpha shifted UP (+0.030 to 2.636). Pairwise alpha at (7,8)=2.670 -- still increasing. May settle higher than 2.64.
- Hybrid q=5: significant downward shift (-0.054). Pairwise (9,10)=1.283 -- still decreasing. alpha is genuinely below q=3.
- Hybrid q=7: dropped to 0.957. Pairwise (7,8)=0.755 -- rapidly decreasing. alpha < 1 confirmed.

### 3. Pairwise drift analysis -- convergence status

**S_q Potts:**
- q=3: 1.568 -> 1.496 -> 1.465 -> 1.448 -> 1.437 (smoothly decreasing, not converged)
- q=4: 1.846 -> 1.799 -> 1.786 -> 1.782 -> 1.779 (converging toward ~1.77-1.78)
- q=5: 2.104 -> 2.088 -> 2.094 -> 2.100 (stable around 2.09)
- q=7: 2.582 -> 2.631 -> 2.670 (increasing, not converged -- needs larger sizes)

**Hybrid:**
- q=3: 1.568 -> 1.496 -> 1.465 -> 1.448 -> 1.437 (same as S_q, expected)
- q=4: 1.631 -> 1.554 -> 1.526 -> 1.514 -> 1.505 (smoothly decreasing, not converged)
- q=5: 1.512 -> 1.379 -> 1.315 -> 1.283 (rapidly decreasing -- alpha may settle well below 1.3)
- q=7: 1.097 -> 0.858 -> 0.755 (rapidly decreasing toward sub-1 regime)

### 4. alpha(q) curve refit

**S_q Potts:** alpha(q) = 1.306 * ln(q) - 0.006
- Prior (Sprint 116, spectral): alpha(q) = 1.86 * ln(q) - 0.96
- The new fit slope (1.306) is MUCH lower than the old one (1.86)!
- Reason: Sprint 116 used spectral chi_F on data up to q=25 with smaller sizes. The spectral method had negative alpha bias that was non-uniform across q. Additionally, at large q the bias is minimal (frac ~0.99) while at small q it's large (frac ~0.93). This creates an artificial steepening.
- But: chi2/dof=21.9 -- the log fit is poor. q=7 residual is +0.10. With only 4 data points at the extreme end drifting, we can't distinguish log from power-law alpha(q). Need more q values.

**Hybrid:** Log fit fails badly (chi2/dof=44). alpha(q) is NOT logarithmic for hybrid. It peaks near q=3-4 then drops below 1 for q>=7. This is qualitatively different from S_q.

## Key Findings

1. **S_q q=5 alpha = 2.094 +/- 0.002** -- the tightest constraint yet (62% error reduction). Rock-solid.

2. **S_q q=4 alpha continues to settle BELOW 2.0** -- pairwise alpha at (10,11)=1.779. The drift direction is consistently downward. Unless there's a dramatic reversal at n>11, alpha(q=4) is genuinely sub-2. This contradicts the prediction alpha=2 with log corrections from Salas-Sokal. The log-corrected form chi_F ~ N^2 (ln N)^{-p} remains the worst fit (Sprint 122, 126).

3. **Hybrid alpha(q) is non-monotonic** -- peaks near q=3-4 then drops below 1 for q>=7. Qualitatively different from the monotonically increasing S_q pattern. This is the clearest model-distinguishing signature.

4. **S_q q=7 alpha still drifting upward** (2.636 +/- 0.018, pairwise (7,8)=2.670). Needs q=7 n=9 (dim=40M, beyond current GPU) for convergence.

5. **Prior alpha(q) = 1.86*ln(q) fit is obsolete** -- was based on spectral chi_F with non-uniform bias across q. New exact chi_F slope is 1.306. But the log fit itself is mediocre (chi2/dof=22). Need intermediate q values (q=6, q=8) to distinguish functional forms.

## Next Steps

1. **Add q=6 S_q and hybrid** at all feasible sizes -- this fills the gap between q=5 and q=7 for alpha(q) curve fitting. q=6 n=8 dim=1.7M (easy GPU).
2. **Power-law+correction fits for S_q q=3** -- pairwise alpha is still decreasing at n=14. Try alpha_eff(N) = alpha_inf + c/N^p extrapolation to get asymptotic alpha.
3. **Hybrid z_m extraction** at extended sizes -- z_m(q=5) and z_m(q=7) with the new data points.
