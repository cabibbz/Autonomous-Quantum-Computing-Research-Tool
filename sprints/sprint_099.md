# Sprint 099 — Complex CFT: Detect Im(c) Oscillations in Casimir Residuals

**Date:** 2026-04-02
**Thread:** New direction — oscillatory corrections from complex CFT
**Status:** Complete (3 experiments)

## Motivation

Complex CFT (Gorbenko, Rychkov & Zan, JHEP 2018) predicts that q>4 Potts models have complex scaling dimensions, producing **spiral RG flow** and **oscillatory corrections** to finite-size scaling. For q=5, c ≈ 1.138 ± 0.021i. Sprint 098 confirmed that Casimir energy tracks Re(c) to <1%. The natural next question: **do the Casimir residuals show oscillations predicted by Im(c)?**

Literature context:
- Gorbenko et al. (arXiv:1807.11512, 1808.04380): walking = flow between complex fixed points → "drifting scaling dimensions"
- arXiv:2502.02001: complex entanglement entropy for non-Hermitian q=5 Potts
- arXiv:2507.14732: spiral RG flow in entanglement spectrum of non-Hermitian q=5 Potts

**Key prediction:** corrections ~cos(ω·ln(N) + φ)/N⁴ with ω = 2·arccosh(√q/2). For q=5: ω=0.962, full period needs N ratio ~684. For q=7: ω=1.567, period needs N ratio ~55.

## Experiment 099a — Dense Casimir Scan

Generated E₀(N) at ALL integer sizes:
- q=2: N=4-14 (11 points, control — real CFT)
- q=5: N=4-10 (7 points, walking)
- q=7: N=4-8 (5 points, broken walking)

**New data (odd sizes not previously measured):**

| q | N | E₀/N | gap×N |
|---|---|------|-------|
| 2 | 5 | -1.147213595500 | 0.791922 |
| 2 | 7 | -1.141994172491 | 0.788711 |
| 2 | 9 | -1.139863387016 | 0.787398 |
| 2 | 11 | -1.138788562121 | 0.786736 |
| 2 | 13 | -1.138171523889 | 0.786355 |
| 5 | 5 | -1.142428477185 | 0.659091 |
| 5 | 7 | -1.133550397616 | 0.644379 |
| 5 | 9 | -1.129947738339 | 0.635172 |

3-param fits (ε + A/N² + B/N⁴):
- q=2: R²=0.99999999974, residuals ~10⁻⁷
- q=5: R²=0.99999998, residuals ~10⁻⁶ (10× larger)
- q=7: R²=0.99999988, residuals ~3×10⁻⁶ (30× larger)

**Residual sign changes:** q=2: 3/10, q=5: 3/6, q=7: 3/4. Residuals grow with q.

## Experiment 099b — Oscillatory vs Monotonic Model Comparison

Compared models for E₀/N:
- M1: ε + A/N² + B/N⁴ (3-param, monotonic)
- M2: ε + A/N² + B/N⁴ + D/N⁶ (4-param, monotonic)
- M3: ε + A/N² + C·cos(ω·ln(N)+φ)/N⁴ (4-param, oscillatory, ω fixed from complex CFT)
- M4: ε + A/N² + B/N⁴ + C·cos(ω·ln(N)+φ)/N⁴ (5-param, mixed)

| q | M2 RSS | M3 RSS | M3/M2 | Best model | C_osc |
|---|--------|--------|-------|------------|-------|
| 2 | 3.4e-19 | 1.7e-15 | 5051 | M2 (monotonic dominates) | -0.078 |
| 5 | 2.5e-13 | 3.8e-13 | 1.54 | M2 > M3, but M4 best (BIC=-220) | -0.208 |
| 7 | 4.4e-13 | 7.0e-13 | 1.61 | M2 (monotonic wins) | -0.295 |

**Oscillatory model DOES NOT beat monotonic** at same parameter count. The 1/N⁶ correction is more important than oscillatory corrections.

For q=5, M4 (mixed 5-param) beats M2 (BIC=-220 vs -209), but with only 7 points and 5 params, overfitting risk is high.

**Pairwise c/Re(c) for q=5 is NON-MONOTONIC:**
(4,5): 1.005, (5,6): 1.002, (6,7): 1.003, (7,8): 1.004, (8,9): 1.006, (9,10): 1.008

First decreases then reverses upward. q=2 is perfectly monotonically decreasing (9/9 decreases).

## Experiment 099c — Separating Casimir from Velocity

**Critical finding:** the pairwise vc (Casimir product v×c) is MONOTONICALLY DECREASING for ALL q including q=5. The non-monotonicity in c/Re(c) comes from the velocity correction: v decreases faster than vc, causing c = vc/v to increase.

Detrended vc (after removing linear 1/N² trend):
- q=2: RMS = 4.75e-5, pattern `+-----++++` (smooth curvature)
- q=5: RMS = 3.97e-4, pattern `-+++--` (8.4× larger, one "wiggle")
- q=7: RMS = 9.96e-4, pattern `-++-` (21× larger)

The RMS scales with q, consistent with larger corrections for larger q. But the detrended sign changes don't match the predicted omega — the observed "frequency" is much higher than complex CFT predicts.

**Richardson extrapolation reveals persistent drift for q>4:**
- q=2: extrapolated vc converges cleanly to 0.4999 → c ≈ 0.500 ✓
- q=5: keeps drifting downward (0.841 → 0.837), not converging
- q=7: drifts faster (0.878 → 0.866)

This persistent drift could be from logarithmic corrections or very slowly varying oscillatory corrections.

## Conclusions

1. **Im(c) oscillations are UNDETECTABLE at accessible exact diag sizes.** The predicted period is 2π/ω = 6.53 in ln(N) for q=5, requiring N span from ~4 to ~2700 for one cycle. Our range (N=4-10) covers only 14% of a period.

2. **Monotonic corrections (1/N⁶) dominate over oscillatory corrections.** M2 beats M3 for all q at same parameter count.

3. **Pairwise c/Re(c) non-monotonicity for q=5 is a VELOCITY effect**, not a Casimir oscillation. The Casimir product vc itself is monotonic.

4. **Detrended vc residuals are 8-21× larger for q>4 than q=2.** Consistent with walking-enhanced corrections but not diagnostic of oscillatory vs monotonic.

5. **Richardson extrapolation drifts for q>4.** This may indicate log corrections or the beginning of an oscillation whose period is too long to resolve.

**What would be needed to detect Im(c) oscillations:**
- DMRG Casimir extraction at N=10-100 with <10⁻⁶ precision
- Or: non-Hermitian formulation (as in arXiv:2502.02001) where complex c enters directly
- Our exact diag sizes are fundamentally too small for this test

**Not novel.** Negative result, but informative for scoping future work.
