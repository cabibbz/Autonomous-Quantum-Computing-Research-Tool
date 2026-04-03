# Sprint 111 — q=10 χ_F Measurement: Linear vs Sublinear α(q) Discrimination

**Date:** 2026-04-02
**Status:** Complete (2 experiments)

## Motivation

Sprint 110 revised α(q) from 0.315q+0.469 to 0.26q+0.81 using walking-only data (q=5-9). Two fits are competitive:
- **Linear:** α = 0.262q + 0.815 (RMS=0.019)
- **Quadratic:** α = -0.009q² + 0.389q + 0.386 (RMS=0.011, ΔAIC=3.3)

At q=10, these predict: linear α≈3.44, quadratic α≈3.28. Difference ~5%.

**Literature check:** No prior χ_F scaling at q=10 S_q Potts. Gorbenko/Rychkov/Zan complex CFT predicts Re(c)=1.584 for q=10. Walking fully broken (c_eff/Re(c)=0.60, Sprint 080).

## Experiment 111a: q=10 χ_F Spectral Decomposition (n=5,6,7)

Dimensions: n=5 (100k, 3s), n=6 (1M, 49s), n=7 (10M GPU, 700s).

| n | dim | gap_m | |me|² | χ_F | frac | time |
|---|-----|-------|-------|-----|------|------|
| 5 | 100,000 | 0.3284 | 269.6 | 500.1 | 1.0000 | 3.2s |
| 6 | 1,000,000 | 0.2482 | 337.1 | 911.9 | 1.0000 | 48.5s |
| 7 | 10,000,000 | 0.1950 | 411.0 | 1544.2 | 1.0000 | 699.5s |

**Key findings:**
- **Single-multiplet dominance universal through q=10** (frac=1.0000 at all sizes)
- Pairwise α: (5,6)→3.295, (6,7)→3.417 — **converging upward** (drift +0.122)
- Global α = 3.349
- z_m = 1.566 (gap closing faster than any prior q)
- β_me = 1.285 (matrix element growth faster than any prior q)
- Decomposition exact: α = β_me + 2z_m - 1 to machine precision

**Comparison to predictions:**
- Linear (0.262q+0.815) predicts 3.435 → residual -0.086
- Quadratic (-0.009q²+0.389q+0.386) predicts 3.376 → residual -0.027
- q=10 last-pair (3.417) is converging toward linear prediction

**Surprise:** n=7 took 700s (over 300s limit). The 10M build (68s) + GPU eigsh (630s) pushed timing. q=10 n=7 is at the practical limit.

## Experiment 111b: α(q) Refit with q=5-10 (6 data points)

### Model comparison (AIC):

| Model | Formula | RMS | AIC | ΔAIC |
|-------|---------|-----|-----|------|
| Power-law | 0.690·q^0.694 | 0.0125 | -48.56 | 0.00 |
| Sqrt | 1.404·√q − 1.045 | 0.0157 | -45.87 | 2.69 |
| Quadratic | -0.0046q²+0.329q+0.579 | 0.0133 | -45.82 | 2.74 |
| Linear | 0.260q+0.827 | 0.0177 | -44.44 | 4.12 |

**Power-law α(q) = 0.69·q^0.69 is AIC-preferred.** The sublinear exponent (0.69) means α grows slower than linear. At q=10, all models give similar values (range 3.35-3.44) but diverge at larger q:
- q=20: linear→6.0, power-law→4.8, sqrt→5.2

### Component fits (q=5-10):
- z_m(q) = 0.083q + 0.734 (Sprint 110: 0.082q + 0.741 — stable)
- β_me(q) = 0.094q + 0.360 (Sprint 110: 0.098q + 0.333 — stable)
- Reconstructed α = 0.260q + 0.827 exactly matches direct linear fit

### Convergence analysis:
All q values show **upward convergence** in pairwise α:
- q=5: drift +0.009/pair (well-converged, 4 pairs)
- q=6: drift +0.022/pair (converging, 3 pairs)
- q=7: drift +0.042/pair (converging, 2 pairs)
- q=10: drift +0.122/pair (early convergence, 2 pairs)

Drift increases with q — larger q values are further from asymptotic α. The last-pair α for q=10 (3.417) is a lower bound.

## Conclusions

1. **α(q=10) ≈ 3.42** (last-pair, lower bound), confirming super-first-order χ_F scaling continues to q=10
2. **Single-multiplet dominance is universal** through q=10 (frac=1.000 at all sizes for all q=5-10)
3. **Power-law α(q) = 0.69·q^0.69 preferred over linear** by ΔAIC=4.1, but the discrimination is weak — all models fit comparably (RMS 0.013-0.018)
4. **Component decomposition stable**: z_m and β_me linear fits unchanged from Sprint 110. The linear reconstruction of α is exact
5. **Walking mechanism persists even where walking breaks down**: χ_F super-scaling (α>2) is present for ALL q≥5 regardless of c_eff/Re(c) walking status. This is the single-multiplet mechanism, independent of the entropy-based walking phenomenon.

## Sprint report files
- `exp_111a_q10_chif.py` — q=10 spectral decomposition
- `exp_111b_alpha_refit_q10.py` — Full α(q) refit
- `results/sprint_111a_q10_chif.json`
- `results/sprint_111b_alpha_refit_q10.json`
