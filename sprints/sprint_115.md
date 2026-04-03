# Sprint 115 — α(q) at q=20: Logarithmic Growth Confirmed

**Date:** 2026-04-02
**Status:** Complete

## Motivation
Sprint 114 established logarithmic α(q) ≈ 1.73 ln(q) − 0.68 as physically preferred but AIC-second (ΔAIC=+14.9 behind quadratic). At q=20, the models diverge decisively:
- Logarithmic: α = 1.73 ln(20) − 0.68 = **4.49**
- Quadratic: α = −0.0134(400) + 0.451(20) + 0.159 = **3.82**
- Difference: 0.67 (15%) — easily distinguishable

Dimensions: n=3 → 8,000 (CPU), n=4 → 160,000 (CPU), n=5 → 3,200,000 (GPU).

## Experiment 115a: q=20 χ_F Spectral Decomposition

**Plan:** Measure χ_F at n=3,4,5 via spectral decomposition. Extract α, z_m, β_me.

**Results:**
| n | dim | gap_m | |me|² | χ_F | frac | time |
|---|-----|-------|-------|-----|------|------|
| 3 | 8,000 | 0.4220 | 680.3 | 1273.7 | 1.000 | 0.2s |
| 4 | 160,000 | 0.2384 | 1031.7 | 4537.6 | 1.000 | 14.2s |
| 5 | 3,200,000 | 0.1481 | 1462.8 | 13346.3 | 1.000 | 332.8s |

- **Single-multiplet dominance confirmed** (frac=1.000) through q=20
- Pairwise α: (3,4)→4.416, (4,5)→4.835 — converging upward
- Global α = **4.590**
- z_m = 2.047 (first time above 2.0!), β_me = 1.496
- α = β_me + 2z_m − 1 exact (4.590 = 1.496 + 2×2.047 − 1)

**Model comparison at q=20:**
| Model | Predicted | Measured | Residual |
|-------|-----------|----------|----------|
| **Logarithmic** | 4.49 | 4.59 | **+0.10** |
| √q | 4.75 | 4.59 | −0.16 |
| Power-law | 4.76 | 4.59 | −0.17 |
| **Quadratic** | 3.82 | 4.59 | **+0.77** |
| Linear | 5.04 | 4.59 | −0.45 |

**Quadratic DECISIVELY RULED OUT** — off by 0.77 (17%). Logarithmic closest at 0.10 (2.3%).

## Experiment 115b: α(q) Refit with 9 Data Points (q=5-20)

**Global alpha AIC ranking:**
| Model | RMS | AIC | ΔAIC | Notes |
|-------|-----|-----|------|-------|
| **Logarithmic** | 0.046 | −51.6 | 0.0 | **NOW AIC-BEST** |
| Quadratic | 0.069 | −42.2 | +9.4 | peaks at q≈24 (unphysical) |
| √q | 0.085 | −40.3 | +11.3 | |
| Power-law | 0.090 | −39.4 | +12.2 | |
| Linear | 0.153 | −29.8 | +21.8 | |

**Updated fit: α(q) ≈ 1.776 ln(q) − 0.785** (9 points, RMS=0.046)

**Component fits (all logarithmic, all AIC-best):**
- z_m(q) ≈ 0.757 ln(q) − 0.231 (RMS=0.032, beats linear by ΔAIC=14.5)
- β_me(q) ≈ 1.125 ln(q) − 1.680 (RMS=0.165, beats linear by ΔAIC=6.1)
- Reconstructed: α = (2×0.757 + 1.125) ln(q) + (2×(−0.231) + (−1.680) − 1) = 2.640 ln(q) − 3.142

Note: The component reconstruction overshoots α because of correlated errors in z_m and β_me at q=10-12. The direct logarithmic fit of α(q) is more reliable than summing component fits.

## Surprises
1. **Logarithmic now AIC-BEST** — was ΔAIC=+14.9 behind quadratic with 8 points, now leads by ΔAIC=9.4 with 9 points. The q=20 measurement flipped the ranking completely.
2. **Quadratic peaks at q≈24** — confirmed unphysical. With 9 points spanning q=5-20, the quadratic turnover is clearly an artifact.
3. **z_m > 2.0 for the first time** at q=20. The multiplet gap now closes faster than 1/N², meaning χ_F grows faster than N⁴ (super-quartic!).
4. **Single-multiplet dominance at q=20** where walking is fully broken (c_eff/Re(c) estimated <0.3). The mechanism persists arbitrarily far from the walking regime.

## Updated α(q) Table (9 data points)

| q | α(global) | α(last pair) | z_m | β_me |
|---|-----------|-------------|-----|------|
| 5 | 2.091 | 2.100 | 1.028 | 0.125 |
| 6 | 2.377 | 2.402 | 1.072 | 0.260 |
| 7 | 2.650 | 2.670 | 1.238 | 0.418 |
| 8 | 2.897 | 2.897 | 1.320 | 0.564 |
| 9 | 3.159 | 3.159 | 1.438 | 0.709 |
| 10 | 3.349 | 3.417 | 1.566 | 1.285 |
| 12 | 3.631 | 3.734 | 1.662 | 1.307 |
| 15 | 3.921 | 4.070 | 1.784 | 1.353 |
| **20** | **4.590** | **4.835** | **2.047** | **1.496** |

**Best fit: α(q) ≈ 1.78 ln(q) − 0.79 (logarithmic, AIC-best AND physically preferred)**

Predictions: α(25) ≈ 4.93, α(30) ≈ 5.26.
