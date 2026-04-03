# Sprint 114 — α(q) at q=15: Logarithmic Growth Emerging

**Date:** 2026-04-02
**Status:** Complete

## Motivation
Sprint 112 established quadratic α(q) as AIC-best (ΔAIC=2.2 over √q). At q=15, models diverge enough to discriminate: quadratic predicts 4.20, √q predicts 4.32. Dimensions: n=4 (51k CPU), n=5 (759k GPU).

## Experiment 114a: q=15 χ_F Spectral Decomposition

**Plan:** Measure χ_F at n=3,4,5 via spectral decomposition. Extract α.

**Results:**
| n | dim | gap_m | |me|² | χ_F | frac | time |
|---|-----|-------|-------|-----|------|------|
| 3 | 3,375 | 0.5247 | 365.1 | 442.1 | 1.000 | 0.1s |
| 4 | 50,625 | 0.3172 | 533.1 | 1324.8 | 1.000 | 2.2s |
| 5 | 759,375 | 0.2107 | 729.5 | 3285.1 | 1.000 | 56.4s |

- **Single-multiplet dominance confirmed** (frac=1.000) through q=15
- Pairwise α: (3,4)→3.815, (4,5)→4.070 — converging upward
- Global α = 3.921
- z_m = 1.784, β_me = 1.353
- α = β_me + 2z_m − 1 exact (3.921 = 1.353 + 2×1.784 − 1)

**ALL prior models overpredict:**
- Quadratic (−0.010q²+0.41q+0.30): predicted 4.20, residual −0.28
- √q (1.35√q−0.91): predicted 4.32, residual −0.40
- Log (2.50 ln(q)−2.60): predicted 4.17, residual −0.25
- Linear (0.260q+0.827): predicted 4.73, residual −0.81

## Experiment 114b: α(q) Refit with q=15

**8 data points** (q=5,6,7,8,9,10,12,15). Model comparison using global α:

| Model | RMS | AIC | ΔAIC | Params |
|-------|-----|-----|------|--------|
| **Quadratic** | 0.014 | −62.2 | 0.0 | −0.0134q²+0.451q+0.159 |
| Logarithmic | 0.041 | −47.3 | +14.9 | 1.726 ln(q)−0.684 |
| √q | 0.082 | −36.0 | +26.2 | 1.150√q−0.392 |
| Power-law | 0.088 | −34.9 | +27.3 | 0.900q^0.556 |
| Linear | 0.129 | −28.8 | +33.5 | 0.185q+1.343 |

**Quadratic overwhelmingly best** (ΔAIC=14.9 vs log, was only 2.2 in Sprint 112). But has a **physically suspect prediction**: the fitted quadratic peaks at q ≈ 17 then DECREASES. α(20)=3.83, α(25)=3.08 — this means the transition would become *less* singular at very large q, which is unphysical.

**Logarithmic is the best physically unbounded model**: α(q) ≈ 1.73 ln(q) − 0.68. Predicts α(20)=4.49, α(25)=4.87. Grows forever but slowly.

**Component analysis:**
- z_m(q) best fit: 0.749 ln(q) − 0.214 (RMS=0.034, beats linear RMS=0.059)
- β_me(q) best fit: 1.308 ln(q) − 2.046 (RMS=0.146, beats linear RMS=0.176)
- Both components are **logarithmic** in q, consistent with α(q) ∝ ln(q)

**Key physical interpretation:** The quadratic is a local approximation that captures curvature within q=5-15 but cannot extrapolate. The logarithmic form α(q) ∝ ln(q) is the most physically plausible functional form:
1. Both components (z_m and β_me) are individually logarithmic
2. It grows without bound (consistent with increasing singularity)
3. It's sublinear (consistent with all data)
4. No maximum or rollover (physical requirement)

## Surprises
- **All 4 prior models overpredict at q=15.** The sublinearity is *stronger* than expected. The q=12 data point happened to fall close to predictions; q=15 breaks them.
- **Quadratic fit peaks at q≈17** — this is an artifact of polynomial fitting, not physics. It strongly suggests logarithmic is the true asymptotic form.
- **Single-multiplet dominance persists at q=15** where walking is fully broken (c_eff/Re(c) estimated <0.4). The fidelity susceptibility mechanism is truly independent of the entropy-based walking regime.

## Updated α(q) Table (8 data points)

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

Best fit: α(q) ≈ 1.73 ln(q) − 0.68 (logarithmic, second-best AIC but physically preferred).
