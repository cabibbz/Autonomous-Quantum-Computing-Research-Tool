# Sprint 112 — q=12 χ_F Spectral: Quadratic α(q) Strongly Preferred

**Status:** Complete (2 experiments)

## Motivation
Sprint 111 found power-law α(q)=0.69·q^0.69 slightly preferred over linear (ΔAIC=4.1), but discrimination was weak at q≤10. At q=12, models diverge: linear predicts 3.95, power-law 3.83, quadratic 3.76. Testing which functional form survives.

## 112a — q=12 χ_F spectral decomposition

**q=12 at n=4,5,6 (dims 21k, 249k, 3.0M).** GPU for n=6 (160s). g_c=1/12.

| n | dim | gap_m | |me|² | χ_F | frac | time |
|---|-----|-------|-------|-----|------|------|
| 4 | 20,736 | 0.3907 | 318.1 | 520.8 | 1.000 | 0.6s |
| 5 | 248,832 | 0.2710 | 422.6 | 1150.9 | 1.000 | 13s |
| 6 | 2,985,984 | 0.1991 | 540.6 | 2273.5 | 1.000 | 160s |

Pairwise α: (4,5)→3.553, (5,6)→3.734 — **converging upward** (consistent with all prior q).
Global α = 3.631. z_m = 1.662, β_me = 1.307.

**Single-multiplet dominance (frac=1.000) universal through q=12.** 11-fold degenerate S₁₂ multiplet captures 100% of χ_F. Spectral gap symmetry-forbidden.

**All prior fits overshoot:**
- Linear (0.260q+0.827): 3.947 → residual -0.316
- Power-law (0.69q^0.69): 3.833 → residual -0.201
- Quadratic (-0.009q²+0.389q+0.386): 3.758 → residual -0.127

## 112b — α(q) refit with 7 data points (q=5-12)

Using last-pair pairwise α (best finite-size estimate):

| q | α (last-pair) |
|---|---------------|
| 5 | 2.100 |
| 6 | 2.402 |
| 7 | 2.670 |
| 8 | 2.897 |
| 9 | 3.159 |
| 10 | 3.417 |
| 12 | 3.734 |

**AIC model comparison:**

| Model | Formula | RMS | AIC | ΔAIC |
|-------|---------|-----|-----|------|
| **Quadratic** | **-0.0103q² + 0.411q + 0.302** | **0.020** | **-49.1** | **0** |
| Sqrt | 1.353√q − 0.914 | 0.026 | -46.9 | +2.2 |
| Logarithmic | 1.895·ln(q) − 0.990 | 0.033 | -43.8 | +5.3 |
| Power-law | 0.742·q^0.656 | 0.033 | -43.6 | +5.5 |
| Linear | 0.237q + 0.986 | 0.054 | -36.9 | +12.2 |

**Quadratic is strongly preferred (ΔAIC=12.2 over linear).** The q=12 data point decisively rules out linear growth. α(q) is sublinear and approaching saturation.

Predictions: quadratic gives α(q=15)≈4.15, α(q=20)≈4.40 — suggesting α may plateau near 4-5.

**Component analysis:** z_m linear fit is 0.098q+0.533 (RMS=0.035), β_me is noisier (RMS=0.12) — q=10,12 components may be unreliable due to only 2-3 sizes each, but overall α is more robust.

## Key Findings

1. **α(q) is sublinear** — quadratic with negative q² term, not linear. Prior Sprint 110 formula 0.262q+0.815 is WRONG at q≥12.
2. **Single-multiplet dominance universal through q=12.** (q-1)-fold degenerate S_q multiplet captures 100% of χ_F at all q tested.
3. **Walking super-scaling persists at q=12** where walking has fully broken down (c_eff/Re(c)≈0.51 extrapolated). Confirms χ_F mechanism is independent of entropy walking.
4. **Best functional form: α(q) ≈ -0.010q² + 0.41q + 0.30** or equivalently √q or ln(q) forms (all sublinear).

## Surprises
- Linear model decisively ruled out (ΔAIC=12.2) — was still viable at q=10 (ΔAIC=4.1)
- All five prior extrapolations overshoot at q=12
- √q form (ΔAIC=+2.2) is nearly as good as quadratic — suggests α(q) ~ √q may be the asymptotic form

[Full data: results/sprint_112a_q12_chif.json, results/sprint_112b_alpha_refit_q12.json]
