# Sprint 046 — ν(q=4) Extraction: Marginal Point Shows BKT-like Behavior

## Motivation
We have ν(q) at three points: q=2 (ν=1, Ising exact), q=3 (ν=5/6, Potts exact), q=5 (ν≈2.0, Sprint 045). The ν(q) curve is non-monotonic — decreases from q=2→3, then increases sharply to q=5. But we're missing q=4, which is the *marginal* point where 2D classical Potts transitions from second-order to first-order. Key questions:

1. Where does ν(4) fall? Between 5/6 and 2? Or is the minimum at q=4, not q=3?
2. Does q=4 show logarithmic corrections (as expected from 2D classical marginality)?
3. Is the ν non-monotonicity smooth or does it jump at q=4?

## Method
- Direct MPS contraction for all-pairs MI (validated pipeline from Sprint 043)
- Three system sizes: n=8, 12, 16
- χ=20 (sufficient for d=4)
- Dense g sweep: 0.20 to 1.05
- Joint (ν, g_c) data collapse optimization + slope ratio analysis

## Experiments

### Exp 046a: q=4 Potts n=8 MI-CV sweep (direct MPS, χ=20)
*Status:* Complete (36.5s)

g=0.70→1.05: CV increases monotonically from 0.427 to 0.738. No crossing visible at single size.

| g | CV | E0 |
|---|----|----|
| 0.70 | 0.4270 | -13.41 |
| 0.75 | 0.4820 | -14.17 |
| 0.80 | 0.5329 | -14.94 |
| 0.85 | 0.5800 | -15.71 |
| 0.88 | 0.6067 | -16.18 |
| 0.90 | 0.6238 | -16.49 |
| 0.93 | 0.6487 | -16.96 |
| 0.95 | 0.6646 | -17.27 |
| 1.00 | 0.7027 | -18.05 |
| 1.05 | 0.7383 | -18.84 |

### Exp 046b: q=4 Potts n=12 MI-CV sweep (direct MPS, χ=20)
*Status:* Complete (73.0s)

n=12 CV systematically HIGHER than n=8 at all g from 0.70 to 1.05. Gap increases with g.

### Exp 046c: q=4 Potts n=16 MI-CV sweep + low-g extension (direct MPS, χ=20)
*Status:* Complete (118.3s)

n=16 > n=12 > n=8 at all g from 0.70 to 1.05. Extended to g=0.40: CV minimum at g≈0.40-0.50.

### Exp 046d: Low-g extension for n=8,12
*Status:* Complete

Found crossings at very low g:
- n=8,12: crossings near g≈0.29, 0.35, 0.41 (noisy, near CV minimum)
- n=12,16: crossing near g≈0.46
- n=8,16: crossing near g≈0.42

Crossings are 50-65% below expected g_c≈0.89 — far more shifted than q=5 (10-16% below).

### Exp 046e: Free data collapse
*Status:* Complete

Free collapse gives ν≈2.18, g_c≈0.53 (transition region) or ν≈3.0, g_c≈0.43 (broad region). But g_c is fitting to the CV minimum, NOT the phase transition.

### Exp 046f: Constrained collapse and slope analysis
*Status:* Complete

**Constrained collapse at g_c=0.89:** ν hits upper bound (10+), quality 10x worse than free fit. Standard FSS collapse fails because curves don't cross near g_c.

**Slope ratios near g_c≈0.89 (the key result):**
| Size pair | Slope ratio | 1/ν exponent | ν estimate |
|-----------|-------------|-------------|------------|
| n=8,12 | 1.235 | 0.520 | **1.92** |
| n=12,16 | 1.112 | 0.370 | **2.71** |
| n=8,16 | 1.371 | 0.455 | **2.20** |

**Critical: ν estimate INCREASES with system size (1.92 → 2.71).** This is the signature of logarithmic corrections or BKT-like behavior where ν is effectively infinite.

## Results

### MI-CV data consolidated (direct MPS, χ=20)

| g | n=8 | n=12 | n=16 | Δ(12-8) | Δ(16-12) |
|---|-----|------|------|---------|----------|
| 0.40 | 0.174 | 0.171 | 0.163 | -0.003 | -0.008 |
| 0.50 | 0.160 | 0.204 | 0.209 | +0.044 | +0.004 |
| 0.60 | 0.302 | 0.352 | 0.352 | +0.049 | -0.000 |
| 0.70 | 0.427 | 0.486 | 0.488 | +0.059 | +0.002 |
| 0.80 | 0.533 | 0.607 | 0.616 | +0.074 | +0.010 |
| 0.85 | 0.580 | 0.663 | 0.678 | +0.083 | +0.015 |
| 0.90 | 0.624 | 0.716 | 0.737 | +0.092 | +0.021 |
| 0.95 | 0.665 | 0.767 | 0.794 | +0.102 | +0.027 |
| 1.00 | 0.703 | 0.815 | 0.849 | +0.113 | +0.033 |
| 1.05 | 0.738 | 0.861 | 0.902 | +0.123 | +0.040 |

Key features:
- CV minimum near g≈0.40-0.50 for all sizes (maximally Democratic ordered phase)
- n=16 > n=12 > n=8 at all g above minimum — NO crossing near g_c
- Gap Δ(16-12) grows steadily: 0.002 (g=0.70) → 0.040 (g=1.05)
- Gap Δ(12-8) grows faster: 0.059 (g=0.70) → 0.123 (g=1.05)

### Derivative at g_c ≈ 0.89
| n | dCV/dg |
|---|--------|
| 8 | 0.855 |
| 12 | 1.050 |
| 16 | 1.180 |

## Conclusions

### 1. ν(q=4) ≥ 2.2, trending toward larger values — BKT-like behavior

The slope ratio method gives ν estimates that INCREASE with system size: 1.92 (n=8,12) → 2.71 (n=12,16). This is qualitatively different from q=2 (converges to 1), q=3 (converges to 5/6), and q=5 (converges to 2). If ν continues to grow, it signals BKT-like behavior where ξ ~ exp(c/√|g-g_c|) rather than power-law ξ ~ |g-g_c|^{-ν}.

### 2. Standard FSS collapse FAILS at q=4

Unlike q=2, q=3, and q=5 where data collapse with power-law scaling works, q=4 shows:
- Free collapse finds g_c=0.53 (fitting CV minimum, not phase transition)
- Constrained collapse at g_c=0.89 requires ν→∞
- No MI-CV curve crossing near g_c

This is consistent with q=4 being the marginal/BKT point of 1D quantum Potts.

### 3. ν(q) picture updated

| q | ν | Method | Nature |
|---|---|--------|--------|
| 2 | 1.0 | Ising exact | Standard second-order |
| 3 | 5/6 ≈ 0.833 | Potts exact | Standard second-order (minimum) |
| 4 | ≥2.2 (diverging?) | Slope ratio, trending up | **Marginal/BKT-like** |
| 5 | ~2.0 | Data collapse (Sprint 045) | Large-ν second-order |

The ν(q) curve has a minimum at q=3, then jumps sharply at q=4 (the marginal point). q=4 may have ν=∞ (BKT), with q=5 returning to finite ν=2.

### 4. Sprint 040 Gell-Mann data partially misleading

Sprint 040 reported MI-CV crossings at g_c≈0.893 for q=4 using Gell-Mann reconstruction. Direct MPS shows NO crossing near g_c — the curves are monotonically ordered (n16>n12>n8). The Gell-Mann method at d=4, while validated at d=2,3, appears to have small systematic errors at d=4 that created artificial crossings.

### Surprises
- ν(4) INCREASES with system size — first hint of BKT in the ν(q) series
- NO crossings near g_c — qualitatively different from all other q
- Gell-Mann reconstruction misleading even at d=4 (validated range was d≤4 but only tested at d=2,3)
- CV minimum at g≈0.4-0.5, far below g_c≈0.89 — ordered phase has complex MI structure
