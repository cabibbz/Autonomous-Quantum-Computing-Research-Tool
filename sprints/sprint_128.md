# Sprint 128 — Fill q=6 gap in alpha(q) + asymptotic extrapolation

## Goal
1. Compute S_q q=6 exact chi_F at sizes n=4-9 (g_c=1/6, up to 10M dim GPU)
2. Find hybrid q=6 g_c via gap*N crossing, then compute chi_F
3. Asymptotic extrapolation: fit alpha_eff(N) = alpha_inf + c/N^p for all q

## Results

### 128a: S_q q=6 exact chi_F
Six sizes computed (n=4-9, max dim 10M):

| n | dim | chi_F | time |
|---|-----|-------|------|
| 4 | 1,296 | 63.5657 | 0.0s |
| 5 | 7,776 | 107.3445 | 0.0s |
| 6 | 46,656 | 164.7249 | 0.4s |
| 7 | 279,936 | 237.0861 | 2.4s |
| 8 | 1,679,616 | 325.8047 | 7.5s |
| 9 | 10,077,696 | 432.2600 | 55.7s |

**alpha = 2.375 +/- 0.006** (R^2 = 0.99998)

Pairwise exponents *increasing*: 2.348 -> 2.349 -> 2.362 -> 2.381 -> 2.400. Same trend as q=5,7 (walking regime, alpha grows with system size).

### 128b: Hybrid q=6 g_c + chi_F
Gap*N crossing scan at n=6,8 with 25 g-values in [0.44, 0.54]:
**g_c(hybrid, q=6) = 0.474** (between g=0.473 and g=0.478)

Chi_F at g_c=0.474 for n=4-8:

| n | chi_F |
|---|-------|
| 4 | 3.4657 |
| 5 | 4.7029 |
| 6 | 5.8891 |
| 7 | 7.0114 |
| 8 | 8.0660 |

**alpha = 1.186 +/- 0.038** (R^2 = 0.9974)

Pairwise exponents *decreasing*: 1.368 -> 1.234 -> 1.132 -> 1.049. Heading toward alpha~1 or below. Consistent with hybrid q>=5 being continuous (non-walking).

### 128c: Asymptotic extrapolation
Fit alpha_eff(N) = alpha_inf + c/N^p to pairwise exponent data.

**Reliable extrapolations:**
- S_q q=3: alpha_inf = 1.412 +/- 0.002 (exact: 7/5=1.400, 0.8% error -- validates method)
- Hybrid q=3: alpha_inf = 1.408 +/- 0.001 (matches S_q q=3, expected since models equal at q=3)
- Hybrid q=4: alpha_inf = 1.515 +/- 0.021

**Unreliable extrapolations** (non-monotonic or too few points):
- S_q q=4: noisy pairwise (1.85, 1.80, 1.84, 1.72, 1.79) -- no clear trend
- S_q q>=5: pairwise increasing, extrapolation diverges upward
- Hybrid q>=5: extrapolation hits zero (alpha_inf = 0), meaning chi_F growth slower than any power law

### alpha(q) refit with q=6

**S_q:** alpha(q) = 1.337 * ln(q) - 0.023
- chi2/dof = 5.0 (improved from 22 without q=6)
- q=6 sits exactly on the fit curve (residual +0.003)
- Quadratic-in-ln(q) slightly better: chi2/dof = 3.7

| q | alpha | predicted | residual |
|---|-------|-----------|----------|
| 3 | 1.468 | 1.446 | +0.022 |
| 4 | 1.794 | 1.830 | -0.036 |
| 5 | 2.139 | 2.129 | +0.011 |
| 6 | 2.375 | 2.372 | +0.003 |
| 7 | 2.584 | 2.578 | +0.006 |

**Hybrid:** alpha(q) clearly non-log. Quadratic in ln(q) fits well:
alpha(q) = -1.532*ln(q)^2 + 4.013*ln(q) - 1.076 (chi2/dof = 1.5)
Peaks at q~3.7, decreases for q>4 toward zero.

## Key Insights

1. **S_q alpha(q) is very close to logarithmic.** q=6 fills the gap and improves the fit dramatically. The formula alpha(q) ~ 1.34*ln(q) works to 2-3% for q=3-7.

2. **Hybrid transitions collapse at large q.** Pairwise alpha at hybrid q=6 already at 1.05 and falling. Extrapolation gives alpha_inf=0 for q>=5. This means chi_F grows slower than any power law -- logarithmically? This is consistent with the hybrid q>=5 transitions being genuinely continuous (BKT-like slow divergence rather than power-law walking).

3. **S_q q=4 remains stubborn.** Pairwise alpha oscillates around 1.79 with no clear monotonic drift. Neither power-law nor log-corrected models fit well. This is the BKT point where corrections are exponentially slow.

4. **Extrapolation validates at q=3.** Getting alpha_inf=1.412 vs exact 1.400 (0.8% error) gives confidence in the method for cases where it converges.
