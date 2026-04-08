# Sprint 116 -- q=25 chi_F: Logarithmic alpha(q) with Subleading Correction

**Status:** Complete

## Motivation
Sprint 115 confirmed logarithmic alpha(q) = 1.78 ln(q) - 0.79 as AIC-best with 9 data points (q=5-20). Prediction for q=25: alpha = 4.94. This is the 10th data point testing whether logarithmic growth continues or saturates.

Dimensions: n=3 (15.6k), n=4 (390k), n=5 (9.8M GPU). All feasible but n=5 pushes GPU limits.

## Experiment 116a -- q=25 chi_F spectral decomposition

**Result:** alpha(q=25) = 5.170 (global), 5.504 (last-pair)

| n | dim | gap_m | |me|^2 | chi_F | dominant_frac | time |
|---|-----|-------|---------|-------|---------------|------|
| 3 | 15,625 | 0.3544 | 1099.0 | 2916.8 | 1.000 | 0.6s |
| 4 | 390,625 | 0.1885 | 1713.6 | 12053.8 | 1.000 | 22.9s |
| 5 | 9,765,625 | 0.1100 | 2489.5 | 41164.2 | 1.000 | 1616s |

- z_m = 2.286 (global), continuing to grow past 2.0
- 100% single-multiplet dominance confirmed at q=25
- Spectral mechanism alpha = beta_me + 2*z_m - 1 verified exactly (pair diff = 0.000000)

Logarithmic prediction was 4.94, measured 5.17 -- residual +0.23 (4.7% underprediction). Suggests subleading correction.

## Experiment 116b -- alpha(q) refit with 10 data points

**Key result: log+loglog model is now AIC-best.**

AIC ranking (global alpha, 10 points):
1. **log+loglog**: alpha(q) = 2.616 ln(q) - 1.773 ln(ln(q)) - 1.258, AIC=-50.13
2. logarithmic: alpha(q) = 1.863 ln(q) - 0.962, AIC=-48.73, dAIC=+1.4
3. sqrt: AIC=-45.90, dAIC=+4.2
4. power_law: AIC=-45.00, dAIC=+5.1
5. quadratic: AIC=-42.45, dAIC=+7.7 (UNPHYSICAL: peaks at q=32.6)
6. linear: AIC=-31.86, dAIC=+18.3

The dAIC=1.4 between log+loglog and pure log is not decisive (need dAIC>4 for strong preference). Both remain viable. The subleading ln(ln(q)) correction is physically plausible but requires more data to confirm.

Component fits:
- z_m(q) = 0.786 ln(q) - 0.290 (logarithmic AIC-best)
- beta_me(q): logarithmic best but noisier (RMS=0.177)

Predictions:
| Model | q=30 | q=50 |
|-------|------|------|
| log+loglog | 5.47 | 6.56 |
| logarithmic | 5.37 | 6.32 |

## Conclusions
1. alpha(q=25) = 5.17 confirms continued growth beyond previous range
2. Pure logarithmic slightly underpredicts at q=25 (+4.7%)
3. log+loglog correction marginally preferred (dAIC=1.4) but not decisive
4. z_m crosses 2.0 for first time at q=20, reaches 2.29 at q=25
5. Single-multiplet dominance universal through q=25
6. n=5 at q=25 (9.8M dim) took 1616s -- approaching GPU memory limit
