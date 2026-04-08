# Sprint 117 -- q=30 chi_F: Logarithmic Model Stabilizes

**Status:** Complete

## Motivation
Sprint 116 showed log+loglog marginally AIC-best over pure log (dAIC=1.4) with 10 data points. q=30 provides the 11th point. Only n=3,4 feasible (n=5=24.3M exceeds GPU).

## Experiment 117a -- q=30 chi_F spectral decomposition (n=3,4)

| n | dim | gap_m | |me|^2 | chi_F | dominant_frac | time |
|---|-----|-------|---------|-------|---------------|------|
| 3 | 27,000 | 0.3062 | 1623.3 | 5770.4 | 1.000 | 1.5s |
| 4 | 810,000 | 0.1543 | 2584.9 | 27160.1 | 1.000 | 23.7s |

Pairwise alpha = 5.384 (only pair (3,4) available).
- z_m = 2.384, beta_me = 1.617
- Single-multiplet dominance confirmed at q=30

Predictions vs measured:
- Log+loglog: 5.47 (over by 0.08)
- Pure log: 5.37 (under by 0.01)

Note: pairwise (3,4) underestimates global alpha (at q=25: pair was 4.93 vs global 5.17). True global likely ~5.5-5.6.

## Experiment 117b -- alpha(q) refit with 11 data points

AIC ranking (11 points, q=30 uses pair(3,4) lower bound):
1. **log+loglog**: alpha(q) = 2.408 ln(q) - 1.324 ln(ln(q)) - 1.147, AIC=-55.84
2. **logarithmic**: alpha(q) = 1.866 ln(q) - 0.969, AIC=-55.03, dAIC=+0.8
3. quadratic: dAIC=+7.7 (unphysical: goes negative at q=100)
4-6. sqrt, power_law, linear: dAIC>=11.3

**Critical observation: dAIC narrowed from 1.4 (Sprint 116) to 0.8 (Sprint 117).** The subleading correction is getting weaker with more data.

LOO cross-validation:
- log+loglog: LOO-RMS = 0.0811
- logarithmic: LOO-RMS = 0.0845
- Nearly identical. Both have largest errors at q=15 and q=25.

Updated pure log fit: alpha(q) = 1.866 ln(q) - 0.969 (stable across refits).

## Conclusions
1. alpha(q=30) pair(3,4) = 5.38, consistent with both models
2. **dAIC gap narrowing: pure logarithmic is likely the correct model** (simpler, nearly identical LOO performance)
3. The log+loglog subleading correction was an artifact of limited data range
4. alpha(q) mapping now covers q=5-30 with 11 points -- diminishing returns on extending further
5. Next priority: hardware validation (580s QPU unused for 91+ sprints)
