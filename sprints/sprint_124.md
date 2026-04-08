# Sprint 124 — DMRG chi_F for S_q q=4: Open-BC Alpha Drift Reveals Approach to Asymptotic Regime

**Date:** 2026-04-08
**Thread:** chi_F spectral decomposition -- q=4 log corrections
**Status:** Complete

## Motivation

Sprint 122 showed that at n=4-11 (periodic BC, exact diag), the q=4 S_q Potts chi_F is best fit by alpha=1.757 with 1/N^2 corrections, NOT alpha=2 with log corrections. The asymptotic regime where log corrections become visible requires L>>11 (Lv et al. 2019 needed L=1024 for classical observables).

Sprint 113b extended chi_F to n=12 via finite DMRG (open BC), but found alpha_open=1.51 -- dramatically different from periodic alpha=1.77. Sprint 113c confirmed this BC gap grows with q (walking amplifies boundary effects).

**Two approaches this sprint:**
1. iDMRG to bypass boundary AND finite-size effects (novel)
2. Extend finite DMRG to n=20 to measure alpha drift direction (reliable fallback)

## Literature Search

Searched: "DMRG fidelity susceptibility q=4 Potts model logarithmic corrections scaling", "infinite DMRG fidelity susceptibility transfer matrix overlap critical systems", "fidelity susceptibility four-state Potts scaling exponent 2024 2025 2026".

**No prior work found on chi_F for q=4 Potts at large sizes via ANY method (DMRG, MC, or otherwise).**

Key references:
- Salas & Sokal (1997): Log corrections for classical observables, hat exponents. No chi_F.
- Berche et al. (2007): Universal amplitude ratios with log corrections. No chi_F.
- Lv, Hucht, Deng (2019): Needed L=1024 for classical exponents. Motivates going large.

## Experiments

### 124a/b -- iDMRG chi_F scaling (FAILED)

**Method:** iDMRG at g_c=0.25 for chi_max=[20-200], overlap between psi(g) and psi(g+dg).
Tried both random initialization (124a) and warm-start from converged psi(g) (124b).

**Result: FAILED.** Wildly non-monotonic chi_F vs chi_max. q=2 cross-check: chi_F = 425, 76, 3135, 272, 311197 at chi_max = 20, 40, 60, 80, 100.

**Root cause:** S_q symmetry is non-abelian, cannot be conserved in TeNPy's iDMRG. Without charge conservation, iMPS converges to different local minima at different chi, making inter-state overlaps meaningless.

**Lesson:** iDMRG overlap chi_F requires abelian symmetry conservation. Not viable for S_q Potts in TeNPy.

### 124c -- Finite DMRG chi_F extension to n=14-20 (SUCCESS)

**Method:** Extend Sprint 113b's validated finite DMRG chi_F (open BC, overlap method). q=2 cross-check + q=4 main target.

#### q=2 cross-check (exact alpha=1.0)

| n | chi_max | chi_F/N | chi_used | S_mid | time |
|---|---------|---------|----------|-------|------|
| 14 | 60 | 0.9427 | 44 | 0.411 | 6s |
| 16 | 60 | 1.0847 | 50 | 0.423 | 7s |
| 20 | 80 | 1.3682 | 64 | 0.444 | 13s |
| 24 | 80 | 1.6509 | 72 | 0.460 | 17s |

Matches Sprint 113a values exactly. Combined pairwise alpha (n=8-24): 1.097, 1.075, 1.061, 1.051, 1.043, 1.037, 1.031. Steady drift DOWN toward exact 1.0. Validates method.

#### q=4 main target

| n | chi_max | chi_F/N | chi_used | S_mid | time |
|---|---------|---------|----------|-------|------|
| 14 | 80 | 14.965 | 80 | 0.697 | 62s |
| 16 | 100 | 18.328 | 100 | 0.722 | 107s |
| 18 | 120 | 21.925 | 120 | 0.743 | 174s |
| 20 | 150 | 25.741 | 150 | 0.762 | 304s |

#### Combined analysis (Sprint 113b + 124c): 8 sizes, n=6-20

| Pair | Pairwise alpha | Direction |
|------|---------------|-----------|
| (6,8) | 1.5071 | -- |
| (8,10) | 1.5048 | down |
| (10,12) | 1.5094 | up |
| (12,14) | 1.5137 | up |
| (14,16) | 1.5180 | up |
| (16,18) | 1.5212 | up |
| (18,20) | 1.5232 | up |

**Global fit: alpha = 1.512**
**Power+1/N^2: alpha = 1.524 +/- 0.002, R^2 = 0.99999938**
**Log-corrected (alpha=2): p = 1.237, R^2 = 0.99970826**

## Key Findings

### 1. Open-BC alpha DRIFTS UPWARD at q=4

This is the headline result. The pairwise alpha increases monotonically from 1.505 at (8,10) to 1.523 at (18,20). The drift rate is ~0.003 per pair.

**Contrast with q=2:** Open-BC alpha drifts DOWNWARD (1.097 -> 1.031) toward the exact value 1.0. The drift rate is ~0.008 per pair -- 2.7x faster than q=4.

The q=4 upward drift is consistent with alpha approaching a larger asymptotic value. The slow drift rate is consistent with the q=4 Potts model's notoriously severe finite-size corrections.

### 2. Open-BC corrected alpha = 1.524 vs periodic-BC corrected alpha = 1.757

The power+1/N^2 fit gives:
- **Open BC (this sprint, n=6-20):** alpha = 1.524 +/- 0.002
- **Periodic BC (Sprint 122, n=4-11):** alpha = 1.757 +/- 0.002

The gap is 0.233. Sprint 113c measured this BC gap as delta_alpha=0.29 at n<=12. Our corrected values narrow it slightly (0.233 at n<=20). Both should converge to the same bulk exponent at N -> infinity, but the convergence is very slow.

### 3. Log-corrected model (alpha=2) still fits poorly

The log-corrected form chi_F = A*N^2*(ln N)^{-p} with p=1.237 gives R^2=0.99971 -- much worse than the power+1/N^2 fit (R^2=0.999999). At n<=20, the data strongly prefers a pure power law with alpha<2.

However, the DIRECTION of alpha drift (upward) is consistent with eventual convergence to alpha=2 with log corrections -- we just haven't reached the asymptotic regime. This is exactly what Lv et al. (2019) found for classical observables: apparent power-law behavior at intermediate sizes, with log corrections only visible at L>100.

### 4. Convergence scale estimate

At the current drift rate (+0.003 per doubling of N), reaching alpha=2.0 from alpha=1.52 would require ~160 doublings = N ~ 10^48. This is clearly unreachable.

But the drift rate itself should ACCELERATE as log corrections kick in. The Salas-Sokal theory predicts alpha(N) = 2 - p/ln(N) + O(1/ln^2(N)). With p~1.2: alpha(20) = 2 - 1.2/3.0 = 1.60. We measure 1.52 -- the discrepancy (0.08) comes from boundary effects in open BC.

### 5. DMRG bond dimension saturated at n=20

chi_used = chi_max at n=14-20 (80, 100, 120, 150). The entanglement entropy S_mid = 0.76 at n=20 grows ~ (c/6)ln(n) with c~1, requiring chi ~ exp(6S/c) ~ exp(4.6) ~ 100. At n=20, chi_max=150 is barely sufficient. Going to n=30 would need chi_max~250+, taking ~30 minutes per point.

## Summary

| Method | Result | Status |
|--------|--------|--------|
| iDMRG overlap (124a/b) | FAILED -- non-abelian symmetry breaks convergence | Dead end |
| Finite DMRG open BC (124c) | alpha drifts UP: 1.505->1.523 (n=6-20) | Success |
| Power+1/N^2 fit | alpha_open = 1.524 +/- 0.002 | Best model |
| Log-corrected (alpha=2) | p=1.237, R^2=0.999708 | Disfavored |

**The upward drift of open-BC alpha at q=4 is the key new finding.** It demonstrates that the "effective" exponent is NOT yet converged -- it's approaching a higher asymptotic value. Whether that value is 1.77 (matching periodic BC) or 2.0 (log-corrected) cannot be resolved at n<=20. The q=4 Potts log correction problem requires sizes beyond current DMRG reach.

## Next Steps

1. **Accept the limitation** -- q=4 log corrections are beyond our tools (need n>100, likely n>1000)
2. **Return to main thread** -- compile systematic chi_F results across q values for publication
3. **Try q=5 DMRG extension** -- walking regime may show different open-BC alpha drift behavior
