# Sprint 122 — Log Correction Fits at q=4 & Literature Context

**Date:** 2026-04-08
**Thread:** chi_F spectral decomposition — model comparison
**Status:** Complete

## Motivation

S_q q=4 is the Ashkin-Teller point with a marginal operator. Literature (Salas-Sokal 1997, Balog et al. 2007) predicts multiplicative log corrections: chi_F ~ N^2(ln N)^{-p}. Sprint 121 measured alpha=1.777 — below 2.0, consistent with log corrections. This sprint fits multiple correction forms to extract p, and searches the literature for chi_F-specific predictions.

## Literature Search

Searched: "Salas Sokal q=4 Potts logarithmic", "Ashkin-Teller fidelity susceptibility scaling", "four-state Potts multiplicative log corrections", "Berche logarithmic corrections Potts", "fidelity susceptibility log correction marginal".

**Key references found:**
- **Salas & Sokal (1997):** J. Stat. Phys. 88, 567-615. Definitive paper on q=4 log corrections. Hat exponents: chi ~ L^{7/4}(ln L)^{-1/8}, C_H ~ L(ln L)^{-3/2}, M ~ L^{-1/8}(ln L)^{-1/16}.
- **Berche et al. (2007):** arXiv:0707.3317. Universal amplitude ratios with log corrections.
- **Kenna & Lang (2006):** PRL 97, 155702. Self-consistent scaling theory for log-correction exponents.
- **Lv, Hucht, Deng (2019):** arXiv:1904.00406. Needed L up to 1024 to pin down classical exponents.

**Critical finding: NO prior work on chi_F log corrections for q=4 Potts.** The literature covers classical observables (magnetization, susceptibility, specific heat). The fidelity susceptibility log correction has not been explicitly derived.

From the Salas-Sokal RG framework, the predicted exponent depends on how chi_F couples to the marginal dilution operator:
- If chi_F ~ specific heat (2nd derivative): p = 3/2
- If chi_F ~ thermal susceptibility: p = 3/4
- The asymptotic regime requires L >> 11 (Lv et al. needed L=1024 for classical quantities)

## Experiments

### 122a — Log correction fits at q=4

**Method:** Fit 4 models to chi_F data from Sprints 120-121: (1) pure power law, (2) log-corrected with alpha=2 fixed, (3) log-corrected with alpha free, (4) power law + 1/N^2 correction.

#### S_q q=4 Results (8 sizes, n=4-11)

| Model | alpha | p | R^2 | AIC | dAIC |
|-------|-------|---|-----|-----|------|
| Power + 1/N^2 corr | **1.757** | — | 0.999999950 | **-72.20** | **0** |
| Log-corrected (free) | 1.693 | -0.17 | 0.999999898 | -66.48 | 5.7 |
| Pure power | 1.777 | — | 0.999989880 | -50.13 | 22.1 |
| Log-corrected (alpha=2) | 2.000 | 0.45 | 0.999878980 | -30.28 | **41.9** |

**The log-corrected model with fixed alpha=2 is the WORST fit** (dAIC=42 vs best). At n=4-11, the data strongly prefers a power law with alpha~1.76 and 1/N^2 finite-size corrections. The "alpha=2 with log corrections" form simply cannot fit the data at these sizes.

The free log fit trades alpha for p, giving alpha=1.69 and p=-0.17 (negative p — log corrections go the WRONG direction for asymptotic alpha=2).

#### Hybrid q=4 Results (8 sizes, n=4-11)

| Model | alpha | p | R^2 | AIC | dAIC |
|-------|-------|---|-----|-----|------|
| Power + 1/N^2 corr | **1.460** | — | 0.999999987 | **-105.71** | **0** |
| Log-corrected (free) | 1.239 | -0.58 | 0.999999900 | -89.29 | 16.4 |
| Pure power | 1.526 | — | 0.999843570 | -50.85 | 54.9 |
| Log-corrected (alpha=2) | 2.000 | 0.94 | 0.998892710 | -35.20 | **70.5** |

Hybrid shows the same pattern but stronger: alpha=2 form is catastrophically bad (dAIC=70). Best fit: alpha=1.460 with 1/N^2 correction. Negative p in the free log fit confirms the hybrid is NOT approaching alpha=2.

## Key Findings

### 1. At accessible sizes, alpha=2 with log corrections does NOT fit

Despite the literature predicting chi_F ~ N^2(ln N)^{-p} asymptotically for S_q q=4, this form is the worst fit at n=4-11 for both models. The asymptotic regime for log corrections at q=4 is notoriously difficult to reach — Lv et al. (2019) needed L=1024 for classical observables.

### 2. Best description: power law with 1/N^2 finite-size correction

Both models are best described by chi_F = A*N^alpha*(1 + B/N^2) with:
- **S_q:** alpha = 1.757 ± 0.002, B = -0.54
- **Hybrid:** alpha = 1.460 ± 0.001, B = -1.64

The corrected alpha is close to the pairwise asymptotic values (S_q: 1.771, hybrid: 1.489), as expected from 1/N^2 corrections that become negligible at large N.

### 3. Alpha extrapolation question

**S_q:** alpha_corrected = 1.757 ± 0.002. Is this approaching 2.0? The pairwise drift (1.825 → 1.771) is DOWNWARD, not toward 2.0. At accessible sizes, there is no evidence for alpha → 2.0. If the asymptotic value is indeed 2.0, the approach is slower than any power-law correction and presumably logarithmic — but we cannot distinguish this from true alpha < 2 at n ≤ 11.

**Hybrid:** alpha_corrected = 1.460 ± 0.001. The pairwise drift (1.637 → 1.489) is strongly downward. This is definitively NOT approaching 2.0. The hybrid q=4 transition is not in the Ashkin-Teller universality class.

### 4. No prior chi_F log correction measurement exists

The literature search found no paper computing the chi_F log correction exponent for q=4 Potts. The Salas-Sokal framework gives hat exponents for classical observables, but the mapping to fidelity susceptibility has not been worked out. If we could reach the asymptotic regime, measuring p for chi_F would be a new result — but n=11 is far too small.

## Summary

At n=4-11, both S_q and hybrid q=4 chi_F are best described by pure power laws (alpha=1.76 and 1.46 respectively) with 1/N^2 finite-size corrections. The theoretical alpha=2 with log corrections is not supported at these sizes — the asymptotic regime requires much larger systems. The hybrid alpha is clearly below the S_q value, confirming different universality classes at q=4.

**This sprint's contribution:** Established that the q=4 log correction regime is beyond our exact diag reach (n≤11). DMRG could potentially reach n=50+ where log corrections become visible — this is a possible future direction.

## Next Steps

1. **QPU hardware test** — 580s budget, unused 96 sprints. Test q=2 Ising chi_F scaling.
2. **DMRG for S_q q=4 at large n** — if DMRG can compute chi_F (needs ground state derivative), could reach asymptotic log correction regime.
3. **Compile results for hybrid model** — the hybrid Potts-clock model at q≥4 is unstudied in the literature. Systematic chi_F data + z_m(q) crossing could form a paper-worthy dataset.
