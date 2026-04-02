# Sprint 075 — q=2 Ly=5 Cylinder & Cylinder Entanglement Entropy

**Date:** 2026-04-02
**Thread:** Cylinder dimensional crossover — completing q=2, new entropy observable

## Motivation

The q=2 Ly convergence series (Ly=2,3,4) shows decelerating approach to 2D: 38.6% → 77.7% → 84.0%. Exponential fit predicted Ly=5: g_c≈0.745 (95% to 2D). Additionally, entanglement entropy has never been measured on cylinder geometries — does c_eff interpolate between 1D and 2D?

## Experiments

### 075a — q=2 Ly=5 Cylinder Gap Crossings
- **Goal:** Measure g_c(q=2, Ly=5) via gap×Lx crossings
- **Prediction:** g_c ≈ 0.745 (exponential fit)
- **Status:** COMPLETE

**Results:** g_c(Ly=5) = 0.701 from (Lx=3, Lx=4) crossing. 86.6% to 2D. Lx=5 (dim=33M) takes 541s/pt — infeasible for full scan.

**Exponential prediction OVERSHOOTS by 5.9%.** Ly=4→5 increment is only +0.013, much smaller than Ly=3→4 (+0.033). Convergence decelerates faster than exponential.

| Ly | g_c | Progress | Increment |
|----|-----|----------|-----------|
| 1 | 0.250 | 0.0% | — |
| 2 | 0.451 | 38.6% | +0.201 |
| 3 | 0.655 | 77.7% | +0.204 |
| 4 | 0.688 | 84.1% | +0.033 |
| 5 | 0.701 | 86.6% | +0.013 |
| 2D | 0.771 | 100% | — |

### 075b — Entanglement Entropy on q=2 Cylinders at Criticality
- **Goal:** Measure S(Lx/2) at g_c for Ly=1,2,3,4 cylinders
- **Status:** COMPLETE

**Results:** c_eff extracted from S ~ (c/6)·ln(Lx) scaling at each Ly:

| Ly | c_eff | Lx values | S range |
|----|-------|-----------|---------|
| 1 | 0.610 | 4,6,8,10,12 | 0.28-0.40 |
| 2 | 0.752 | 4,6,8 | 0.31-0.40 |
| 3 | 0.817 | 4,6 | 0.32-0.38 |
| 4 | N/A (1 pt) | 4 | 0.36 |

c_eff grows with Ly: 0.61→0.75→0.82. BUT this is a finite-size overshoot effect, NOT approach to higher c. 2D Ising also has c=0.5. The area-law contribution from the Ly bonds being cut masks the logarithmic term at small Lx, inflating the apparent c_eff.

**S per bond cut DECREASING:** 0.285 (Ly=1) → 0.157 (Ly=2) → 0.107 (Ly=3) → 0.089 (Ly=4). Individual bonds become less entangled as Ly grows — approaching area-law behavior expected in 2D.

### 075c — Convergence Analysis
- **Goal:** Compare exponential vs power-law fits; cross-analyze g_c and entropy
- **Status:** COMPLETE

**Power-law fits BETTER than exponential:**
- Exponential: g_c = 0.771 - 1.316·exp(-Ly/1.383), RMS = 0.025
- Power-law: g_c = 0.771 - 1.285/Ly^2.03, RMS = 0.016

The ~1/Ly² convergence is consistent with the expected finite-width correction for quasi-1D systems approaching 2D. The exponential model from Sprint 073 is superseded — it was fit with fewer points.

**Updated predictions (power-law):**
| Ly | g_c (predicted) | Progress |
|----|----------------|----------|
| 6 | 0.737 | 93.5% |
| 8 | 0.752 | 96.4% |
| 10 | 0.759 | 97.7% |

## Key Findings

1. **g_c(q=2, Ly=5) = 0.701** — 86.6% to 2D, Ly=4→5 increment only +0.013
2. **Power-law 1/Ly² convergence, NOT exponential** — exponential overestimates convergence rate. RMS 0.016 vs 0.025.
3. **c_eff grows on cylinders (0.61→0.82)** but this is FSS overshoot from area-law contamination, not a physical c increase
4. **S per bond cut decreases monotonically** (0.285→0.089) — area-law regime emerging

## Surprises

- Exponential prediction overshoots by 5.9% (0.745 vs 0.701) — convergence slowing faster than expected
- Ly=1→2 and Ly=2→3 increments are nearly identical (+0.201, +0.204) — first two steps are linear!
- Abrupt deceleration at Ly=3→4 (increment drops 6x) — possible crossover from linear to power-law regime
- Lx=5 at Ly=5 (dim=33M) takes 541s/pt — beyond exact diag feasibility
- c_eff at Ly=1 (0.61) already 22% above exact c=0.5 — small-Lx overshoot is significant

[Results: results/sprint_075a_cylinder_ly5_q2.json, results/sprint_075b_cylinder_entropy.json, results/sprint_075c_convergence.json]
