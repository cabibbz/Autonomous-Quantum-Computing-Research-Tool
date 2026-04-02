# Sprint 073 — Cylinder Ly Convergence: q=3 Ly=3 & q=2 Ly=4

**Date:** 2026-04-02
**Status:** Complete

## Motivation

Sprint 072 established Ly=2 cylinder dataset for q=2,3,5 and Ly=3 for q=2 only. Two key gaps:
1. q=3 Ly=3 cylinder — does q=3 converge to 2D faster or slower than q=2?
2. q=2 Ly=4 cylinder — does the convergence accelerate? (77.7% at Ly=3, predict ~85-90% at Ly=4)

## Experiments

### 073a — q=3 Ly=3 Cylinder Gap Crossings

Dimensions: 3^{3·Lx} = 27^Lx. Lx=3: 19683, Lx=4: 531k, Lx=5: 14.3M (GPU).
Target: g_c(q=3, Ly=3) via gap×Lx crossings.
Expected: between Ly=2 (0.565) and 2D (1.267).

**Results:**
- Lx=5 (dim=14.3M): 251s/pt — infeasible for scan. Only Lx=3,4 available.
- Single crossing pair (3,4): **g_c = 0.797**
- gap×Lx at crossing: 2.187
- Ly convergence for q=3: 1D (0.333) → Ly=2 (0.565, 24.8%) → Ly=3 (0.797, 49.7%) → 2D (1.267)
- **q=3 converges MUCH slower than q=2** (49.7% vs 77.7% at Ly=3)

### 073b — q=2 Ly=4 Cylinder Gap Crossings

Dimensions: 2^{4·Lx} = 16^Lx. Lx=3: 4096, Lx=4: 65536, Lx=5: 1M, Lx=6: 16M.
Target: g_c(q=2, Ly=4) and Ly convergence rate.
Expected: ~85-90% progress toward 2D g_c=0.771.

**Results:**
- Lx=6 (dim=16.8M): 268s/pt — infeasible. Lx=3,4,5 scanned.
- Two crossing pairs: (3,4)=0.683, (4,5)=0.693. Mean **g_c = 0.688**.
- gap×Lx at crossings: 2.276 (3,4), 2.444 (4,5) — converging toward ~2.4 (consistent with Ly=3 value)
- Ly convergence for q=2: Ly=1 (0%), Ly=2 (38.6%), Ly=3 (77.7%), **Ly=4 (84.0%)**, 2D (100%)
- Exponential fit: g_c(Ly) = 0.771 - 1.678·exp(-Ly/1.20). Predicts g_c(Ly=5) = 0.745.
- **Convergence decelerating**: Ly=2→3 jump was +0.204, Ly=3→4 jump only +0.033.

### 073c — Cross-q Ly Convergence Analysis

Compare convergence rates for q=2 and q=3 as function of Ly.
Does q=3 approach 2D faster or slower?

**Results:**
- **q=3 converges 1.9x SLOWER than q=2** to 2D. Exponential decay length: B(q=3)=3.03 vs B(q=2)=1.60.
- q=3 at Ly=3: 49.7% of way to 2D. q=2 at Ly=3: 77.7%. At same Ly, q=3 is ~28% behind.
- q=3 has larger 2D/1D g_c ratio (3.80 vs 3.08) — more ground to cover AND slower convergence.
- Convergence ratio (deficit reduction per Ly step): q=2 irregular (0.61, 0.36, 0.72), q=3 steady (~0.70).
- Free g_c(2D) fit for q=3 from 3 Ly points is underdetermined (only 3 data points for 3 parameters).
- Exponential fit with known g_c(2D)=1.267 predicts q=3 reaches 80% at Ly≈5 (vs q=2 at Ly≈3.5).

## Key Findings

1. **q-dependent Ly convergence rate.** Larger q converges more slowly to 2D. The 1D→2D crossover requires more circumference for higher-q models. Decay length B ∝ q^α with α ≈ 1.6 (from 2 points).

2. **q=2 Ly=4: 84% to 2D.** Convergence decelerating as expected for exponential approach. The Ly=3→4 jump (+0.033) is 6x smaller than Ly=2→3 (+0.204). Gap×Lx at crossing converging toward ~2.4.

3. **q=3 needs larger cylinders.** At Ly=3, q=3 has only 1 crossing pair (Lx=3,4) due to dim=27^Lx scaling. Ly=4 would need 81^Lx: Lx=3 (531k), Lx=4 (43M) — borderline GPU feasible.

4. **Convergence is NOT universal in Ly.** Different q values have qualitatively different Ly-scaling. This is expected: correlation length ξ at criticality depends on q, and the Ly→2D crossover occurs when Ly ≫ ξ.

## Surprises

- q=3 at Ly=3 is only halfway to 2D — much slower than naive extrapolation from q=2
- q=2 Ly=4 convergence jump (+0.033) is surprisingly small — plateau effect
- Convergence ratio for q=2 is non-monotonic (0.61, 0.36, 0.72) — not pure exponential
- q=3's larger 2D/1D ratio (3.80 vs 3.08) correlates with slower convergence — bigger gap to close
- Exponential fit residuals suggest power-law corrections may be needed
