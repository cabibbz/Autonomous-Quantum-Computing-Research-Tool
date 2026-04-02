# Sprint 073 — Cylinder Ly Convergence: q=3 Ly=3 & q=2 Ly=4

**Date:** 2026-04-02
**Status:** In progress

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

**Results:** [pending]

### 073c — Cross-q Ly Convergence Analysis

Compare convergence rates for q=2 and q=3 as function of Ly.
Does q=3 approach 2D faster or slower?

**Results:** [pending]
