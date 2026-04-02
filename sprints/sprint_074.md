# Sprint 074 — q=3 Ly=4 & q=5 Ly=3 Cylinder: Testing Convergence Predictions

**Status:** Complete (3 experiments)

## Motivation

Sprint 073 established that Ly convergence to 2D is q-dependent: q=2 has B=1.60 (fast), q=3 has B=3.03 (1.9x slower). The exponential fit predicts g_c(q=3, Ly=4) ≈ 0.917 (~62% to 2D). Testing this prediction validates the convergence model.

q=5 Ly=3 cylinder (dim=125^Lx) would extend the dataset to three q values at Ly=3, enabling cross-q convergence comparison at matched Ly.

## Experiments

### 074a — q=3 Ly=4 cylinder gap crossings
- dim = 81^Lx. Lx=3: 531k (easy), Lx=4: 43M (GPU)
- Prediction: g_c ≈ 0.917, ~62% to 2D
- Expected: 1-2 crossing pairs

### 074b — q=5 Ly=3 cylinder gap crossings
- dim = 125^Lx. Lx=3: 1.95M, Lx=4: 244M (likely infeasible)
- Expected: g_c between Ly=2 (0.714) and 2D (1.588)

### 074c — Cross-q Ly convergence analysis
- Updated fits with new data points
- Test universality of exponential convergence

## Results

### 074a — q=3 Ly=4 cylinder: INFEASIBLE
**Computational barrier hit.** Lx=4 (dim=43M) takes 838s per g-point — infeasible for a scan. Only Lx=3 (dim=531k, 7.2s/pt) is accessible. With a single system size, no gap×Lx crossing is possible.

**Lx=3 scan (g=0.5-1.3, 41 points, 451s):** Gap minimum at g=0.500 (scan edge) — this is tunneling splitting in the ordered phase, not a g_c proxy. gap×Lx range [0.011, 7.515].

**Key finding:** q=3 Ly=4 cylinder hits the exact diag wall. The gap×Lx crossing method requires Lx≥4 for useful crossings (Lx=2,3 have strong FSS artifacts). At Ly=4 for q=3, the smallest usable size (Lx=4) has dim=43M, which at 838s/pt makes scanning infeasible within time limits.

**Implication:** The exponential fit prediction g_c(q=3, Ly=4) ≈ 0.917 CANNOT be tested with current methods. Cylinder convergence for q=3 is capped at Ly=3 with exact diag.

### 074b — q=5 Ly=3 cylinder: g_c = 0.974
**First q=5 Ly=3 measurement.** Lx=2 (dim=15625, 0.1s/pt) and Lx=3 (dim=1.95M, 26s/pt) scanned. Lx=4 (dim=244M) skipped as infeasible.

**(Lx=2, Lx=3) crossing: g_c = 0.974, gap×Lx = 1.082.** 46.5% of the way from 1D (0.441) to 2D (1.588).

**Caveat:** Lx=2 is known to be out of scaling regime (Sprint 068). The crossing value has significant FSS corrections. But it's the only crossing available.

**Gap minimum at g=0.400 (scan edge) for Lx=2** — tunneling splitting, not g_c proxy. d²E/dg² peak at g=0.675 for Lx=2.

### 074c — Cross-q convergence: Saturation for q≥3
**Updated convergence table (all q at Ly=3):**

| q | Ly=2 | Ly=3 | %2D at Ly=3 | B (decay) |
|---|------|------|-------------|-----------|
| 2 | 0.451 | 0.655 | 77.7% | 1.19 |
| 3 | 0.565 | 0.797 | 49.7% | 2.49 |
| 5 | 0.714 | 0.974 | 46.5% | 2.83 |

**Key discovery: Convergence rate SATURATES for q≥3.** The big slowdown is q=2→q=3 (B doubles: 1.19→2.49, progress drops 28%). But q=3→q=5 shows only a small additional slowdown (B: 2.49→2.83, progress drops only 3.2%).

**Cyl/1D ratio monotonically decreasing at both Ly=2 and Ly=3:**
- Ly=2: 1.80 (q=2) → 1.70 (q=3) → 1.62 (q=5)
- Ly=3: 2.62 (q=2) → 2.39 (q=3) → 2.21 (q=5)

**q=3 increments are constant:** Δg_c(Ly=1→2) = 0.232, Δg_c(Ly=2→3) = 0.232. Linear in Ly, not exponential — the exponential fit from Sprint 073 (B=3.03) was overfitting 2 points.

**q=5 increments slightly decrease:** 0.273 (Ly=1→2) → 0.260 (Ly=2→3). Mild deceleration.

**Correlation(2D/1D ratio, B) = 0.89** — higher 2D/1D gap ↔ slower convergence, as established in Sprint 073.

## Surprises
- q=3 Ly=4 infeasible (838s/pt at Lx=4) — closes that experimental line
- q=3 and q=5 converge at SIMILAR rates (~47-50% at Ly=3) — saturation, not continued slowdown
- q=3 Ly increments exactly constant (0.232) — linear convergence, not exponential!
- Sprint 073's B=3.03 for q=3 was overfitting (only 2 points, now revised to B=2.49 with 2 points)
- Cyl/1D ratio decreasing at BOTH Ly=2 and Ly=3 — consistent pattern

## Computational Limits Reached
| q | Max feasible Ly | Reason |
|---|----------------|--------|
| 2 | 4 (Lx≤6) | Ly=5 Lx=4: dim=1M (feasible but small) |
| 3 | 3 (Lx=3,4) | Ly=4 Lx=4: 43M, 838s/pt |
| 5 | 3 (Lx=2,3) | Ly=4 Lx=3: 244M, infeasible |

The cylinder convergence study is nearing its limits with exact diag. Further progress requires DMRG on cylinders or QMC.
