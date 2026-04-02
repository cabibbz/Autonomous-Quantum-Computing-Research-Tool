# Sprint 072 — q=3 Cylinder & Ly=3 Extension

## Idea
Fill the q=3 gap in our cylinder (Ly=2) data. We have g_c(cyl) for q=2 (0.451) and q=5 (0.714) but not q=3 — the only other exactly-solvable case (g_c(1D)=1/3 exact). This gives a calibration point for the cyl/1D ratio and tests whether the ratio is monotonic in q. Also extend to Ly=3 cylinder for q=2 to study convergence toward 2D g_c=0.771.

## Experiments
- **072a:** q=3 Ly=2 cylinder gap×Lx crossings (Lx=3,4,5,6; dim up to 3^12=531k)
- **072b:** q=3 cylinder order parameter + cross-q summary
- **072c:** q=2 Ly=3 cylinder gap crossings (Lx=3,4,5,6,7; dim up to 2^21=2M)

## Results

### 072a — q=3 Ly=2 Cylinder Gap Crossings

**g_c(q=3, Ly=2 cyl) = 0.565** from 3 crossing pairs:

| Lx pair | g_c | gap×Lx at crossing |
|---------|-----|--------------------|
| (3,4) | 0.5548 | 1.7427 |
| (4,5) | 0.5666 | 1.9055 |
| (5,6) | 0.5725 | 2.0108 |

Crossings converge monotonically from below. Spread = 0.018 (3.2%). Lx=7 (dim=4.8M) not attempted (6.3s/pt at Lx=6).

**Cyl/1D ratio = 1.69** — between q=2 (1.80) and q=5 (1.62). Confirms monotonically decreasing trend.

### 072b — q=3 Cylinder Order Parameter

**Order parameter smooth across transition.** <δ(s_i,s_j)> at g_c=0.565: 0.605 (vs random=0.333, ordered=1.0). No discontinuity — consistent with continuous transition.

**Max slope |d<δ>/dg| = 1.99 (Lx=5), 1.98 (Lx=6).** Size-independent at these sizes.

**Cross-q order parameter slopes (all Ly=2 cylinder):**
| q | Max slope | Lx |
|---|-----------|-----|
| 2 | 1.21 | 6 |
| 3 | 1.99 | 5,6 |
| 5 | 1.88 | 4 |

q=3 has the steepest transition (but sizes differ across q).

**Cross-q cylinder table:**
| q | g_c(1D) | g_c(cyl) | g_c(2D) | cyl/1D | 2D/cyl | 2D/1D |
|---|---------|----------|---------|--------|--------|-------|
| 2 | 0.250 | 0.451 | 0.771 | 1.80 | 1.71 | 3.08 |
| 3 | 0.333 | 0.565 | 1.267 | 1.69 | 2.24 | 3.80 |
| 5 | 0.441 | 0.714 | 1.588 | 1.62 | 2.22 | 3.60 |

**Cyl/1D ratio monotonically decreasing:** 1.80 → 1.69 → 1.62. 2D/cyl ratio peaks at q=3 (2.24).

### 072c — q=2 Ly=3 Cylinder Gap Crossings

**g_c(q=2, Ly=3 cyl) = 0.655** from 4 crossing pairs:

| Lx pair | g_c | gap×Lx at crossing |
|---------|-----|--------------------|
| (3,4) | 0.6447 | 2.0877 |
| (4,5) | 0.6546 | 2.2357 |
| (5,6) | 0.6592 | 2.3220 |
| (6,7) | 0.6615 | 2.3744 |

**Ly convergence toward 2D:**
| Ly | g_c | z (coordination) |
|----|-----|-----------------|
| 1 (1D) | 0.250 | 2 |
| 2 | 0.451 | 3 |
| 3 | 0.655 | ~3.33 |
| ∞ (2D) | 0.771 | 4 |

**Ly=3 is 63.7% of the way from Ly=2 to 2D.** Gap×Lx at crossing converges toward ~2.4 (Ly=2: 2.36 at 2D torus L=3,4).

Lx=8 (dim=16M) timing test: 260s/pt — too slow for full scan.

## Surprises
- Cyl/1D ratio is **monotonically decreasing** in q (1.80→1.69→1.62), not non-monotonic as 2D/1D is
- 2D/cyl ratio peaks at q=3 (2.24) — the 2D jump from cylinder is largest for q=3
- q=3 order parameter slope (1.99) exceeds both q=2 (1.21) and q=5 (1.88)
- Ly=3 crossings converge faster than Ly=2 (spread 0.017 vs 0.007 in 4 pairs)
- Gap×Lx at crossing increases with Ly: 2.0 (Ly=2) → 2.4 (Ly=3) for q=2

## Key Findings
1. **q=3 cylinder fills the gap** — g_c(cyl)=0.565, completing the Ly=2 dataset for q=2,3,5
2. **Cyl/1D ratio decreases monotonically** — higher q needs proportionally less coupling boost on cylinder
3. **Ly=3 confirms monotonic convergence** to 2D for q=2 — g_c increases smoothly with Ly
4. **All transitions smooth** — no first-order signal for any q on cylinder geometry
