# Sprint 040 — q=4 Potts MI-CV: Marginal Transition Shows Crossing, Not Dome

**Date:** 2026-04-01
**Status:** Complete (3 experiments: 040a2, 040b/b2, 040c)

## Idea

The 2D classical q-state Potts model has second-order transitions for q≤4 and first-order for q>4. At q=4, the transition is marginal — at the boundary, with logarithmic corrections (BKT-like). The quantum 1D q=4 Potts (clock) model should inherit this marginal character.

**Question:** Does q=4 Potts MI-CV show crossing curves (second-order) or a dome (BKT-like marginal)?

**Prediction:** Crossing curves, but with weaker/slower convergence due to logarithmic corrections.

**Technical:** d=4 per site requires 15 generalized Gell-Mann matrices (SU(4) generators). DMRG with d=4 is ~3x slower than d=3 per point.

## Experiments

### 040a/040a2: q=4 Potts MI-CV at n=8

chi=15 sufficient (validated: chi=15 and chi=40 give identical CV=0.7273 at g=1.0). ~6.5s/point.

| g/J | n=8 CV |
|-----|--------|
| 0.50 | 0.088 |
| 0.80 | 0.372 |
| 0.90 | 0.547 |
| 0.95 | 0.639 |
| 1.00 | 0.727 |
| 1.05 | 0.810 |
| 1.10 | 0.887 |
| 1.30 | 1.118 |

Monotonically increasing with g, same qualitative shape as q=3 Potts.

### 040b/040b2: q=4 Potts MI-CV at n=12

chi=20, ~18s/point. Five key g values:

| g/J | n=8 CV | n=12 CV | Direction |
|-----|--------|---------|-----------|
| 0.80 | 0.372 | 0.323 | ↓ ordered |
| 0.90 | 0.547 | 0.551 | ↑ disordered |
| 0.95 | 0.639 | 0.693 | ↑ disordered |
| 1.00 | 0.727 | 0.838 | ↑ disordered |
| 1.10 | 0.887 | 1.099 | ↑ disordered |

**Crossing confirmed!** CV decreases with n on the ordered side (g=0.8) and increases on the disordered side (g≥0.9). The crossing point is at g_c ≈ 0.893 (linear interpolation between g=0.8 and g=0.9).

### 040c: Cross-model comparison

**Crossing point comparison across q:**
| Model | q | g_c (n=8,12) | Self-dual |
|-------|---|-------------|-----------|
| TFIM | 2 | ~0.93 | h=1.0 |
| Potts | 3 | 0.923 | g=1.0 |
| Potts | 4 | 0.893 | g=1.0 |

Crossing point shifts further below self-dual with increasing q. Larger finite-size effect at larger q.

**Transition slope at g=1.0:**
| n | q=3 slope | q=4 slope |
|---|-----------|-----------|
| 8 | 2.27 | 1.72 |
| 12 | 3.98 | 2.71 |

q=4 slope LOWER than q=3 at same n. Slope ratio n=12/n=8: q=3 gives 1.75, q=4 gives 1.58.

**q=4 CV is systematically lower than q=3 above transition:**
| g/J | q=3 CV | q=4 CV | q4/q3 ratio |
|-----|--------|--------|-------------|
| 0.50 | 0.069 | 0.088 | 1.27 |
| 0.80 | 0.341 | 0.372 | 1.09 |
| 0.90 | 0.577 | 0.547 | 0.95 |
| 0.95 | 0.705 | 0.639 | 0.91 |
| 1.00 | 0.826 | 0.727 | 0.88 |
| 1.05 | 0.933 | 0.810 | 0.87 |

q4/q3 ratio crosses 1.0 near g=0.85, then q=4 has MORE uniform MI above transition. Larger d with more internal degrees of freedom creates more evenly distributed correlations in the disordered phase.

## Key Insights

### 1. q=4 Potts Shows Crossing Curves, Not Dome
Despite being at the marginal boundary (q=4 is where 2D Potts transitions from second-order to first-order), MI-CV shows the same crossing signature as q=3 and TFIM. The marginal/BKT-like character does NOT manifest as a dome at n=8,12.

### 2. Logarithmic Corrections May Be Hidden at Small n
The q=4 transition has ν=2/3 (exact) with multiplicative logarithmic corrections. At n=8,12, these corrections may be subdominant to the leading power-law behavior that produces crossings. Testing at n=16-32 would reveal whether crossings persist or develop into a dome.

### 3. q=4 Is More Uniform Than q=3 Above Transition
Counter-intuitive: larger local dimension (d=4 vs d=3) gives LOWER MI-CV in the disordered phase. More internal degrees of freedom distribute correlations more evenly. The q4/q3 ratio decreases systematically from 1.27 (deep ordered) to 0.87 (disordered).

### 4. Crossing Point Shifts Further From Self-Dual With q
g_c(n=8,12): 0.923 (q=3) → 0.893 (q=4). Finite-size effects grow with q, consistent with larger corrections at the marginal point.

## Surprises
- **No dome at q=4** — The prediction that marginal transitions show dome signature was falsified at these sizes. Crossings dominate.
- **q=4 slope LOWER than q=3** — Expected opposite (smaller ν → steeper). The logarithmic corrections at q=4 suppress the slope at small n.
- **q4/q3 ratio < 1 above transition** — Larger d means MORE uniform MI, not less. The disordered phase "democratizes" better with more states.

## Next Steps
1. q=4 at n=16 — does the crossing persist or start showing logarithmic broadening?
2. q=5 Potts — first-order transition. Should show step function (no crossings), confirming MI-CV classification across the q=4 boundary.
3. Universal scaling functions F(x) for q=2,3,4 on same axes — test whether the collapsed curves are visually distinct.
