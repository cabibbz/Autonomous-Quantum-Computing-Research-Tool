# Sprint 041 — q=5 Clock MI-CV: Crossings Persist, But Shifted Dramatically

**Status:** Complete (4 experiments)

## Motivation

q=5 is definitively above the marginal q=4 boundary where 2D Potts transitions change from second-order to first-order. The prediction: q=5 MI-CV should show NO crossings (either step function or BKT dome).

**Important caveat:** We're using TeNPy's `ClockChain`, which is the q-state clock model, NOT the q-state Potts model. For q=2 (Ising) and q=3, clock ≡ Potts. But for q≥4, the clock model differs — it has nearest-neighbor Z_q coupling cos(2π(s_i-s_j)/q) rather than Potts δ(s_i,s_j). For q≥5, the clock model is expected to have an intermediate "floating" BKT phase.

## Predictions (ALL WRONG)
1. If clock model: expect BKT dome (no crossings) — **WRONG**
2. If Potts-like: expect step function (first-order, no crossings) — **WRONG**
3. Either way: NO crossings expected (unlike q≤4) — **WRONG**

## Experiment 041a — q=5 n=8 Timing Test

- SU(5) has 24 generators (d²-1), 576 correlation pairs per site pair
- ~16.5s/point at n=8 with chi=15 (12.8s in correlation computation)
- Compare: q=4 was ~6.5s/point — q=5 is 2.5x slower as expected

## Experiment 041b — q=5 n=8 Full Sweep

| g/J | q=5 n=8 CV | q=4 n=8 CV |
|-----|-----------|-----------|
| 0.50 | 0.1409 | 0.088 |
| 0.80 | 0.3616 | 0.372 |
| 0.90 | 0.4384 | 0.547 |
| 0.95 | 0.4790 | 0.639 |
| 1.00 | 0.5213 | 0.727 |
| 1.05 | 0.5651 | 0.810 |
| 1.10 | 0.6103 | 0.887 |
| 1.30 | 0.7941 | 1.118 |

q=5 CV is systematically lower than q=4 above g≈0.85. Slope at g=1.0 is HALF of q=4 (0.86 vs 1.72).

## Experiment 041c — q=5 n=12 (Key Points)

| g/J | n=8 CV | n=12 CV | Change |
|-----|--------|---------|--------|
| 0.50 | 0.1409 | 0.1206 | -0.020 ↓ |
| 0.80 | 0.3616 | 0.3766 | +0.015 ↑ |
| 1.00 | 0.5213 | 0.5794 | +0.058 ↑ |
| 1.10 | 0.6103 | 0.7017 | +0.091 ↑ |

**SURPRISE: Crossing exists at g ≈ 0.67.** CV decreases with n at g=0.5 (ordered) but increases at g≥0.8. Interpolated crossing at g≈0.673.

## Experiment 041d — Cross-q Analysis

### Crossing point trend with q (n=8,12)

| q | Model | g_c (crossing) | Shift from q-1 |
|---|-------|----------------|-----------------|
| 2 | Ising (=clock) | 0.930 | — |
| 3 | Potts (=clock) | 0.923 | -0.007 |
| 4 | Clock | 0.893 | -0.030 |
| 5 | Clock | 0.673 | **-0.220** |

The q=5 crossing shifts by **0.22** — an order of magnitude larger jump than q=3→4 (0.030). The trend is dramatically non-linear.

### Slope at g=1.0 (n=8)

| q | Slope |
|---|-------|
| 4 | 1.72 |
| 5 | 0.86 |

q=5 slope is exactly half of q=4.

### CV at g=1.0

| q | n=8 | n=12 | Growth |
|---|-----|------|--------|
| 4 | 0.727 | 0.838 | +15% |
| 5 | 0.521 | 0.579 | +11% |

q=5 grows slower — weaker size dependence consistent with BKT crossover.

## Key Findings

1. **q=5 clock model STILL shows crossing curves** — prediction that crossings would disappear was wrong
2. **Crossing shifts dramatically to g≈0.67** — far below self-dual, 10x larger shift than q=3→4
3. **Slope halves from q=4 to q=5** — transition becomes gentler, consistent with BKT character
4. **CV systematically lower at q=5** — larger local dimension distributes correlations more evenly
5. **Clock ≠ Potts distinction matters:** The q≥4 clock model has different universality from Potts, explaining why crossings persist (second-order/BKT transition rather than first-order)

## Surprises
- ALL three predictions were wrong — crossings persist at q=5
- The g_c shift is 10x larger from q=4→5 than q=3→4 (0.22 vs 0.03)
- Slope exactly halves — too clean to be coincidence?
- The crossing is deep in what would be the "ordered phase" for q≤4

## Interpretation

The ClockChain is NOT the Potts model for q≥4. The 1D quantum q-state clock model for q=5 likely has a BKT-like transition (or two BKT transitions with an intermediate floating phase). At n=8,12, the crossings suggest the ordered→floating transition still looks second-order-like from MI-CV's perspective, but with dramatically reduced slope and shifted crossing point. Testing q≥5 Potts (not clock) would require a custom model with Kronecker-delta coupling.
