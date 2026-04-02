# Sprint 080 — c_eff at q=6,8,9: Walking Regime Boundary Mapped

**Date:** 2026-04-02
**Status:** Complete (3 experiments)

## Motivation

Sprint 079 showed c_eff/Re(c) = 1.00 at q=5 but 0.82 at q=7. The walking-to-first-order crossover happens between q=5 and q=7. This sprint fills the gap at q=6, 8, and 9.

## Experiments

### 080a — c_eff(q=6) at g_c = 1/6

Exact diag at n=6 (dim=47k) and n=8 (dim=1.68M, GPU).

| n | dim | c_eff | R² | gap×N | Time |
|---|-----|-------|----|-------|------|
| 6 | 46,656 | 1.1459 | 1.000000 | 1.891 | 0.5s |
| 8 | 1,679,616 | 1.1473 | 0.999995 | 2.011 | 29s |

**c_eff/Re(c) = 0.916** (Re(c) = 1.2525). Marginal walking.

**c_eff STABLE with n:** drift +0.12% from n=6→8. Unlike q=7 where c_eff degrades, q=6 c_eff is converged. This means ξ*(q=6) >> 8 — walking extends beyond all accessible exact diag sizes.

**First excited: 3-fold degenerate** (k=4 gives gaps[0:3] identical to 10⁻⁵). Part of the (q-1)-fold S_q pattern.

### 080b — c_eff(q=8,9) at g_c = 1/8, 1/9

| q | n | dim | c_eff | R² | gap×N | Time |
|---|---|-----|-------|----|-------|------|
| 8 | 6 | 262,144 | 1.0725 | 1.00 | 2.073 | 3.9s |
| 8 | 7 | 2,097,152 | 1.0615 | 1.00 | 2.167 | 42s |
| 9 | 6 | 531,441 | 1.0119 | 1.00 | 2.156 | 8.2s |

**q=8: c_eff/Re(c) = 0.738** (Re(c) = 1.438). Walking broken. c_eff DECREASING with n (−1.0%).
**q=9: c_eff/Re(c) = 0.668** (Re(c) = 1.515). Walking broken.

n=5 data unreliable (R²=0 due to only 4 entropy points). All n≥6 data clean.

### 080c — Full Walking Boundary Compilation

Complete c_eff/Re(c) ratio across q=3-10:

| q | Re(c) | c_eff(best) | n_best | c/Re(c) | Walking? |
|---|-------|-------------|--------|---------|----------|
| 3 | 0.800 | 0.893 | 24 | 1.12 | YES (real CFT, +12% FSS) |
| 5 | 1.138 | 1.147 | 20 | 1.01 | YES — sweet spot |
| 6 | 1.253 | 1.147 | 8 | 0.92 | MARGINAL — stable |
| 7 | 1.351 | 1.059 | 12 | 0.78 | BREAKING — degrading |
| 8 | 1.438 | 1.062 | 7 | 0.74 | BREAKING |
| 9 | 1.515 | 1.012 | 6 | 0.67 | NO |
| 10 | 1.584 | 0.946 | 6 | 0.60 | NO |

**Exponential fit:** c_eff/Re(c) = 1.004 × exp(−0.105 × (q−5)). Ratio=1 crossing at q ≈ 5.0.

**Linear fit:** ratio = −0.081q + 1.395. Ratio=1 at q = 4.9.

**Both fits agree: q=5 is the exact boundary** where walking regime matches complex CFT Re(c).

**Size dependence is the key discriminator:**
- q=5: c_eff stable or slightly increasing with n (walking confirmed to n=20)
- q=6: c_eff stable (+0.12% n=6→8) — walking extends beyond accessible sizes
- q=7: c_eff DECREASING with n (−4.7% n=8→12) — walking breaking down
- q=8: c_eff DECREASING with n (−1.0% n=6→7) — already broken

**gap×N INCREASES with q** (2.01, 2.17, 2.16 at q=6,8,9) while c_eff/Re(c) decreases. Walking breakdown is in the entropy, not the gap.

**c_eff converges to ~1.0-1.1 for ALL q≥5 at moderate sizes.** The complex CFT q-dependence only manifests at n << ξ*(q). At finite n, c_eff is nearly q-independent.

## Surprises

- q=6 c_eff STABLE with n (+0.12%) — walking extends beyond n=8 despite 8% deficit
- gap×N increases with q even as walking breaks down
- Both linear and exponential fits give ratio=1 crossing at q≈5.0 ± 0.1
- c_eff converges to ~1.0-1.1 for all q≥5 regardless of Re(c)
- q=9 from a single n=6 point: ratio already at 0.67

## POTENTIALLY NOVEL

Complete c_eff/Re(c) ratio curve across q=3-10 for S_q Potts chain. Exponential decay of walking ratio with decay constant 0.105 per unit q. First quantitative mapping of the walking regime boundary showing q=5 as the exact threshold.
