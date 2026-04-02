# Sprint 076 — S_q Potts vs Hybrid: Direct Universality Class Comparison at q=5

**Date:** 2026-04-02
**Status:** Complete (3 experiments)

## Motivation

We've established that the Potts-clock hybrid model has continuous transitions for q≥5 (Sprints 065-066), while the literature predicts S_q Potts is first-order for q>4 (Gorbenko, Rychkov & Zan, JHEP 2018). However, we've never directly built and compared the S_q Potts model. This sprint provides head-to-head evidence that the hybrid defines a genuinely new universality class by showing qualitatively different gap behavior at the same q.

**Key difference in Hamiltonians:**
- **Hybrid (our model):** H = -Σ δ(s_i,s_j) - g(X + X†) — Z_q symmetric field
- **S_q Potts:** H = -Σ δ(s_i,s_j) - g Σ_{k=1}^{q-1} X^k — full S_q symmetric field

For q=2,3 these are equivalent (X† = X for q=2; X² = X† for q=3). For q≥4 they differ.

## Experiments

### 076a — S_q Potts gap crossings at q=5

**Verification:** At q=3, S_q Potts ≡ hybrid (gaps match to 10 decimal places). Correct.

**S_q Potts field:** Σ_{k=1}^{q-1} X^k = (all-ones matrix) - I. Strength (q-1)=4 per unit g, vs hybrid's 2cos(2π/5)≈0.618. Ratio 6.47x.

**Results (periodic BC, n=4,6,8):**

| Model | (4,6) g_c | (6,8) g_c | gap×N at crossing |
|-------|-----------|-----------|-------------------|
| S_q Potts | 0.2016 | 0.2007 | 0.699 → 0.673 |
| Hybrid | 0.4387 | 0.4379 | 0.464 → 0.459 |

**S_q Potts g_c ≈ 0.200** — crossings converge tightly (-0.4% drift). **Hybrid g_c ≈ 0.438** — also tight convergence (matches known 0.441 within 0.7%).

**Surprise: S_q Potts ALSO shows clean gap crossings at n≤8.** If first-order, the crossing should eventually fail at large n (gap closes exponentially, not as 1/N). But at n≤8, the correlation length ξ may be larger than N, masking the first-order nature. The gap×N at crossing DECREASES for S_q (0.699→0.673, -3.7%) but barely for hybrid (0.464→0.459, -1.1%). This faster drift may signal first-order.

### 076b — Conformal tower comparison at q=5, n=4,6,8

**Qualitatively different degeneracy structures:**

| Level | S_q Potts | Hybrid | Interpretation |
|-------|-----------|--------|----------------|
| x₁-x₄ | 4-fold degenerate (x=1.000) | 2+2 split (x=1.000, 2.413) | S_q has S₅ spin degeneracy; hybrid has Z₅ conjugate pairs |
| x₅ | singlet (x≈4.51) | singlet (x≈7.65) | Different scaling dimension |
| x₆-x₉ | 4-fold (x≈8.16) | 3+2+1 split | S_q symmetry combines; Z_q resolves |

**S_q Potts has EXACT 4-fold degeneracy.** x₁=x₂=x₃=x₄=1.0000 at ALL sizes. This is the S₅ permutation symmetry: the q-1=4 degenerate spin field components. The gap to the NEXT level (x₅≈4.51) is huge.

**Hybrid has 2-fold Z₅ conjugate pairs.** x₁=x₂=1.000 (one conjugate pair), x₃=x₄=2.413 (second pair). The Z₅ symmetry resolves what S₅ keeps degenerate.

**x₃ is the KEY discriminator.** S_q: x₃=1.000 (degenerate with x₁). Hybrid: x₃=2.413 (+141% different). This is the spin field splitting that proves different symmetry classes.

**Both models show clean conformal tower convergence.** S_q x₅ converges: 4.81→4.64→4.51. Hybrid x₃ converges: 2.44→2.42→2.41. Both are genuine CFTs at these sizes — the first-order nature of S_q Potts is not visible in the tower at n≤8.

**Surprise: S_q Potts spectrum at n≤8 looks like a CLEAN CFT, not first-order.** Converging tower ratios, exact degeneracies, smooth size dependence. Either (a) ξ_first_order >> 8 so CFT behavior is visible, or (b) the first-order prediction may need revisiting for 1D quantum chain at q=5.

### 076c — Δ·N scaling comparison (first-order vs continuous)

**q=3 control (both models identical):** Δ·N decreases smoothly: 0.746→0.735→0.729→0.725→0.722. Power-law (Δ~N^{-1.03}) fits 52x better than exponential. z_eff=1.03 ≈ 1. Confirms continuous transition.

**q=5 comparison:**

| n | S_q Δ·N | change | Hybrid Δ·N | change |
|---|---------|--------|------------|--------|
| 4 | 0.6713 | — | 0.4607 | — |
| 6 | 0.6507 | -3.07% | 0.4589 | -0.39% |
| 8 | 0.6393 | -1.74% | 0.4586 | -0.07% |

**Hybrid is NEARLY CONVERGED** at n=8 (total drift 0.46%, z_eff=1.007). S_q Potts still drifting (total drift 4.8%, z_eff=1.071). Both fit power-law better than exponential at n≤8.

**S_q Potts FSS corrections are 10x LARGER than hybrid.** The changes per step are -3.07% → -1.74% for S_q vs -0.39% → -0.07% for hybrid. Same direction (decelerating decrease) but vastly different magnitudes.

**q=7 S_q Potts: DRAMATIC gap collapse.** Δ·N drops 69% from n=4 to n=6 (0.175 → 0.055). This is NOT power-law behavior. Either (a) g_c estimate is wrong, or (b) the first-order nature is already visible at q=7. If the estimate g_c≈0.12 is approximately correct, this is a first-order signal.

**Δ₂/Δ₁ = 1.0000 EXACTLY for BOTH models at ALL sizes.** Conjugate pair degeneracy is exact regardless of whether the field is Z_q or S_q symmetric. The first excited state is always a degenerate pair.

## Summary

**First head-to-head comparison of S_q Potts vs Potts-clock hybrid at q=5.** Three independent probes confirm different universality classes:

1. **Different g_c:** S_q Potts g_c=0.200, Hybrid g_c=0.438 (field strength ratio 6.47x accounts for factor ~2.2x in g_c)
2. **Different degeneracy structure:** S_q has 4-fold (S₅ symmetric) first excited state. Hybrid has 2-fold (Z₅ conjugate pairs) with second pair at x₃=2.41
3. **Different FSS corrections:** S_q corrections 10x larger than hybrid. S_q Δ·N drifts 4.8% over n=4-8 vs hybrid's 0.46%

**At n≤8, BOTH models look like CFTs.** Neither shows definitive first-order signatures (no exponential gap closing, no discontinuities). This is consistent with the complex CFT picture (Gorbenko et al.): S_q Potts q=5 first-order transition has ξ >> 8, so CFT-like behavior is visible at accessible sizes.

**q=7 S_q Potts shows potential first-order signal.** The 69% Δ·N drop from n=4→6 is incompatible with power-law. Needs confirmation with correct g_c.

**The clearest difference is the spectrum, not the gap scaling.** The 4-fold vs 2-fold degeneracy pattern unambiguously distinguishes S₅ from Z₅ symmetry. The FSS corrections are quantitatively different but qualitatively similar (both decrease). Only at larger q (q≥7) might gap scaling alone discriminate.

## Surprises

- S_q Potts at q=5 looks like a genuine CFT at n≤8 — no first-order signature visible
- Degeneracy structure is the SHARPEST discriminator (no fit needed, exact integers)
- S_q field = (all-ones) - I, an elegant identity: sum of all cyclic shifts minus identity
- q=7 S_q gap collapse is dramatic but only 2 sizes (inconclusive)
- Both models have Δ₂/Δ₁ = 1.0000 exactly — Z_q conjugate pairs always degenerate even with S_q field
- q=3 verification perfect: 10-digit match confirms code correctness

[Full data: results/sprint_076a_sq_potts_gap_q5.json, results/sprint_076b_tower_comparison.json, results/sprint_076c_gap_scaling_comparison.json]
