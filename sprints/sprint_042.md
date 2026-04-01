# Sprint 042 — True q=5 Potts MI-CV: Second-Order, Not First-Order

**Status:** Complete (5 experiments)

## Motivation

Sprint 041 showed clock model crossings persist at q=5, but ClockChain ≠ Potts for q≥4. The 2D classical correspondence predicts q>4 Potts is first-order. We built a custom Potts model with Kronecker-delta coupling to test: does true Potts q=5 show step function (first-order) or crossings (second-order)?

## Physical difference: Clock vs Potts

- **Clock**: H_bond = cos(2π(s_i-s_j)/q) — graduated penalty for misalignment
- **Potts**: H_bond = δ(s_i,s_j) — binary same/different (all misalignments equal)
- Equivalent for q=2,3. Different universality for q≥4.

## Prediction: WRONG
q=5 Potts should show step function (first-order, no crossings). **It shows crossings instead.**

## Experiment 042a — Timing Test

Custom PottsChain model built: q projectors P_a for Kronecker-delta coupling + shift operator X for transverse field. Runtime ~17s/point at n=8 (same as clock d=5).

First surprise: g=0.8 gives CV=1.49 (vs clock 0.36). Potts transition is at much lower g.

## Experiment 042b,c — n=8 Full Profile

| g/J | Potts CV | Clock CV |
|-----|----------|----------|
| 0.10 | 0.012 | — |
| 0.30 | 0.157 | — |
| 0.35 | 0.295 | — |
| 0.40 | 0.547 | — |
| 0.45 | 0.839 | — |
| 0.50 | 1.065 | 0.141 |
| 0.80 | 1.492 | 0.362 |
| 1.00 | — | 0.521 |

Potts transition region: g=0.30-0.50 (CV rises 0.16→1.06). Much steeper than clock.

## Experiment 042d — n=12 Crossing Test (KEY RESULT)

| g/J | n=8 CV | n=12 CV | Change |
|-----|--------|---------|--------|
| 0.30 | 0.157 | 0.134 | -0.023 ↓ (ordered) |
| 0.40 | 0.547 | 0.516 | -0.031 ↓ (ordered) |
| 0.50 | 1.065 | 1.360 | +0.296 ↑ (disordered) |

**CROSSING at g_c ≈ 0.41** (interpolated). This is the second-order signature (CV decreases in ordered phase, increases in disordered).

**First-order prediction (step function, no crossing) is WRONG.**

## Experiment 042e — Potts vs Clock Comparison

| Property | Potts q=5 | Clock q=5 |
|----------|-----------|-----------|
| Crossing g_c | 0.41 | 0.67 |
| Transition type | Second-order (crossing) | Second-order (crossing) |
| n=8 slope | 4.54 | 0.80 |
| n=12 slope | 6.13 | 1.01 |
| Slope ratio n12/n8 | 1.35 | 1.27 |

Potts transition is **5.7x steeper** than clock but at a much lower g value.

## Key Findings

1. **True q=5 Potts shows CROSSING curves** — not first-order step function
2. **Potts crossing at g_c ≈ 0.41**, much lower than clock's 0.67
3. **Potts slope 5.7x steeper than clock** — stronger, more Ising-like transition
4. **Both models show second-order MI-CV signature** despite different universality
5. **The 2D classical "q>4 → first-order" does NOT apply to 1D quantum Potts with transverse field**

## Surprises

- Main prediction WRONG: Potts q=5 is NOT first-order in 1D quantum
- Potts transition 5.7x steeper than clock — not just shifted, fundamentally different character
- Potts ordered phase much more Democratic (CV=0.012 at g=0.1 vs clock 0.14 at g=0.5)
- Slope growth ratio similar (Potts 1.35 vs Clock 1.27) despite very different absolute slopes

## Interpretation

The 2D classical Potts model has a first-order transition for q>4. But the quantum-classical correspondence maps 1D quantum to 2D classical with ANISOTROPIC couplings. The transverse field introduces a different universality from the isotropic 2D model. In the anisotropic limit (weak transverse field), the 1D quantum Potts transition remains second-order because the correlation length in the "imaginary time" direction diverges differently.

The Potts-Clock distinction manifests as:
- Different transition location (g_c = 0.41 vs 0.67)
- Very different slopes (4.5 vs 0.8)
- Different ordered-phase CV (0.01 vs 0.14)
- But the SAME qualitative signature (crossings = second-order)

MI-CV correctly identifies both as second-order but reveals quantitative differences consistent with different universality classes.
