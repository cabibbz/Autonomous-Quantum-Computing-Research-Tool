# Sprint 042 — True q=5 Potts MI-CV: First-Order Transition Test

**Status:** In progress

## Motivation

Sprint 041 showed clock model crossings persist at q=5. But ClockChain uses cos(2π(s_i-s_j)/q) coupling, NOT Potts Kronecker-delta δ(s_i,s_j). For q≥4 these are different models with different universality.

The 1D quantum q-state Potts model (with true Kronecker-delta coupling) should undergo a first-order transition for q>4 (by 2D classical correspondence). Our MI-CV classification (Sprint 037) predicts first-order → step function, no crossings.

**Key test:** Clock q=5 shows crossings (Sprint 041). True Potts q=5 should show NO crossings (step function). If confirmed, MI-CV cleanly separates clock from Potts universality.

## Physical difference: Clock vs Potts

Clock bond energy: cos(2π(a-b)/q) — graduated penalties for misalignment
Potts bond energy: δ(a,b) — binary same/different

For q=2,3: equivalent (up to additive/multiplicative constants)
For q≥4: genuinely different models, different universality classes

## Predictions
1. q=5 Potts MI-CV shows step function (no crossings) — first-order transition
2. Ordered-phase CV ≈ 0 (like FM phase in XXZ first-order, Sprint 037)
3. Transition steeper than clock model but WITHOUT crossing curves

## Experiment 042a — Custom Potts Model Timing Test
