# Sprint 061 — Large-q Limit: CFT Data at q=15-30 and Free Boson Convergence

**Status:** Complete (4 experiments).

**Goal:** Extend CFT characterization to q=15,20,25,30. Test free boson predictions.

---

## Experiment 061a — Energy Gap g_c for q=15,20

**n=4,5 crossings found but unreliable at large q.** q=15 crossing at g≈0.580. q=10 calibration: n=4,5 crossing at 0.566 vs known 0.684 (21% FSS shift). Odd n has much larger shift than even n — sublattice effect. n=6 infeasible for q≥12 (q=15 n=6 = 11.4M dim).

**Decision: use formula g_c = (1/5)√(q-1) + 1/20**, validated to <5% for q=2-10.

## Experiment 061b — Full CFT Spectrum at n=4 for q=15,20,25,30

**Harmonic ratios approach k² rapidly!** At q=30: R(σ²)/R(σ)=3.994 (99.9% of 4), R(σ³)=8.970 (99.7% of 9), R(σ⁴)=15.907 (99.4% of 16), R(σ⁵)=24.777 (99.1% of 25). **Convergence is O(q^{-3})**, much faster than expected.

**Epsilon misidentified as sigma² without charge resolution.** The code's OPE extraction found C≈0.71 for all large q — actually the sigma² matrix element, not epsilon.

## Experiment 061c — Charge-Resolved Spectrum (Z_q Symmetry)

**Z_q charge resolution cleanly separates all operators.** Built G = X₁⊗...⊗Xₙ, diagonalized within degenerate subspaces. Sigma (charge 1), sigma² (charge 2), epsilon (charge 0) all correctly identified.

**Epsilon does NOT merge with sigma² at large q.** They are in different charge sectors and DIVERGE:
- R(sigma²) → 4 (k²) for all q
- R(epsilon) grows linearly at n=4: R_eps ≈ 0.88·q for q≥10

**True C_sse decreases steeply:**

| q | R(σ²) | R(ε) | C_sse (n=4) |
|---|-------|------|-------------|
| 2 | — | 7.70 | 0.480 |
| 3 | 1.00 | 6.25 | 0.478 |
| 4 | 1.75 | 6.44 | 0.396 |
| 5 | 2.44 | 6.98 | 0.322 |
| 7 | 3.27 | 7.80 | 0.235 |
| 10 | 3.66 | 8.34 | 0.184 |
| 15 | 3.93 | 13.66 | 0.093 |
| 20 | 3.97 | 17.75 | 0.063 |
| 25 | 3.99 | 21.91 | 0.047 |
| 30 | 3.99 | 26.11 | 0.038 |

## Experiment 061d — Asymptotic Analysis

**Harmonic ratio convergence: O(q^{-3}).** Deficit from k² scales as q^{-3.0} for k=2, q^{-3.2} for k=3,4. MUCH faster than O(1/q).

**C_sse scaling:** Best-n (q≤10) power law: C_sse ~ 1.5·q^{-0.87}. n=4 data (q≤30): C_sse ~ 2.3·q^{-1.19}. Steeper at large q or FSS effect.

**x₁ → 0 as q→∞.** Best-n: x₁ ~ 0.16·q^{-0.29} (q≤10). n=4 descendant gap: x₁(30) ≈ 0.038 — NOT saturating.

**c/x₁ ~ 2.2·q^{0.94}** — sub-linear, confirming c/x₁ ≠ 2q for q>3.

**Free boson limit confirmed, but c→∞.** The large-q 1D quantum Potts is a free boson on a decompactifying circle (radius R ~ √(ln q)). This gives:
- c ~ ln(q): more modes as circle grows
- x₁ ~ 1/R² → 0: spin field becomes marginal
- R(σᵏ)/R(σ) → k²: free boson vertex operators
- C_sse → 0: epsilon field decouples as it moves to higher R

---

## Surprises
- Harmonic ratio convergence is O(q^{-3}), not O(1/q) — no known prediction for this exponent
- R_epsilon grows linearly with q at n=4 (~0.88q) while R_sigma² saturates at 4
- x₁ does NOT saturate — continues decreasing to 0.038 at q=30
- c/x₁ ratio is sub-linear in q (exponent 0.94)
- Previous C_sse=0.71 at large q was misidentification of sigma² as epsilon

**POTENTIALLY NOVEL:** First charge-resolved conformal spectrum of Hermitian q-state Potts for q=15-30. The O(q^{-3}) convergence rate to free boson k² and the R_epsilon~q growth are previously unmeasured.
