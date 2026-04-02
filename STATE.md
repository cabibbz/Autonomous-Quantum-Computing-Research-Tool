# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 082 — Spin-spin correlator in walking regime. x_σ(q) ≈ 0.13 nearly universal (q=2-8). Conformal form exact to 0.04%. No oscillatory corrections from complex exponents. Walking breakdown is an ENTROPY phenomenon — invisible in correlators and x_σ. First velocity extraction: v(q) decreases from 1.0 (q=2) to 0.66 (q=8).

## Active Research Thread
**S_q Potts walking regime: correlator analysis complete.**

x_σ from periodic chain correlators (Sprint 082c):

| q | x_σ | v(q) | gap×N | Walking? |
|---|-----|------|-------|----------|
| 2 | 0.123 | 1.02 | 0.786 | real CFT |
| 3 | 0.132 | 0.89 | 0.733 | real CFT |
| 5 | 0.136 | 0.75 | 0.639 | YES |
| 7 | 0.132 | 0.68 | 0.560 | BREAKING |
| 8 | 0.130 | 0.66 | 0.536 | NO |

Key insight: Walking vs non-walking is invisible in x_σ and correlator form at accessible sizes. Only c_eff (entropy) and velocity distinguish them. This means the gap-entropy decoupling (Sprints 079-081) extends to correlators too.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU, 57 sprints since last use. Best candidate: encode q=2 S_q Potts at n=6-8 on qubits, measure gap×N and/or correlator. Compare to exact diag predictions.
2. **Velocity v(q) from Casimir energy** — Independent v extraction from E₀/N at two sizes. Cross-validate Sprint 082 v(q) curve. If they match, v(q) is established.
3. **Correlator at larger sizes (DMRG periodic)** — DMRG with periodic BC (expensive but clean) for q=5 at n=16-24. Need 5+ data points to detect oscillatory corrections. Requires DMRG-X or period-BC DMRG.

## What's Been Ruled Out
- Oscillatory correlator corrections at r ≤ 7 (Im(x_σ) too small)
- Open-BC raw power law for x extraction (inflates η by 5×)
- x_σ as walking-regime discriminator (nearly constant across q)
- q=6 as second walking case (c_eff drops 2.9% n=8→12)
- Sharp walking-to-first-order transition (smooth crossover)

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=5 (fast, n≤24), q=6 (slow, n≤12), q=7 (n≤12)
- Periodic exact diag correlator: q=2 n≤14, q=3 n≤10, q=5 n≤8, q=7 n≤7
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
