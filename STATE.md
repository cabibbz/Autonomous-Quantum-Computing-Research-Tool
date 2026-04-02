# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 083 — Casimir energy cross-validates complex CFT. c_implied from E₀ matches Re(c) to ±3% for ALL q=2-8, even where entropy c_eff deviates 40%. Walking breakdown is EXCLUSIVELY an entropy phenomenon. v(q) from Casimir agrees with gap/correlator to <3%.

## Active Research Thread
**S_q Potts walking regime: energy vs entropy hierarchy established.**

Three observables probe different aspects of complex CFT:

| Observable | Formula | Sees Re(c)? | Walking signature? |
|------------|---------|-------------|-------------------|
| Casimir E₀ | ε_∞ - πvc/(6N²) | YES (all q) | NO |
| Correlator G(r) | chord^{-2x_σ} | N/A (x_σ) | NO |
| Entropy S | (c/6)·ln(N) | q≤5 only | YES (q>5 deviates) |

The reduced density matrix is uniquely sensitive to walking breakdown. Ground state energy, gap, and correlators all behave as exact CFT.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU, 58 sprints since last use. Best candidate: encode q=2 S_q Potts at n=6-8 on qubits, measure gap×N and/or correlator. Compare to exact diag predictions.
2. **Entanglement spectrum at walking boundary** — The entropy fails but what about the full entanglement spectrum? Compare ρ_A eigenvalues for q=5 (walking) vs q=7 (broken). Does the spectrum shape change, or just the overall scale?
3. **Casimir higher-order corrections** — Fit E₀/N to include 1/N⁴ term. Does the correction coefficient show walking signature where 1/N² doesn't?

## What's Been Ruled Out
- Oscillatory correlator corrections at r ≤ 7 (Im(x_σ) too small)
- Open-BC raw power law for x extraction (inflates η by 5×)
- x_σ as walking-regime discriminator (nearly constant across q)
- Casimir energy as walking discriminator (follows Re(c) for all q)
- q=6 as second walking case (c_eff drops 2.9% n=8→12)
- Sharp walking-to-first-order transition (smooth crossover)

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=5 (fast, n≤24), q=6 (slow, n≤12), q=7 (n≤12)
- Periodic exact diag correlator: q=2 n≤14, q=3 n≤10, q=5 n≤8, q=7 n≤7
- Casimir energy: periodic chain, same size limits as above
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
