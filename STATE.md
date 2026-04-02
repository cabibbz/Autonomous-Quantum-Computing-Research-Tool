# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 084 — Entanglement spectrum at walking boundary. Computed full entanglement spectrum for q=2-8 at g_c=1/q. Walking breakdown = entropy concentration in (q-1)-fold degenerate first excited entanglement level. Level 1 absorbs 48% (q=2) → 69% (q=8) of entropy. Entanglement gap INCREASES with q. Tail entropy correlates with c_eff/Re(c) (r=0.80).

## Active Research Thread
**S_q Potts walking regime: microscopic mechanism identified.**

Walking breakdown is entropy concentration, not a spectral phenomenon:

| q | Level 1 %S | Tail %S | c_eff/Re(c) | Δξ |
|---|-----------|---------|-------------|-----|
| 2 | 47.8% | 19.1% | 1.00 | 1.03 |
| 5 | 62.7% | 15.2% | 1.01 | 1.69 |
| 7 | 66.8% | 13.5% | 0.78 | 1.98 |
| 8 | 69.1% | 11.7% | 0.74 | 2.13 |

The (q-1)-fold degenerate first excited entanglement level absorbs progressively more entropy weight. Energy/gap/correlator observables see only the lowest levels (perfectly conformal). Entropy sums over all levels → sensitive to redistribution.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU, 59 sprints since last use. Encode q=2 S_q Potts at n=6-8 on qubits, measure entanglement spectrum or gap×N.
2. **Entanglement spectrum scaling with n** — Does the level 1 concentration GROW with n at q=7? If so, this directly explains why c_eff drifts downward. Use DMRG for n=10-16.
3. **Rényi entropies at walking boundary** — S_α = (1/(1-α)) ln(Σ λ_i^α). Different α weight the spectrum differently. α→∞ probes λ_max only (should track Re(c)). α→0 counts effective dimension. Map S_α(q) across walking boundary.

## What's Been Ruled Out
- Oscillatory correlator corrections at r ≤ 7 (Im(x_σ) too small)
- Open-BC raw power law for x extraction (inflates η by 5×)
- x_σ as walking discriminator (nearly constant across q)
- Casimir energy as walking discriminator (follows Re(c) for all q)
- Entanglement gap as walking discriminator (INCREASES with q)
- q=6 as second walking case (c_eff drops 2.9% n=8→12)
- Sharp walking-to-first-order transition (smooth crossover)

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=5 (fast, n≤24), q=6 (slow, n≤12), q=7 (n≤12)
- Periodic exact diag correlator: q=2 n≤14, q=3 n≤10, q=5 n≤8, q=7 n≤7
- Entanglement spectrum: same size limits as exact diag
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
