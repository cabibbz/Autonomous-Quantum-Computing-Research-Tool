# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 090 — q=4 entanglement spectrum confirms M/[(q-1)/q] crossover at q≈4. Three experiments: (1) DMRG q=4 n=8-24 gives M/[(q-1)/q] = 1.003 at n=16, crossing 1.0 between n=16-20. (2) Tail exponent b=2.024 — universal b≈2.0 confirmed for q=2,3,4,5,7 (mean 2.018 ± 0.029). (3) q_cross converges from 4.49 (n=8) to 3.92 (n=24) — approaches q=4 in thermodynamic limit.

## Active Research Thread
**Entanglement spectrum as walking/CFT probe — democracy index crossover confirmed.**

Complete M/[(q-1)/q] at n=16 DMRG (+ exact diag for q=6,8):

| q | (q-1)/q | M/[(q-1)/q] | type | tail b |
|---|---------|-------------|------|--------|
| 2 | 0.500 | 1.352 | real | 1.98 |
| 3 | 0.667 | 1.086 | real | 2.01 |
| 4 | 0.750 | 1.003 | boundary | 2.02 |
| 5 | 0.800 | 0.964 | walking | 2.01 |
| 6 | 0.833 | 0.910* | marginal | — |
| 7 | 0.857 | 0.930 | broken | 2.07 |
| 8 | 0.875 | 0.895* | broken | — |

(*exact diag small n, not DMRG)

q_cross converges to q≈4.0 as n→∞. Universal tail exponent b=2.018±0.029.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unused for 65 sprints. Measure entanglement entropy at g_c on real hardware (q=2, n=6-8). Prediction: S matches simulator to ~10%.
2. **Thermodynamic limit of q_cross** — DMRG at n=32,40 for q=3,4,5 to see if q_cross → 4.0 exactly. Currently 3.92 at n=24.
3. **Entanglement Hamiltonian H_E across walking boundary** — H_E = -log(ρ_A) structure. Do the entanglement energies show the same crossover? Is H_E more/less local at walking vs broken?

## What's Been Ruled Out
- Tail weight growth as walking-specific (universal n^2 for all q, Sprints 088-090)
- Sprint 087 exponents 1.72, 2.52 (real-space fitting artifact)
- q=7 b≈3.0 exponent (pre-asymptotic, converges to b≈2.1, Sprint 089a)
- Dynamic transfer rates as walking discriminator (not monotonic, Sprint 089b)
- α=3 as Re(c) recovery probe (single-size extraction artifact, Sprint 086)
- c_∞ (min-entropy) as Re(c) recovery probe (degrades like c_1)
- Open-BC CC profile fit for Rényi c_α extraction (R²<0.9 at n≤12)
- x_σ as walking discriminator (nearly constant across q)
- Casimir energy as walking discriminator (follows Re(c) for all q)
- Entanglement gap as walking discriminator (INCREASES with q at fixed n)
- q=6 as second walking case (c_eff drops 2.9% n=8→12)

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=2,3 (very fast, n≤24+), q=4 (moderate, n≤24 in 255s), q=5 (fast, n≤24), q=7 (n≤16)
- Periodic exact diag correlator: q=2 n≤14, q=3 n≤10, q=5 n≤8, q=7 n≤7
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
