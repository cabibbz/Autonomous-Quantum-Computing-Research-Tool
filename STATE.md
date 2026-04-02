# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 093 — BW entanglement temperature profile. Three experiments: (1) Position-dependent BW coefficients for q=2,3,5: sin-envelope shape confirmed (R²>0.99), field/bond ratio ≠ g_c. (2) Global BW Frobenius R² across q=2-7: degrades monotonically with q, 34× amplification at q=5 from nA=3→4 vs 1.8× at q=2. (3) BW residual spatial profile: bulk-concentrated for real CFT (q=2,3), uniform for walking (q=5).

## Active Research Thread
**BW entanglement Hamiltonian structure across walking/CFT boundary.**

Sprint 093 established the SPATIAL structure of BW — the sin-envelope works for all q, but the residual has a q-dependent spatial distribution. Key new finding: walking makes BW corrections uniform in space (boundary enrichment ≈ 1.0), while real CFT concentrates corrections in the bulk (enrichment ≈ 0.4).

BW Frobenius R² at nA=4:
| q | R² | 1-R² | residual profile |
|---|-----|------|-----------------|
| 2 | 0.9995 | 5e-4 | bulk-concentrated (0.42×) |
| 3 | 0.9992 | 8e-4 | bulk-concentrated (0.47×) |
| 4 | 0.9984 | 2e-3 | — |
| 5 | 0.9655 | 3e-2 | uniform (0.98×) |

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unused for 68 sprints. Measure entanglement entropy S(n/2) at g_c on real hardware for q=2 n=6-8. Prediction: S matches simulator to ~10%. Would validate all finite-size results.
2. **BW residual operator decomposition at larger nA** — At q=2 nA=6, residual is 19% of H_E. What operators dominate the residual? Are they still MM/DFD (Sprint 092) or do new types emerge? Could identify the CFT operators responsible for BW breakdown.
3. **nA scaling of BW R²** — Map R²(nA) for q=2-5 with nA=3-7 (q=2) and nA=3-5 (q=3-5). Extract the walking amplification exponent: does 1-R² ~ nA^gamma with gamma(q) depending on q? Would quantify the walking-BW connection.

## What's Been Ruled Out
- 093a: nA=3 is too small for meaningful BW envelope analysis (only 2 bonds)
- 093b: 2-parameter fit (separate coupling/field alpha) is WORSE — BW's single envelope is the right decomposition
- Walking boundary doesn't change BW operator types (Sprint 092) — it changes amplitudes and spatial profile (Sprint 093)
- All previously ruled-out items still apply

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=2,3 (very fast, n≤24+), q=4 (moderate, n≤24), q=5 (fast, n≤24), q=7 (n≤16)
- Periodic exact diag correlator: q=2 n≤14, q=3 n≤10, q=5 n≤8, q=7 n≤7
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
