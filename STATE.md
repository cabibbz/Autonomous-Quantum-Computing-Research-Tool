# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 092 — Operator decomposition of entanglement Hamiltonian H_E. Three experiments: (1) Pauli decomposition at q=2 nA=3-6: non-BW dominated by ZXZ (3-body) and YY (2-body), growing from 0.005% to 6.3% with nA. (2) Clock-shift decomposition at fixed nA=3 for q=2-5: non-BW flat at 0.005-0.009%, all q. (3) nA=4 for q=2,3: same operator types (MM and DFD) dominate non-BW for both q. Walking boundary doesn't change operator TYPES, only amplitudes.

## Active Research Thread
**Entanglement Hamiltonian operator content across walking/CFT boundary.**

Key finding: BW corrections are dominated by Mixed operators (M = X^a·Z^b, generalizing Pauli Y) and 3-body DFD/DMD chains. These represent φ·π (position-momentum) entanglement correlations absent from the physical Hamiltonian. No qualitative change at q=4 walking boundary — just amplitude growth.

non-BW operator weight scaling with nA (q=2):
| nA | non-BW weight | dominant |
|----|--------------|----------|
| 3  | 0.005%       | DFD      |
| 4  | 0.008%       | MM, DFD  |
| 5  | 2.8%         | ZZZ, ZX  |
| 6  | 6.3%         | ZXZ, ZZZ |

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unused for 67 sprints. Measure entanglement entropy S(n/2) at g_c on real hardware for q=2 n=6-8. Prediction: S matches simulator to ~10%. Would validate all finite-size results.
2. **nA scaling of operator types q=2 nA=5-7** — Do the MM and DFD operators that dominate at nA=4 continue to grow, or do new types emerge at large nA? Use Pauli basis (q=2) for maximum nA reach.
3. **Entanglement temperature profile** — BW predicts H_E = 2π Σ_i β(i) h_i with β(i) = sin-envelope. Do the BW coefficients we extracted in 092 follow this envelope? Would test BW at the operator level, not just fidelity.

## What's Been Ruled Out
- 092b: q-dependence of non-BW weight at nA=3 (essentially flat, <2× variation)
- 092c: Qualitative difference in non-BW operator types between q=2 and q=3 (same MM/DFD)
- Walking boundary introducing new operator types (it amplifies existing ones)
- All items from Sprint 091 ruled-out list still apply

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=2,3 (very fast, n≤24+), q=4 (moderate, n≤24), q=5 (fast, n≤24), q=7 (n≤16)
- Periodic exact diag correlator: q=2 n≤14, q=3 n≤10, q=5 n≤8, q=7 n≤7
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
