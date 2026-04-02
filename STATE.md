# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 091 — Entanglement Hamiltonian H_E across walking boundary. Three experiments: (1) BW fidelity for q=2-7 at g_c=1/q periodic — confounded by different nA. (2) Fixed nA=4 comparison: non-Potts fraction grows 0.03% (q=2) → 2.7% (q=5) exponentially. Biggest jump at q=3→4 (11.7×). NNN/3-body add <0.2%. (3) nA scaling: q=5 non-Potts grows 1.7× faster than q=2. BW alpha ≈ 3.2 for q≥4.

## Active Research Thread
**Entanglement Hamiltonian structure across walking/CFT boundary.**

Fixed nA=4 BW comparison (091b):

| q | Potts NN locality | non-Potts | BW fidelity | BW alpha |
|---|---|---|---|---|
| 2 | 99.971% | 0.029% | 0.9997 | 2.390 |
| 3 | 99.916% | 0.081% | 0.9994 | 2.401 |
| 4 | 99.000% | 0.955% | 0.9942 | 3.228 |
| 5 | 97.075% | 2.729% | 0.9832 | 3.195 |

Biggest non-Potts jump: q=3→4 (11.7×), coincides with real-to-complex CFT boundary.
q=5 non-Potts grows 1.7× faster with nA than q=2.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unused for 66 sprints. Measure entanglement entropy at g_c on real hardware (q=2, n=6-8). Prediction: S matches simulator to ~10%.
2. **Large-nA BW comparison q=2 vs q=5** — DMRG at n=16-24 (nA=8-12) to see if the 1.7× slope ratio holds. Would confirm walking amplifies BW corrections at accessible DMRG sizes.
3. **Non-Potts operator identification** — What ARE the non-Potts operators in H_E? For q=2 (Pauli basis available): decompose the BW residual into Pauli terms. Identify which local operators dominate.

## What's Been Ruled Out
- 091a direct comparison across q (confounded by different nA)
- NNN/3-body Potts operators as source of BW deviation (<0.2%)
- BW fidelity as walking discriminator (dominated by nA, not q)
- (Prior sprints) All items from Sprint 090 ruled-out list still apply

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=2,3 (very fast, n≤24+), q=4 (moderate, n≤24), q=5 (fast, n≤24), q=7 (n≤16)
- Periodic exact diag correlator: q=2 n≤14, q=3 n≤10, q=5 n≤8, q=7 n≤7
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
