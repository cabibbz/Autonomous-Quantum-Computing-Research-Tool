# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 094 — BW R² vs subsystem size nA for q=2-5. Four experiments: (1) q=2 nA=3-7: threshold at nA=5 (19× jump). (2) q=3 nA=3-6, q=4 nA=3-5: threshold at nA=5 for both (261× for q=3, 178× for q=4). (3) q=5 nA=3-4: threshold at nA=4 (33×), nA=5 infeasible (24GB memory). (4) Compilation: B(q) = 0.48q + 1.09 (R²=0.999). Walking shifts nA* down by 1 site.

## Active Research Thread
**BW entanglement Hamiltonian structure across walking/CFT boundary.**

Sprint 094 mapped the full nA-dependence of BW fidelity. Key result: BW breakdown has THRESHOLD behavior (not power law). The threshold nA*(q) decreases with q: q=2,3,4 → nA*=5; q=5 → nA*=4. The exponential rate B scales linearly with q. Walking amplification is nA-dependent: mild at nA=3 (3.5×), dramatic at nA=4 (64×).

Complete BW R² table (1-R²):
| q | nA=3 | nA=4 | nA=5 | nA=6 | nA=7 |
|---|------|------|------|------|------|
| 2 | 2.9e-4 | 5.3e-4 | 1.0e-2 | 2.4e-1 | 3.6e-1 |
| 3 | 5.6e-4 | 8.4e-4 | 2.2e-1 | 4.4e-1 | — |
| 4 | 8.1e-4 | 1.8e-3 | 3.1e-1 | — | — |
| 5 | 1.0e-3 | 3.4e-2 | — | — | — |

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unused for 69 sprints. Measure entanglement entropy S(n/2) at g_c on real hardware for q=2 n=6-8. Prediction: S matches simulator to ~10%. Would validate all finite-size results.
2. **BW threshold mechanism** — WHY does BW collapse at nA*≈5 for q=2-4 and nA*≈4 for q=5? Is it the number of entanglement spectrum levels exceeding BW-compatible subspace? Check: at nA*, does the tail weight (beyond (q-1) multiplet) suddenly dominate?
3. **DMRG BW at large n** — Use DMRG to get ρ_A at n=20-24 with nA=3-5 (subsystem of larger chain, NOT equal partition). Does the threshold persist when n≫nA? If BW improves, the threshold is a finite-n/nA artifact.

## What's Been Ruled Out
- Power-law fit to 1-R²(nA): doesn't fit (R²=0.42-0.80). Threshold/two-regime behavior.
- q=5 nA=5 exact diag: infeasible (lil_matrix construction >24GB at dim=9.8M, not eigsolver).
- All previously ruled-out items still apply.

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=2,3 (very fast, n≤24+), q=4 (moderate, n≤24), q=5 (fast, n≤24), q=7 (n≤16)
- Periodic exact diag correlator: q=2 n≤14, q=3 n≤10, q=5 n≤8, q=7 n≤7
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
