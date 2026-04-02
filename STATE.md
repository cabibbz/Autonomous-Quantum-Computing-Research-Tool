# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 095 — BW R² vs n at fixed nA: threshold persistence test. Three experiments: (1) q=2 periodic nA=3-6 with n up to 18: BW accuracy flat in n for nA≤4, collapses at nA=6 regardless of n. (2) q=5 periodic nA=3 n=7,8: walking amplification 10-15× independent of n; DMRG too expensive (chi≥125 needed). (3) Compilation: equal-bipartition (nA/n=0.5) has 2-18× penalty, peaking at nA=5. Sprint 094's nA=5 "threshold" was artifact — real threshold at nA=6.

## Active Research Thread
**BW entanglement Hamiltonian structure across walking/CFT boundary.**

Sprint 095 answered the key question: BW threshold depends on nA (absolute subsystem size), NOT nA/n ratio. The threshold is a UV lattice effect. Equal bipartition has an anomalous penalty (2-18×).

Corrected BW R² table (1-R², best non-equal-bipartition value):
| q | nA=3 | nA=4 | nA=5 | nA=6 |
|---|------|------|------|------|
| 2 | 1.5e-4 | 2.5e-4 | 5.6e-4 | **1.7e-1** |
| 5 | 1.7e-3 | 3.4e-2* | — | — |

*nA=4 q=5 only measured at equal bipartition (S094), likely 5-10× better at n≫nA.

Walking amplification (q=5/q=2) at nA=3: ~10-15×, ratio-independent.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **BW threshold mechanism** — WHY does BW collapse at nA=6 for q=2? Is it when the entanglement spectrum tail weight exceeds some critical fraction? Check: at nA=6 (threshold), how does the non-BW operator content compare to nA=5 (safe)?
2. **Hardware validation** — 580s QPU unused for 70 sprints. Measure entanglement entropy S(n/2) at g_c on real hardware for q=2 n=6-8. Or: BW fidelity at nA=3 on hardware (simplest prediction: R²>0.999).
3. **q=3,4 ratio dependence** — Verify equal-bipartition penalty for q=3,4. Predict: nA=5 penalty is ~18× for q=2, should be larger for q=3,4 (closer to threshold). Quick experiment: periodic q=3 n=11-16 with nA=5.

## What's Been Ruled Out
- Equal-bipartition artifact as sole cause of BW threshold: NO, nA=6 threshold persists at nA/n=0.33.
- Sprint 094 nA=5 threshold: ARTIFACT of equal bipartition (18× penalty). True nA=5 accuracy ~5e-4 for q=2.
- DMRG for BW at q=5: needs chi≥q^nA=125+, too expensive for large n.
- Power-law fit to 1-R²(nA): doesn't fit. Threshold/two-regime behavior.
- All previously ruled-out items still apply.

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=2,3 (very fast, n≤24+), q=4 (moderate, n≤24), q=5 (fast, n≤24, but chi≤60)
- Periodic exact diag: q=2 n≤18, q=3 n≤10, q=5 n≤8
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
