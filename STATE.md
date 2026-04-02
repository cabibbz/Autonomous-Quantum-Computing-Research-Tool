# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 097 — H_E Operator Compactness: BW Breakdown is Fundamental. Three experiments: (1) Full Pauli decomposition q=2 n=14 nA=3-6: H_E described by exactly 2nA-1 BW operators pre-threshold, explodes to ~1000 at nA=6. (2) Non-BW structure at threshold: diffuse (PR=357), 743 operators for 90% of non-BW weight. XZ-type mixed operators dominate (44%). (3) q=3 n=10 clock-shift basis: same pattern. Non-Potts fraction even larger (22.4% vs 16.3%).

## Active Research Thread
**BW entanglement Hamiltonian structure across walking/CFT boundary.**

Sprint 097 closed the BW mechanism story: BW is the BEST compact approximation of H_E. Pre-threshold, 2nA-1 operators suffice (>99.9%). At threshold, H_E develops diffuse non-Potts content that cannot be captured by any compact ansatz.

BW thread summary (Sprints 091-097):
- BW accuracy depends on nA (UV lattice effect), not nA/n ratio
- BW breaks when non-Potts (mixed clock-shift) operators exceed ~16-22% of H_E
- Non-Potts threshold coincides with BW threshold (sharp, 244× jump)
- Adding more Potts operators (any range) barely helps (3-5% improvement)
- Non-BW content is intrinsically diffuse (PR=357) — no compact fix exists
- Universal across q=2,3. Walking shifts threshold down (q=5: nA*≈4, q=2: nA*≈6)

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU budget unused for 72 sprints. Strongest prediction: BW R²>0.999 at nA=3 on real hardware for q=2. Measure S(nA) at g_c on real QPU. Or: entanglement spectrum on hardware.
2. **New direction: 2D or higher-q BW** — BW story for 1D S_q Potts is complete. Move to 2D (Ly=2 cylinder) or test BW for hybrid model. Different operator algebras may show different thresholds.
3. **Harden Casimir finding (Sprint 083)** — c_implied/Re(c) ≈ 1.00 for ALL q, but only 3 data points at q≥7. Extend with GPU.

## What's Been Ruled Out
- "BW + corrections" ansatz: RULED OUT. Non-BW content is diffuse (PR=357), not compact.
- Entanglement spectrum as BW predictor: spectrum is smooth across threshold (Sprint 096).
- Power-law fit to 1-R²(nA): doesn't fit. Threshold/two-regime behavior.
- DMRG for BW at q=5: needs chi≥q^nA=125+, too expensive.
- All previously ruled-out items still apply.

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=2,3 (fast, n≤24+), q=4 (moderate, n≤24), q=5 (fast, n≤24, chi≤60)
- Periodic exact diag: q=2 n≤18, q=3 n≤10, q=5 n≤8
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
