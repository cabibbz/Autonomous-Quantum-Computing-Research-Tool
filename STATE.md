# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 099 — Complex CFT: Im(c) Oscillation Detection (NEGATIVE RESULT). Dense Casimir scan at all integer N for q=2,5,7. Oscillatory model doesn't beat monotonic 1/N⁶. Predicted period (6.5 in ln N for q=5) requires N~4-2700 — our sizes cover 14% of a cycle. Pairwise non-monotonicity in c/Re(c) for q=5 traced to velocity, not Casimir. Richardson extrapolation drifts for q>4.

## Active Research Thread
**Casimir/walking thread COMPLETE.** Two confirmed novel findings:
1. Casimir energy obeys Re(c) across walking boundary (16× more consistent than entropy)
2. Walking breakdown is exclusively an entropy phenomenon at finite size

**BW thread COMPLETE.** Im(c) oscillation attempt: undetectable at accessible sizes.

Need a NEW research direction.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unused for 74 sprints. Strongest prediction: Casimir E₀ at g_c for q=2 on QPU (5-qubit periodic Ising). Can we see πvc/(6N²) correction? Also: BW R²>0.999 at nA=3 on real hardware.
2. **Entanglement Hamiltonian in 2D** — BW story complete in 1D. Test on Ly=2 cylinder (DMRG accessible). Does BW threshold shift? Is non-Potts content different?
3. **DMRG Casimir precision** — Can DMRG extract E₀(N) precisely enough to detect Im(c) oscillations at N=10-100? Would need <10⁻⁶ relative precision. Test with q=2 (known answer) first.

## What's Been Ruled Out
- Im(c) oscillation detection via exact diag: RULED OUT (Sprint 099) — sizes too small
- All BW correction approaches: RULED OUT (Sprint 097)
- Casimir hardening: DONE, confirmed (Sprint 098)
- All previously ruled-out items still apply

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=7
- Dense Casimir data: q=2 N=4-14 (all), q=5 N=4-10 (all), q=7 N=4-8 (all)
- S_q Potts DMRG: q=2,3 (fast, n≤24+), q=5 (n≤24, chi≤60)
- IBM QPU: 580s remaining
