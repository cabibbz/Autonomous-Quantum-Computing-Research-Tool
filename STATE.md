# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 045 — ν(q) Extraction: ν(q=5)≈2.0, universality class changes with q

## Active Research Thread
MI-CV as universal phase transition classifier + universality class extraction via data collapse. Major finding: ν(q) is NON-MONOTONIC — decreases from q=2 (ν=1) to q=3 (ν=5/6), then INCREASES sharply to q=5 (ν≈2). This is opposite to mean-field prediction. Large ν at q=5 means transition is "almost first-order" — the tendency toward first-order manifests as diverging ν rather than discontinuity.

Also discovered: Sprint 042d Gell-Mann MI reconstruction was wrong for q=5 (errors up to 11x). Direct MPS contraction (Sprint 043) is essential for d≥5.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **ν(q=7) or ν(q=10) extraction** — Does ν keep increasing with q? Need larger n for d≥7 (n=12 DMRG at d=7 takes >120s, may need lower chi or Z_q conserve). If ν→∞ as q→∞, that's the 1D quantum analog of first-order.
2. **q=4 ν extraction** — q=4 is marginal in 2D classical (logarithmic corrections). What is ν(4) in 1D quantum? Should be between ν(3)=5/6 and ν(5)=2. If ν(4)≈1, the ν non-monotonicity starts exactly at q=4.
3. **Recheck Sprint 042 conclusion with direct MPS** — The q=5 Potts "5.7x steeper than clock" claim used Gell-Mann MI for Potts. Need to verify clock vs Potts slope comparison with direct MPS.

## What's Been Ruled Out
- Small-scale QEC active correction: exhaustively proven to fail for [[5,1,3]] (Sprints 026-028)
- Pauli fraction as BW metric: fails at large n (Sprint 035)
- NN-only MI-CV: dominated by boundary effects (Sprint 036)
- BW sin_inv envelope: wrong at finite size (Sprint 035)
- Simple power-law fit for ν at n≤32: corrections to scaling too large (Sprint 037)
- 1D quantum Potts first-order for ANY q: all tested q (2-20) show second-order (Sprint 043)
- Clock ≡ Potts for q≥4: false, models differ (Sprint 041)
- χ=10 for d≥7: NOT converged, 25% CV inflation (Sprint 044)
- **Gell-Mann MI for q≥5: WRONG, up to 11x error (Sprint 045)**
- **ν decreasing with q: WRONG, ν increases for q>3 (Sprint 045)**

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- Direct MPS contraction MI: all-pairs MI for ANY d (Sprint 043) — ONLY reliable method for d≥5
- Custom PottsChain model: Kronecker-delta coupling for any q (Sprint 042)
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
- g_c scaling law: g_c(q) = 0.87*(q-3)^(-0.85) for q≥4 (Sprint 044)
- Data collapse pipeline: joint (ν, g_c) optimization (Sprint 045)
