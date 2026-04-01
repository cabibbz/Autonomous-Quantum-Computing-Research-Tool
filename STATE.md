# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 047 — ν(q=7) Extraction: Crossings Return, q=4 Confirmed as BKT Peak

## Active Research Thread
MI-CV as universal phase transition classifier + universality class extraction via data collapse. Major finding: q=7 shows **crossing curves** (unlike q=4), and ν(q=7) ≈ 0.5 (near mean-field). This confirms q=4 as the unique BKT-like peak in the ν(q) curve.

Updated ν(q) picture: q=2 (1.0) → q=3 (5/6, minimum) → q=4 (≥2.2, BKT peak, NO crossings) → q=5 (2.0) → q=7 (≈0.5, crossings return). q=4 is the ONLY value without MI-CV crossings.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **ν(q=10) extraction** — Does ν continue decreasing toward mean-field (1/2)? If ν(10)≈0.5, confirms large-q → mean-field. Quick: reuse Sprint 043 g_c=0.246, need MI-CV at n=8,12.
2. **BKT test for q=4 refined** — With q=7 confirming q=4 is special, test BKT scaling directly: ξ ~ exp(c/√|g-g_c|). Need finer g-grid near g_c≈0.89.
3. **q=4 entanglement entropy peak** — Verify g_c(4) independently via half-chain entropy. MI-CV failed at q=4 (no crossings), but entropy peak should still work.

## What's Been Ruled Out
- Small-scale QEC active correction: exhaustively proven to fail for [[5,1,3]] (Sprints 026-028)
- Pauli fraction as BW metric: fails at large n (Sprint 035)
- NN-only MI-CV: dominated by boundary effects (Sprint 036)
- BW sin_inv envelope: wrong at finite size (Sprint 035)
- Simple power-law fit for ν at n≤32: corrections to scaling too large (Sprint 037)
- 1D quantum Potts first-order for ANY q: all tested q (2-20) show second-order (Sprint 043)
- Clock ≡ Potts for q≥4: false, models differ (Sprint 041)
- χ=10 for d≥7: NOT converged, 25% CV inflation (Sprint 044)
- Gell-Mann MI for q≥4: unreliable, systematic errors (Sprints 045-046)
- ν monotonically increasing for q>3: WRONG, ν peaks at q=4 then decreases (Sprint 047)
- Standard FSS collapse for q=4: fails, ν→∞ when g_c constrained (Sprint 046)

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- Direct MPS contraction MI: all-pairs MI for ANY d (Sprint 043) — ONLY reliable method for d≥4
- Custom PottsChain model: Kronecker-delta coupling for any q (Sprint 042)
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
- g_c scaling law: g_c(q) = 0.87*(q-3)^(-0.85) for q≥4 (Sprint 044)
- Data collapse pipeline: joint (ν, g_c) optimization (Sprint 045)
