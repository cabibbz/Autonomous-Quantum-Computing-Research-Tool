# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 046 — ν(q=4) Extraction: Marginal Point Shows BKT-like Behavior

## Active Research Thread
MI-CV as universal phase transition classifier + universality class extraction via data collapse. Major finding: ν(q=4) ≥ 2.2 and INCREASING with system size (1.92 → 2.71 from n=8,12 to n=12,16). This suggests BKT-like behavior where ν is effectively infinite. Standard FSS collapse fails at q=4 — curves don't cross near g_c.

Updated ν(q) picture: q=2 (1.0) → q=3 (5/6, minimum) → q=4 (≥2.2, diverging?) → q=5 (2.0). The marginal point q=4 separates standard second-order (q≤3) from large-ν second-order (q≥5), possibly via a BKT-like mechanism.

Also: Sprint 040 Gell-Mann data at d=4 was partially misleading (showed crossings at g_c≈0.89 that don't exist with direct MPS). Gell-Mann now considered unreliable for d≥4.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **ν(q=7) or ν(q=10) extraction** — Does ν decrease back from q=5's ν=2? If ν(7)<ν(5), it confirms q=4 is a peak. Need χ≥20, may need n=8,12 only (n=16 at d=7 is slow).
2. **BKT test for q=4** — Test BKT scaling: ξ ~ exp(c/√|g-g_c|). Plot ln(dCV/dg) vs 1/√(g-g_c) at g_c=0.89. If linear, confirms BKT. Need more g points near g_c.
3. **Verify g_c(4) with entanglement entropy** — Half-chain entropy peak should mark g_c independently of MI-CV. Compare peak position with g_c≈0.89 from scaling law.

## What's Been Ruled Out
- Small-scale QEC active correction: exhaustively proven to fail for [[5,1,3]] (Sprints 026-028)
- Pauli fraction as BW metric: fails at large n (Sprint 035)
- NN-only MI-CV: dominated by boundary effects (Sprint 036)
- BW sin_inv envelope: wrong at finite size (Sprint 035)
- Simple power-law fit for ν at n≤32: corrections to scaling too large (Sprint 037)
- 1D quantum Potts first-order for ANY q: all tested q (2-20) show second-order (Sprint 043)
- Clock ≡ Potts for q≥4: false, models differ (Sprint 041)
- χ=10 for d≥7: NOT converged, 25% CV inflation (Sprint 044)
- Gell-Mann MI for q≥5: WRONG, up to 11x error (Sprint 045)
- **Gell-Mann MI for d=4: small systematic errors create artificial crossings (Sprint 046)**
- ν decreasing with q: WRONG, ν increases for q>3 (Sprint 045)
- **Standard FSS collapse for q=4: fails, ν→∞ when g_c constrained to 0.89 (Sprint 046)**

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- Direct MPS contraction MI: all-pairs MI for ANY d (Sprint 043) — ONLY reliable method for d≥4
- Custom PottsChain model: Kronecker-delta coupling for any q (Sprint 042)
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
- g_c scaling law: g_c(q) = 0.87*(q-3)^(-0.85) for q≥4 (Sprint 044) — but g_c(4) input may be off
- Data collapse pipeline: joint (ν, g_c) optimization (Sprint 045)
