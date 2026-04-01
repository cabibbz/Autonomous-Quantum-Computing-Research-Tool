# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 048 — ν(q=10) Extraction: Dead-Pair Bias Invalidates Raw MI-CV at Large d

## Active Research Thread
MI-CV as universal phase transition classifier. **Major methodological finding**: raw MI-CV has a dead-pair bias at large d that creates spurious crossing absence. At q=10 (d=10), n=8 has 25% dead pairs vs n=12 has 17%, systematically deflating n=12 CV. Filtered MI-CV (excluding near-zero MI pairs) shows n=12 > n=8 everywhere — same BKT-like pattern as q=4.

Sprint 043's "q=10 crossing confirmed" is INVALIDATED (χ=10 artifact + dead-pair bias).

Filtered ν(q) picture may need revision: q=4 AND q=10 both BKT-like, with standard FSS at q=5,7 between them.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **Half-chain entropy FSS** — MI-CV is unreliable at large d due to χ non-convergence. Half-chain entropy S(n/2) requires only the bipartition spectrum, not all-pairs MI. Test at q=4,5,7,10 as an alternative FSS probe for ν extraction.
2. **Re-examine q=7 with filtered MI-CV + χ check** — d=7 requires χ ≈ 49 for MI convergence. Sprint 047 used χ=20. Quick: rerun single point at χ=40,60 to test MI convergence.
3. **Correlation length from transfer matrix** — Extract ξ directly from MPS transfer matrix eigenvalues. ξ ∝ |g-g_c|^(-ν) gives ν without all-pairs MI. Works at any d and χ.

## What's Been Ruled Out
- Small-scale QEC active correction: exhaustively proven to fail for [[5,1,3]] (Sprints 026-028)
- Pauli fraction as BW metric: fails at large n (Sprint 035)
- NN-only MI-CV: dominated by boundary effects (Sprint 036)
- BW sin_inv envelope: wrong at finite size (Sprint 035)
- Simple power-law fit for ν at n≤32: corrections to scaling too large (Sprint 037)
- 1D quantum Potts first-order for ANY q: all tested q (2-20) show second-order (Sprint 043)
- Clock ≡ Potts for q≥4: false, models differ (Sprint 041)
- χ=10 for d≥7: NOT converged (Sprint 044)
- Gell-Mann MI for q≥4: unreliable, systematic errors (Sprints 045-046)
- ν monotonically increasing for q>3: WRONG, ν peaks at q=4 (Sprint 047)
- Standard FSS collapse for q=4: fails, ν→∞ when g_c constrained (Sprint 046)
- **Raw MI-CV for d≥10: dead-pair bias makes size comparison unreliable (Sprint 048)**
- **Sprint 043 q=10 crossing: INVALIDATED (χ=10 + dead-pair artifact) (Sprint 048)**

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- Direct MPS contraction MI: all-pairs MI for ANY d (Sprint 043) — ONLY reliable method for d≥4
- Custom PottsChain model: Kronecker-delta coupling for any q (Sprint 042)
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
- g_c scaling law: g_c(q) = 0.87*(q-3)^(-0.85) for q≥4 (Sprint 044)
- Data collapse pipeline: joint (ν, g_c) optimization (Sprint 045)
- **Dead-pair filtering for MI-CV at large d (Sprint 048)**
