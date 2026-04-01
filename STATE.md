# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 049 — Correlation Length & Entropy FSS: q=3 Potts Critical Point WRONG

## Active Research Thread
**CRITICAL FINDING: q=3 Potts g_c ≈ 0.33, NOT 1.0.** Exact diag confirmed S=0.047 at g=1.0 n=8 — deep disordered phase. Entropy peaks at g≈0 (S=ln3, GHZ) and drops sharply by g≈0.30-0.35. Central charge c(g) from pairwise n=8,12 estimates crosses exact c=4/5 near g≈0.33.

**Impact: Sprints 038-048 Potts results (g_c, scaling law, ν extraction) were at the WRONG g.** MI-CV crossings near g≈1.0 are a disordered-phase crossover, not the phase transition.

**TFIM results validated.** Central charge c = 0.52 → 0.50 (5 sizes). Correlation length method works but finite-size limited (ν_eff = 0.62-0.72). Entropy FSS unsuitable for ν (logarithmic singularity).

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **Find true g_c for q=3-10 Potts** — Use entropy peak / central charge method at each q. Exact diag for small n, DMRG for n=12-16. This is prerequisite for ALL Potts physics.
2. **Re-examine MI-CV at true g_c** — Do MI-CV crossings appear near the real critical points? If yes, prior qualitative signatures (crossing=2nd order, etc.) are rescued.
3. **iDMRG for correlation length** — Infinite MPS avoids finite-size saturation. TeNPy supports iDMRG. Would give ξ→∞ at criticality without boundary effects.

## What's Been Ruled Out
- Small-scale QEC active correction: fails for [[5,1,3]] (Sprints 026-028)
- Pauli fraction as BW metric: fails at large n (Sprint 035)
- NN-only MI-CV: boundary effects (Sprint 036)
- BW sin_inv envelope: wrong at finite size (Sprint 035)
- 1D quantum Potts first-order for ANY q: all tested q show second-order (Sprint 043)
- Clock ≡ Potts for q≥4: false (Sprint 041)
- χ=10 for d≥7: NOT converged (Sprint 044)
- Gell-Mann MI for q≥4: unreliable (Sprints 045-046)
- Raw MI-CV for d≥10: dead-pair bias (Sprint 048)
- Entropy FSS for ν extraction: logarithmic singularity defeats power-law collapse (Sprint 049)
- Naive ξ extraction for ν: finite-size saturation gives ν_eff << ν_true (Sprint 049)
- **g_c = 1.0 for q=3 Potts with our Hamiltonian: WRONG, g_c ≈ 0.33 (Sprint 049)**
- **All prior Potts g_c values (Sprints 038-048): measured at wrong g (Sprint 049)**

## Key Tools Available
- Exact diag: n≤10 (n≤8 for q=3, n≤6 for q=4)
- DMRG (TeNPy): n>10, 1D only. TFIChain has Z₂ conservation. PottsChain has NO conservation.
- Direct MPS contraction MI: all-pairs MI for ANY d — ONLY reliable for d≥4
- Custom PottsChain model: Kronecker-delta coupling for any q
- IBM QPU: 580s remaining
- Central charge extraction: S vs ln(n) at g_c — validated for TFIM (Sprint 049)
