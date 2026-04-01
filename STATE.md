# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 044 — g_c Scaling Law: g_c ≈ 0.87*(q-3)^(-0.85) for q≥4, verified at q=7

## Active Research Thread
MI-CV as universal phase transition classifier. Extracted the scaling law g_c(q) for 1D quantum Potts. Best fit: g_c = 0.87*(q-3)^(-0.85) for q≥4, with g_c=1 for q=2,3 (self-duality). Pole at q=3 reflects self-duality breaking. Blind prediction of g_c(7)=0.263 confirmed within 1.6% (measured 0.259). Large-q regime (g_c flattening) begins around q≈7.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **Exponent 0.85 ≈ 5/6? Connection to Potts universality** — Is the g_c scaling exponent exactly ν_Potts=5/6? This would be a deep connection between the scaling of g_c with q and the critical exponents of the q=3 Potts universality class.
2. **q=10 data collapse with direct MPS** — With the scaling law established, use the predicted g_c(q) to do data collapse at q=10 and extract ν(q). Does ν depend on q?
3. **2D quantum Potts on ladder** — 1D quantum is always second-order. Does adding a second chain (ladder geometry) restore first-order for q>4? Tests whether spatial dimensionality controls transition order.

## What's Been Ruled Out
- Small-scale QEC active correction: exhaustively proven to fail for [[5,1,3]] (Sprints 026-028)
- Pauli fraction as BW metric: fails at large n (Sprint 035)
- NN-only MI-CV: dominated by boundary effects (Sprint 036)
- BW sin_inv envelope: wrong at finite size (Sprint 035)
- Simple power-law fit for ν at n≤32: corrections to scaling too large (Sprint 037)
- ν=0.755 from crossing points: artifact (Sprint 038)
- Fixed g_c collapse: misleading due to finite-size g_c shift (Sprint 039)
- n=24 qutrit DMRG without symmetry: too slow (>60s per point) (Sprint 039)
- 1D quantum Potts first-order for ANY q: all tested q (2-20) show second-order (Sprint 043)
- Clock ≡ Potts for q≥4: false, models differ (Sprint 041)
- χ=10 for d≥7: NOT converged, 25% CV inflation (Sprint 044)
- 1/√q and simple exponential for g_c(q): too slow/fast decay (Sprint 044)

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- Direct MPS contraction MI: all-pairs MI for ANY d (Sprint 043)
- Custom PottsChain model: Kronecker-delta coupling for any q (Sprint 042)
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
- g_c scaling law: g_c(q) = 0.87*(q-3)^(-0.85) for q≥4 (Sprint 044)
