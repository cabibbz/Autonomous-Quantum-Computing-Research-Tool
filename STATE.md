# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 043 — First-Order Search: Large q Potts (q=10, q=20) — all second-order

## Active Research Thread
MI-CV as universal phase transition classifier. Tested 1D quantum Potts at q=10 (crossing confirmed at g_c≈0.246) and q=20 (identical to q=10 at χ=10, continuous entropy). The 2D classical "q>4 → first-order" rule does NOT apply in 1D quantum for any tested q. Also developed direct MPS contraction for ρ_ij, enabling MI-CV at arbitrarily large d.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **g_c scaling law** — Plot g_c vs q for q=2,3,4,5,10 and fit. Is g_c ~ 1/q, 1/√q, or 1/log(q)? This would predict g_c for arbitrary q and connect to the quantum-classical mapping anisotropy.
2. **q=10 data collapse** — With n=8,12 data at crossing, attempt data collapse to extract ν. Does ν match any known universality class, or does it depend on q?
3. **2D quantum Potts on ladder** — 1D quantum is always second-order. Does adding a second chain (ladder geometry) restore first-order for q>4? This would pinpoint whether spatial dimensionality or quantum kinematics controls transition order.

## What's Been Ruled Out
- Small-scale QEC active correction: exhaustively proven to fail for [[5,1,3]] (Sprints 026-028)
- Pauli fraction as BW metric: fails at large n (Sprint 035)
- NN-only MI-CV: dominated by boundary effects (Sprint 036)
- BW sin_inv envelope: wrong at finite size (Sprint 035)
- Simple power-law fit for ν at n≤32: corrections to scaling too large (Sprint 037)
- ν=0.755 from crossing points: artifact (Sprint 038)
- Fixed g_c collapse: misleading due to finite-size g_c shift (Sprint 039)
- n=24 qutrit DMRG without symmetry: too slow (>60s per point) (Sprint 039)
- q=4 marginal showing dome: shows crossing at n=8,12 (Sprint 040)
- q=5 clock showing no crossings: still shows crossings at g≈0.67 (Sprint 041)
- ClockChain as Potts proxy for q≥4: clock ≠ Potts at q≥4 (Sprint 041)
- q=5 Potts first-order: shows second-order crossings at g≈0.41 (Sprint 042)
- **q=10 Potts first-order: shows second-order crossings at g≈0.25 (Sprint 043)**
- **q=20 Potts first-order: continuous entropy, universal with q=10 (Sprint 043)**
- **Gell-Mann MI for d≥10: too slow (9801 correlations). Use direct MPS contraction instead (Sprint 043)**

## Key Tools Available
- Exact diag: n≤10
- DMRG (TeNPy): n>10, 1D systems only
- Correlation-function MI reconstruction: all-pairs MI from MPS (d=2-5 tested, Gell-Mann)
- **Direct MPS contraction MI: all-pairs MI from MPS for ANY d (Sprint 043, d=10-20 tested, 1.4s for n=8)**
- Custom PottsChain model: Kronecker-delta coupling for any q (Sprint 042)
- IBM QPU: 580s remaining, ibm_kingston (Heron, 156 qubits)
- SU(d) generalized Gell-Mann basis: implemented for d=2-5
