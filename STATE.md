# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 050 — True Potts Critical Points via Self-Duality & MI-CV Vindication

## Active Research Thread
**Self-duality resolves q=3 Potts g_c: g_c = J/3 = 0.333 exactly.** Derived from Kramers-Wannier duality of our Hamiltonian H = -Jδ - g(X+X†). For q=2,3, X+X† spans all q-1 generators → self-dual. For q≥4, X+X† misses intermediate powers → self-duality broken, g_c unknown.

**MI-CV crossings CONFIRMED at true g_c.** n=8,12 crossing at g≈0.26 (finite-size shifted below g_c=1/3). The crossing signature from Sprints 038-039 is rescued — qualitative conclusions (crossings=2nd order) are correct, just at wrong g values.

**q=4 Potts g_c still unknown.** Pseudo-critical at g≈0.34 (n=8), but self-duality broken. Need n=16+ DMRG (very slow without symmetry conservation at d=4).

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston, 18 circuits)
- Remaining: 580s

## Top 3 Next Experiments
1. **q=3 ν extraction at true g_c** — Data collapse with n=8,12 MI-CV near g_c=1/3. Does ν=5/6 survive? Need n=16 MI-CV (slow DMRG).
2. **q=4,5 g_c via alternative methods** — Use energy gap or order parameter (magnetization) instead of entropy. Exact diag gap at n=6,8 is fast and well-defined.
3. **Hardware test of Potts MI-CV** — q=2 Ising MI-CV can be tested on IBM QPU (n=6-8 qubits). Strongest simulator prediction: crossing at g≈0.20 (Potts convention) = g≈0.80 (TFIChain).

## What's Been Ruled Out
- Small-scale QEC active correction: fails for [[5,1,3]] (Sprints 026-028)
- Pauli fraction as BW metric: fails at large n (Sprint 035)
- BW sin_inv envelope: wrong at finite size (Sprint 035)
- Clock ≡ Potts for q≥4: false (Sprint 041)
- χ=10 for d≥7: NOT converged (Sprint 044)
- Gell-Mann MI for q≥4: unreliable (Sprints 045-046)
- Raw MI-CV for d≥10: dead-pair bias (Sprint 048)
- Entropy FSS for ν extraction: log singularity (Sprint 049)
- **All Sprints 038-048 Potts g_c values: WRONG (Sprint 049)**
- **g_c scaling law g_c(q) = 0.87(q-3)^{-0.85}: INVALIDATED (Sprint 049)**
- **Self-duality for q≥4 with our Hamiltonian: BROKEN (Sprint 050)**
- **Central charge from n≤12 for q=3 Potts: massive overshoot (Sprint 050)**

## Key Tools Available
- Exact diag: n≤10 (n≤8 for q=3, n≤6 for q=4, n≤4 for q=5)
- DMRG (TeNPy): n>10, 1D only. TFIChain has Z₂ conservation. PottsChain has NO conservation (slow).
- Direct MPS contraction MI: all-pairs MI for ANY d (reliable for d≥4, works for d=3 too)
- Custom PottsChain model: Kronecker-delta coupling for any q
- IBM QPU: 580s remaining
- Self-duality: g_c = J/q for q=2 (0.25) and q=3 (0.333) with our Potts Hamiltonian
