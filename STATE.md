# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 126 — Corrected spectral chi_F with factor-2 fix. Key finding: spectral method has systematic negative alpha bias of 0.005–0.038 because non-dominant fraction grows with size. Exact chi_F (finite-difference) is the gold standard — no truncation artifacts. Extracted authoritative exponents for both S_q and hybrid models at q=3,4,5,7 with proper error bars. Model comparison survives all corrections.

## CRITICAL: Use Exact chi_F Going Forward
Spectral chi_F (even with factor-2 correction) captures only the dominant state. The non-dominant states are at high eigenvalue indices. This causes alpha to be underestimated by 0.005–0.038. Use finite-difference exact chi_F for all future exponent claims.

## CRITICAL: Model Identity Correction (April 2026 Audit)
**All experiments from Sprint 076 onward use the STANDARD S_q Potts model**, not the "Potts-clock hybrid." Sprints 119-122 compare both models. See KNOWLEDGE.md.

## Active Research Thread
**Methodology validated, exponents corrected.** The authoritative alpha values (from exact chi_F) are:
- S_q: q=3→1.48, q=4→1.80, q=5→2.09, q=7→2.61 (all walking, z_m>1)
- Hybrid: q=3→1.48, q=4→1.56, q=5→1.44, q=7→1.03 (continuous for q≥5, z_m<1)
Next: extend to larger sizes with exact chi_F for tighter constraints.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — BLOCKED on credentials (qiskit-ibm.json empty)

## Top 3 Next Experiments
1. **Extend exact chi_F to GPU sizes** — q=5 n=10, q=4 n=11 for tighter alpha with more data points
2. **Refit alpha(q) curve** using exact chi_F values — update logarithmic fit from prior sprints
3. **Hybrid z_m crossing precision** — use exact chi_F to refine q_cross = 3.58±0.04

## What's Been Ruled Out
- Spectral chi_F as primary exponent method — systematic negative bias (Sprint 126)
- "frac=1.000" single-multiplet dominance — actual frac is 0.92-0.99 (Sprint 125)
- k-truncation as source of spectral error — non-dominant states at high eigenvalue indices (Sprint 125-126)
- iDMRG overlap chi_F for S_q Potts — non-abelian symmetry breaks convergence (Sprint 124)
- Hybrid chi_F super-scaling — alpha < 2 for all q>=5
- alpha(q=4)=2.0 for S_q at n<=20 — measured 1.80 periodic (Sprint 126)
- VQE for critical GS — fidelity ~50% at n=4 (Sprint 123b)

## Key Tools Available
- **Exact chi_F** (finite-difference, periodic BC): q=2 n<=18, q=3 n<=14, q=4 n<=11, q=5 n<=10, q=7 n<=8
- Spectral chi_F (dominant state only): same sizes, systematically low by 0.5-4%
- chi_F DMRG (open BC): q=2 n<=24, q=4 n<=20
- **NEW**: hamiltonian_utils.py (shared builders), fss_utils.py (error bars via lmfit)
- GPU eigsh practical limit: ~10M dim
- IBM QPU: 580s remaining (credentials needed)
