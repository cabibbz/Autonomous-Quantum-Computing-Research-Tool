# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 128 -- Filled q=6 gap in alpha(q) curve. S_q q=6 alpha=2.375+/-0.006 (6 sizes, n=4-9, 10M GPU). Hybrid q=6 g_c=0.474, alpha=1.186+/-0.038 (pairwise dropping to 1.05). alpha(q) log fit chi2/dof improved from 22 to 5. Asymptotic extrapolation validated at q=3 (alpha_inf=1.412 vs exact 1.400).

## CRITICAL: Use Exact chi_F Going Forward
Spectral chi_F has systematic negative alpha bias of 0.005-0.038 (Sprint 126). Use finite-difference exact chi_F for all exponent claims.

## CRITICAL: Model Identity Correction (April 2026 Audit)
**All experiments from Sprint 076 onward use the STANDARD S_q Potts model**, not the "Potts-clock hybrid." Sprints 119-122 compare both models. See KNOWLEDGE.md.

## Active Research Thread
**chi_F exponents across q.** Authoritative alpha values (exact chi_F, Sprints 127-128):

| q | Hybrid alpha | S_q alpha | # sizes (S_q) | Pairwise drift (S_q) |
|---|-------------|-----------|---------------|---------------------|
| 3 | 1.481+/-0.014 | 1.468+/-0.012 | 6 (n=4-14) | Decreasing |
| **4** | **1.549+/-0.012** | **1.794+/-0.011** | 6 (n=4-11) | Oscillating ~1.79 |
| 5 | 1.352+/-0.043 | **2.139+/-0.019** | 5 (n=4-10) | Increasing |
| **6** | **1.186+/-0.038** | **2.375+/-0.006** | 6 (n=4-9) | Increasing |
| 7 | 0.971+/-0.058 | 2.584+/-0.015 | 4 (n=4-8) | Noisy |

alpha(q) S_q = 1.337*ln(q) - 0.023, chi2/dof=5.0

Hybrid alpha(q) peaks at q~4, drops to ~1 at q=7. Extrapolation gives alpha_inf=0 for q>=5 (sub-power-law).

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s -- BLOCKED on credentials (qiskit-ibm.json empty)

## Top 3 Next Experiments
1. **Add q=8 S_q** -- next gap in alpha(q), test quadratic-in-ln(q) vs pure log. q=8 n=7 dim=2.1M (GPU)
2. **S_q q=4 log corrections** -- fit chi_F = A*N^2*(ln N)^{-p} directly (Salas-Sokal prediction p=3/2)
3. **Hybrid q=5,6 sub-power-law test** -- fit chi_F ~ (ln N)^beta to check if growth is logarithmic

## What's Been Ruled Out
- Spectral chi_F as primary exponent method -- systematic negative bias (Sprint 126)
- alpha(q) = 1.86*ln(q) - 0.96 -- was biased by spectral method (Sprint 127)
- S_q q=4 alpha=2.0 at accessible sizes -- oscillating around 1.79 (Sprint 128)
- iDMRG overlap chi_F for S_q Potts -- non-abelian symmetry breaks convergence (Sprint 124)
- VQE for critical GS -- fidelity ~50% at n=4 (Sprint 123b)
- Hybrid q>=5 alpha extrapolation -- diverges to zero, need alternative model (Sprint 128)

## Key Tools Available
- **Exact chi_F** (finite-difference, periodic BC): q=3 n<=14, q=4 n<=11, q=5 n<=10, q=6 n<=9, q=7 n<=8
- chi_F DMRG (open BC): q=2 n<=24, q=4 n<=20
- hamiltonian_utils.py, fss_utils.py, gpu_utils.py (auto GPU for dim>50k)
- GPU eigsh practical limit: ~10M dim (~55s per eigsh call at 10M)
- IBM QPU: 580s remaining (credentials needed)
