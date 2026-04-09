# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 127 -- Extended exact chi_F to GPU sizes (up to 10M-dim Hilbert space). 12 new data points across both models. Error bars reduced 17-62%. S_q q=5 alpha=2.094+/-0.002 is now rock-solid. S_q q=4 alpha continues downward drift (1.795, pairwise 1.779 at largest sizes). Hybrid alpha(q) confirmed non-monotonic.

## CRITICAL: Use Exact chi_F Going Forward
Spectral chi_F has systematic negative alpha bias of 0.005-0.038 (Sprint 126). Use finite-difference exact chi_F for all exponent claims.

## CRITICAL: Model Identity Correction (April 2026 Audit)
**All experiments from Sprint 076 onward use the STANDARD S_q Potts model**, not the "Potts-clock hybrid." Sprints 119-122 compare both models. See KNOWLEDGE.md.

## Active Research Thread
**Exponents with extended data.** Authoritative alpha values (exact chi_F, Sprint 127):
- S_q: q=3->1.468+/-0.012, q=4->1.795+/-0.007, q=5->2.094+/-0.002, q=7->2.636+/-0.018
- Hybrid: q=3->1.468+/-0.012, q=4->1.545+/-0.014, q=5->1.382+/-0.030, q=7->0.957+/-0.067
- alpha(q) log fit slope dropped from 1.86 to 1.306 (prior was biased by spectral method)
- S_q q=4 pairwise alpha at (10,11)=1.779 -- settling below 2.0
- Hybrid alpha(q) non-monotonic, peaks at q~3-4 then drops

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s -- BLOCKED on credentials (qiskit-ibm.json empty)

## Top 3 Next Experiments
1. **Add q=6 S_q and hybrid** -- fill gap in alpha(q) curve. q=6 n=8 dim=1.7M (easy GPU)
2. **Asymptotic extrapolation for S_q q=3,4** -- fit alpha_eff(N) = alpha_inf + c/N^p to get converged alpha
3. **Hybrid z_m at extended sizes** -- extract z_m from spectral decomposition at new GPU sizes

## What's Been Ruled Out
- Spectral chi_F as primary exponent method -- systematic negative bias (Sprint 126)
- alpha(q) = 1.86*ln(q) - 0.96 -- was biased by spectral method, actual slope ~1.31 (Sprint 127)
- S_q q=4 alpha=2.0 at accessible sizes -- pairwise drifting down toward 1.78 (Sprint 127)
- iDMRG overlap chi_F for S_q Potts -- non-abelian symmetry breaks convergence (Sprint 124)
- VQE for critical GS -- fidelity ~50% at n=4 (Sprint 123b)

## Key Tools Available
- **Exact chi_F** (finite-difference, periodic BC): q=3 n<=14, q=4 n<=11, q=5 n<=10, q=7 n<=8
- chi_F DMRG (open BC): q=2 n<=24, q=4 n<=20
- **NEW**: hamiltonian_utils.py, fss_utils.py, gpu_utils.py (auto GPU for dim>50k)
- GPU eigsh practical limit: ~10M dim (~42s per eigsh call)
- IBM QPU: 580s remaining (credentials needed)
