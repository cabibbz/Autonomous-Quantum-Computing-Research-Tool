# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 128 -- Five sub-experiments. (a) S_q q=6: alpha=2.375+/-0.006. (b) Hybrid q=6: g_c=0.474, alpha=1.186. (c) Extrapolation fixed: DB-sourced, inf-guarded. S_q q=4 alpha_inf=1.771. (d) q=4 log corrections: Power+1/N^2 wins (dAIC=41 over Salas-Sokal). Salas-Sokal p=3/2 rejected (measured p=0.41). Best alpha=1.760. (e) Hybrid log test: power wins q<=5, log wins q>=6 (dAIC=20). Hybrid thread closed.

## CRITICAL: Use Exact chi_F Going Forward
Spectral chi_F has systematic negative alpha bias of 0.005-0.038 (Sprint 126). Use finite-difference exact chi_F for all exponent claims.

## CRITICAL: Model Identity Correction (April 2026 Audit)
**All experiments from Sprint 076 onward use the STANDARD S_q Potts model**, not the "Potts-clock hybrid." Sprints 119-122 compare both models. See KNOWLEDGE.md.

## Active Research Thread
**S_q chi_F exponents — deepening, not widening.** Authoritative alpha values (exact chi_F):

| q | S_q alpha | Extrap alpha_inf | Pairwise drift | Status |
|---|-----------|------------------|---------------|--------|
| 3 | 1.468+/-0.012 | 1.407+/-0.001 | Decreasing | Converged (exact 7/5) |
| **4** | **1.795+/-0.007** | **1.771+/-0.001** | Decreasing | Key target |
| 5 | 2.094+/-0.002 | 2.093+/-0.015 | Stable | Converged |
| 6 | 2.375+/-0.006 | --- (divergent) | Increasing | Needs more sizes |
| 7 | 2.636+/-0.018 | --- (3 pts) | Increasing | Needs more sizes |

**S_q q=4 is the publishable question.** Power+1/N^2 correction gives alpha=1.760. Log-corrected alpha=2 is WORST fit (dAIC=41). Salas-Sokal p=3/2 rejected. Pairwise converging monotonically to ~1.77. This is new data — no prior chi_F log correction measurement exists.

**Hybrid thread CLOSED (Sprint 128e).** Power-law for q<=4, logarithmic for q>=6, marginal at q=5. Consistent with BKT mechanism at large q. No further compute.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s -- BLOCKED on credentials (qiskit-ibm.json empty)

## Top 3 Next Experiments
1. **S_q q=4 at larger sizes via DMRG** -- extend exact chi_F beyond n=11 to test whether alpha_inf=1.77 holds or eventually converges to 2.0 with log corrections. DMRG open-BC drifts upward (1.505->1.523 at n=20), so periodic BC larger sizes would be definitive.
2. **Add q=8 S_q** -- fill gap in alpha(q), test quadratic-in-ln(q) (chi2/dof=0.06) vs pure log (chi2/dof=25). q=8 n=7 dim=2.1M (GPU)
3. **S_q q=4: test alpha=2 with sub-leading corrections** -- try chi_F = A*N^2*(ln N)^{-p}*(1+c/N^2) (4-param model combining log and power corrections). Needs n>=12 for 4 params — DMRG required.

## What's Been Ruled Out
- Spectral chi_F as primary exponent method -- systematic negative bias (Sprint 126)
- alpha(q) = 1.86*ln(q) - 0.96 -- was biased by spectral method (Sprint 127)
- S_q q=4 alpha=2 with Salas-Sokal log corrections at accessible sizes -- dAIC=41 worse (Sprint 128d)
- S_q q=4 alpha=2 at n<=11 -- pairwise converging to 1.77, not 2.0 (Sprint 128c/d)
- iDMRG overlap chi_F for S_q Potts -- non-abelian symmetry breaks convergence (Sprint 124)
- VQE for critical GS -- fidelity ~50% at n=4 (Sprint 123b)
- **Hybrid chi_F as active thread** -- closed Sprint 128e. Log wins q>=6, power wins q<=4. Story complete.

## Key Literature for Next Sprint
- **Salas & Sokal (1997):** q=4 Potts log corrections, p=3/2 prediction
- **Balog et al. (2007):** Refined log exponents for q=4 Potts
- Our chi_F rejects Salas-Sokal at accessible sizes. Key question: is alpha truly 1.77, or does it cross over to 2.0 at L>>20?

## Key Tools Available
- **Exact chi_F** (periodic BC): q=3 n<=14, q=4 n<=11, q=5 n<=10, q=6 n<=9, q=7 n<=8
- chi_F DMRG (open BC): q=2 n<=24, q=4 n<=20
- hamiltonian_utils.py, fss_utils.py, gpu_utils.py (auto GPU for dim>50k)
- GPU eigsh practical limit: ~10M dim (~55s per eigsh call at 10M)
- IBM QPU: 580s remaining (credentials needed)
