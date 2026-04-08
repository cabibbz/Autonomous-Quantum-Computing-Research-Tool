# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 118 -- q=4 chi_F extended to 8 sizes (n=4-11). **alpha(q=4) = 1.771 +/- 0.001 (converged).** Differs from exact q=4 prediction alpha=2.0 by 11.5%. Likely due to **multiplicative logarithmic corrections at the marginal Ashkin-Teller point** (q=4 is exactly at the second-order/first-order boundary and has known log corrections that modify effective exponents at finite size).

## ⚠ CRITICAL: Model Identity Correction (April 2026 Audit)
**All experiments from Sprint 076 onward use the STANDARD S_q Potts model**, not the "Potts-clock hybrid." Code audit confirmed: all scripts use `for k in range(1, q)` (field = Σ X^k), not X+X†. The model is NOT novel — it is the same model studied by Gorbenko-Rychkov-Zan, Ma & He, Tang et al. See KNOWLEDGE.md for full details. The novel contributions are the PROBES (chi_F spectral decomposition, systematic chi_F scaling), not the model itself.

## Active Research Thread
**S_q Potts walking regime: chi_F findings need reframing.** alpha(q) data for q=5-30 (11 pts) is new. Spectral decomposition/selection rule mechanism (Sprint 107) is the strongest novel result. Casimir vs entropy comparison is useful but expected. Model identity confusion in old knowledge files being corrected.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s -- UNUSED FOR 93 SPRINTS

## Top 3 Next Experiments
1. **Validate alpha(q) at larger sizes** -- Use QMC or transfer matrix methods (quantum-classical mapping) to verify chi_F exponents at n=20-50 for q=5-7. Current 3-size fits at large q are unreliable.
2. **Derive logarithmic alpha(q) theoretically** -- The spectral decomposition (alpha = beta_me + 2*z_m - 1) is empirically strong but needs a theory for why z_m(q) and beta_me(q) are logarithmic.
3. **Revisit hybrid model** -- The Potts-clock hybrid (Sprints 033-075) IS potentially novel (continuous transitions for q>4, no floating phase). Run chi_F spectral decomposition on the HYBRID to compare with S_q results.

## What's Been Ruled Out
- alpha(q=4) = 2.0 -- measured 1.77, likely log corrections at marginal point (Sprint 118)
- alpha(q) quadratic, linear, sqrt, power-law -- dAIC>=7.7 (Sprint 117)
- DMRG open-BC chi_F -- boundary effects (Sprint 113)
- "Our model is novel" for sprints 076+ -- it's the standard S_q Potts (April 2026 audit)
- Plus all prior ruled-out items

## Key Tools Available
- chi_F spectral (periodic BC): q=2 n<=18, q=3 n<=14, q=4 n<=11, q=5 n<=10, q=6 n<=9, q=7 n<=8, q=8-9 n<=7, q=10 n<=7, q=12 n<=6, q=15 n<=5, q=20 n<=5, q=25 n<=5, q=30 n<=4
- GPU eigsh practical limit: ~10M dim
- IBM QPU: 580s remaining
