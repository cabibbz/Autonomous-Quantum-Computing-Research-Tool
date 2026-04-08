# Current State -- Rewrite this completely each sprint

## Last Sprint
Sprint 119 — chi_F spectral decomposition on the **hybrid Potts-clock model** at q=3,5,7,10. **Single-multiplet dominance is universal (frac=1.000 in both models).** Spectral gap symmetry-forbidden in both. Decomposition formula alpha = beta_me + 2*z_m - 1 exact in both. BUT: hybrid alpha **decreases** with q (1.41→0.95→0.55) while S_q alpha **increases** (1.40→2.09→2.65). Walking super-scaling is S_q-specific, not universal. z_m < 1 in hybrid (gap closes slower than 1/N) vs z_m > 1 in S_q.

## ⚠ CRITICAL: Model Identity Correction (April 2026 Audit)
**All experiments from Sprint 076 onward use the STANDARD S_q Potts model**, not the "Potts-clock hybrid." Sprint 119 returns to the hybrid model for comparison. See KNOWLEDGE.md for full details.

## Active Research Thread
**chi_F spectral decomposition as universal diagnostic.** The mechanism is model-independent (selection rule, single-multiplet, exact formula). The exponents (z_m, beta_me) encode model-specific physics: walking (z_m>1, beta_me growing) vs continuous (z_m<1, beta_me→0). This cleanly separates the two universality classes.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s — UNUSED FOR 94 SPRINTS

## Top 3 Next Experiments
1. **Stress-test hybrid g_c accuracy at q=7,10** — rapidly dropping pairwise alpha could be g_c error. Re-verify with chi_F peak scan (small δg sweep).
2. **Hybrid q=4 chi_F** — at the marginal Ashkin-Teller point, both models converge. Test where the alpha divergence starts (q=4 or q=5).
3. **Theoretical understanding** — why does z_m cross 1 at the walking boundary? The spectral decomposition gives a clean observable to match against CFT/complex CFT predictions.

## What's Been Ruled Out
- Hybrid chi_F super-scaling — hybrid alpha < 2 for all q≥5 (Sprint 119)
- alpha(q=4) = 2.0 — measured 1.77 (Sprint 118)
- alpha(q) quadratic, linear, sqrt, power-law for S_q — dAIC>=7.7 (Sprint 117)
- DMRG open-BC chi_F — boundary effects (Sprint 113)
- "Our model is novel" for sprints 076+ — it's the standard S_q Potts (April 2026 audit)

## Key Tools Available
- chi_F spectral (periodic BC): q=2 n<=18, q=3 n<=14, q=4 n<=11, q=5 n<=10, q=6 n<=9, q=7 n<=8, q=8-9 n<=7, q=10 n<=7, q=12 n<=6, q=15 n<=5, q=20 n<=5, q=25 n<=5, q=30 n<=4
- Hybrid chi_F spectral: q=3 n<=12, q=5 n<=10, q=7 n<=8, q=10 n<=7
- GPU eigsh practical limit: ~10M dim
- IBM QPU: 580s remaining
