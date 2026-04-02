# Sprint 102 — Fidelity Susceptibility Across the Walking Boundary

**Goal:** Extract correlation length exponent ν from fidelity susceptibility (χ_F) scaling at g_c = 1/q for S_q Potts chains with q=2,3,4,5,7. Fidelity susceptibility probes the ground state wavefunction overlap directly — independent of entropy. Key questions:
1. Does χ_F_max scaling yield ν consistent with prior gap-based ν measurements?
2. Does the χ_F peak width scale differently for walking (q=5) vs real CFT (q=2,3)?
3. Can χ_F distinguish continuous from weakly first-order behavior?

**Method:** χ_F(g) = (2/N)·(1 - |⟨ψ(g)|ψ(g+δg)⟩|) / δg² for S_q Potts periodic chains. Scan g near g_c = 1/q. Extract peak height and width scaling with N. FSS: χ_F_max ~ N^α where α = 2/ν - 1 (continuous) or α = 2 (first-order in 1D).

**Literature:** Fidelity susceptibility as quantum phase transition probe: Zanardi & Paunković (PRA 2006), Gu (Int. J. Mod. Phys. B 2010). Not measured for q>4 Potts walking regime. Original plan (entanglement asymmetry, Ares et al. 2023) gives ΔS_A = 0 trivially for Z_q-symmetric ground states.

---

## Experiment 102a — χ_F scan near g_c

**Setup:** S_q Potts periodic chain at g_c = 1/q. δg = 10⁻⁴. Scan g over window centered at g_c. Optimized: build H_coup and H_field once, then H(g) = H_coup + g·H_field.

**Sizes:** q=2: n=6-12. q=3: n=6-10. q=5: n=6,8. q=7: n=6.

**Results:**

| q | n | dim | g_peak | g_c | shift | χ_F_max | FWHM | time(s) |
|---|---|-----|--------|-----|-------|---------|------|---------|
| 2 | 6 | 64 | 0.4250 | 0.500 | -0.075 | 0.73 | 0.240 | 0.1 |
| 2 | 8 | 256 | 0.4550 | 0.500 | -0.045 | 0.95 | 0.225 | 0.1 |
| 2 | 10 | 1024 | 0.4700 | 0.500 | -0.030 | 1.19 | 0.210 | 0.1 |
| 2 | 12 | 4096 | 0.4850 | 0.500 | -0.015 | 1.43 | 0.180 | 0.6 |
| 3 | 6 | 729 | 0.3111 | 0.333 | -0.022 | 4.09 | 0.133 | 0.1 |
| 3 | 8 | 6561 | 0.3222 | 0.333 | -0.011 | 6.12 | 0.089 | 0.9 |
| 3 | 10 | 59049 | 0.3222 | 0.333 | -0.011 | 8.27 | 0.078 | 45.7 |
| 5 | 6 | 15625 | 0.2000 | 0.200 | 0.000 | 36.29 | 0.038 | 6.7 |
| 5 | 8 | 390625 | 0.2000 | 0.200 | 0.000 | 66.18 | 0.023 | 128.4 |
| 7 | 6 | 117649 | 0.1429 | 0.143 | 0.000 | 166.88 | 0.017 | 52.2 |

**Observations:**
- Peak location converges to g_c from below for q=2,3 (finite-size shift). For q=5,7 the peak is EXACTLY at g_c — no finite-size shift.
- χ_F at n=6 grows dramatically with q: 0.73 → 4.09 → 36.29 → 166.88

---

## Experiment 102b — FSS analysis + q=2 n=14

**Additional data:** q=2 n=14 (dim=16384): g_peak=0.4875, χ_peak=1.67, FWHM=0.163.

**FSS: χ_F_max ~ N^α → ν = 2/(α+1):**

| q | ν_exact | α_measured | ν_fidelity | ratio |
|---|---------|------------|------------|-------|
| 2 | 1.000 | 0.981 | 1.009 | 1.009 |
| 3 | 0.833 | 1.379 | 0.841 | 1.009 |
| 5 | 0.83 (gap) | 2.089 | 0.648 | 0.781 |

**Pairwise ν from consecutive sizes:**
- q=2 (10,12): 0.996, (12,14): 0.993 — converging to 1.0 ✓
- q=3 (6,8): 0.833 — exact at first pair!
- q=5 (6,8): 0.648 — 22% below gap-based ν=0.83

---

## Experiment 102c — q=4 crossover

**q=4 results:** n=6: χ_peak=13.98, n=8: χ_peak=22.75.
**q=4 FSS:** α=1.69, ν=0.743 (exact: 2/3=0.667, ratio 1.114).

**Complete scaling exponent α(q):**

| q | α_expected (from known ν) | α_measured | deviation |
|---|--------------------------|------------|-----------|
| 2 | 1.00 | 0.98 | -2% |
| 3 | 1.40 | 1.38 | -1% |
| 4 | 2.00 | 1.69 | -16% (log corrections) |
| 5 | 1.41 (if ν=0.83) | 2.09 | **+48%** |

**χ_F at n=6 vs q — exponential growth:**
- Fit: χ_F(n=6) ~ exp(1.06·q - 1.96)
- Growth rate: ~2.88× per unit q

---

## Key Findings

1. **Fidelity susceptibility validates real CFT exponents.** q=2: ν=1.009 (exact 1.0). q=3: ν=0.841 (exact 0.833). Agreement <1%.

2. **Walking regime gives anomalous χ_F scaling.** q=5: α=2.09, giving ν_eff=0.648 instead of gap-derived 0.83. This is NOT a small FSS correction — the discrepancy is 22%, far beyond the 1% seen at q=2,3.

3. **q=5 scaling exponent α≈2 matches first-order prediction.** In 1D, first-order transitions give χ_F ~ N² (α=2). The measured α=2.09 at q=5 is consistent with incipient first-order character in the walking regime.

4. **χ_F grows exponentially with q.** At fixed n=6: ~2.88× per unit q. The ground state wavefunction becomes extremely sensitive to coupling near g_c as q increases.

5. **No finite-size shift for q≥5.** Peak location sits exactly at g_c=1/q. For q=2,3, the peak approaches g_c as ~1/N. This zero shift is consistent with the self-duality fixing g_c exactly.

6. **Observable-dependent exponents in walking regime.** This is the THIRD instance: entropy sees c_eff≠Re(c), Casimir energy tracks Re(c) exactly, and now fidelity susceptibility gives ν_eff≠ν_gap. Walking makes different observables "see" different critical exponents.

**POTENTIALLY NOVEL:** First fidelity susceptibility measurement across the walking boundary in S_q Potts chains. Discovery that χ_F scaling exponent shifts from continuous (α≈2/ν-1) to first-order-like (α≈2) at the walking threshold. Another instance of observable-dependent apparent exponents in walking regime.

---

## Surprises
- Entanglement asymmetry is trivially zero for symmetric ground states (pivoted to fidelity)
- q=5 α≈2 (first-order) despite gap/correlators showing continuous behavior
- Zero finite-size peak shift for q≥5 vs systematic shift for q≤3
- χ_F at n=6 spans 228× range across q=2-7

## Files
- exp_102a_fidelity_suscept.py — χ_F scan for q=2,3,5,7
- exp_102b_fidelity_fss.py — FSS analysis + q=2 n=14
- exp_102c_fidelity_q4_crossover.py — q=4 crossover
- results/sprint_102a_fidelity_suscept.json
- results/sprint_102b_fidelity_fss.json
- results/sprint_102c_fidelity_q4.json
