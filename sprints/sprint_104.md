# Sprint 104 — Energy-Entropy Hierarchy Universality Test: J1-J2 Chain

**Goal:** Test whether the Casimir-energy/entropy decoupling (confirmed novel in S_q Potts, Sprint 098) is universal — does it appear in a completely different model?

**Result: NEGATIVE.** The hierarchy is NOT universal. It is specific to the walking mechanism.

---

## Model

H = J1 Σ S_i·S_{i+1} + J2 Σ S_i·S_{i+2}, spin-1/2, periodic BC.

Three regimes tested:
- **XX (Δ=0, J2=0):** c=1, v=1, no log corrections. Clean CFT control.
- **Heisenberg (Δ=1, J2=0):** c=1, v=π/2, mild log corrections (marginal umklapp operator).
- **J1-J2 BKT (Δ=1, J2=0.2412):** c=1, v≈1.177, strong BKT log³ corrections.

## Literature

- Okamoto & Nomura (1992): J2c/J1 = 0.2412 ± 0.0001
- BKT entropy corrections: c_eff = 1 + A/ln³(N) (arXiv:2410.01683)
- **No prior comparison of Casimir vs entropy convergence at this BKT point.**

---

## Exp 104a — Casimir Energy vs Entropy at Three Regimes

Exact diag in S_z=0 sector, periodic BC, N=8-20.

| Model | v (exact/gap) | c_eff(18,20) | c_Cas(18,20) | |Δc_eff| | |Δc_Cas| | Winner |
|-------|--------------|-------------|-------------|---------|---------|--------|
| XX | 1.000 (exact) | 0.994 | 1.006 | 0.64% | 0.65% | tied |
| Heisenberg | π/2 (exact) | 0.999 | 1.014 | 0.07% | 1.38% | entropy 21× |
| BKT | 1.177 (gap) | 0.996 | 1.000 | 0.37% | 0.03% | Casimir 14× |

**Surprise:** Hierarchy REVERSES between models. Entropy wins at Heisenberg, Casimir wins at BKT.

## Exp 104b — Velocity Extraction and Multi-Size Fits

Independent velocity from S_z=0→S_z=1 gap (scaling dim x=1/(4K)):
- XX: v_gap = 0.501·2 = 1.002 (K=1, x=1/4 → v=2·gap·N/π). Matches v=1. ✓
- Heisenberg: v_gap = 1.387 (K=1/2, x=1/2 → v=gap·N/π). Off 12% from π/2 due to log corrections.
- BKT: v_gap = 1.177 (K=1/2, x=1/2). Smoothly converging, log corrections minimal at BKT point.

Multi-size fit E₀/N = a + b/N² + d/N⁴:
- XX: c_Cas(fit) = 0.9998 (0.02% off)
- Heisenberg: c_Cas(fit) = 1.0052 (0.5% off — log corrections contaminate fit)
- BKT: c_Cas(fit) = 0.9985 (0.15% off)

## Exp 104c — Cross-Model Magnitude Comparison

**The decisive test: how big are the deviations?**

| | Max |Δc_eff| | Max |Δc_Cas| | Hierarchy |
|---|---|---|---|
| **Potts S_q walking** (q=7-8) | **26%** | 1.4% | Casimir 16× better |
| **J1-J2 chain** (all regimes) | **0.6%** | 1.4% | Model-dependent |

**Potts entropy deviations are 41× LARGER than J1-J2.**

Heisenberg c_Cas deviation follows 1/ln(N): fit c_Cas - 1 = 0.283/ln(N) − 0.084. Confirmed logarithmic.

## Why the Hierarchy Differs

Three distinct mechanisms:

1. **Marginal operator logs (Heisenberg):** The umklapp operator generates 1/ln(N) corrections to the Casimir 1/N² coefficient. These are large (1.4% at N=20) and slow to decay. Entropy formula S = c/3·ln(N) absorbs the log correction into the non-universal constant c', so c_eff converges faster. → **Entropy wins.**

2. **BKT log³ corrections (J1-J2 at J2c):** At the BKT critical point, the leading correction to entropy goes as 1/ln³(N) ≈ 4% at N=20. Casimir energy corrections are suppressed because the marginal operator is exactly at the boundary of relevance. → **Casimir wins.**

3. **Walking (Potts q>4):** Entanglement spectrum reorganization — weight flows from tail into (q-1)-fold multiplet. This creates O(1) entropy deviations (up to 26%) while energy observables track Re(c) via lowest levels only. → **Casimir wins by O(1) amount.** This is qualitatively different from mechanisms 1 and 2.

## Conclusions

1. **The energy-entropy hierarchy is NOT universal.** Direction depends on the type of finite-size correction.
2. **The Potts walking hierarchy is UNIQUE** in producing O(1) entropy deviations while Casimir stays accurate. The J1-J2 chain shows only O(1%) differences between methods.
3. **Scope of Sprint 098 finding:** "Walking transitions in permutation-symmetric models have energy-entropy decoupling where energy tracks complex CFT while entropy does not." This is NOT a general CFT property.
4. **New insight:** Marginal operator log corrections preferentially pollute Casimir (not entropy), while BKT corrections preferentially pollute entropy (not Casimir). The direction depends on whether corrections multiply the 1/N² Casimir term or add to the ln(N) entropy term.

**This is a strong NEGATIVE result that constrains the scope of our confirmed-novel Casimir finding.** The finding remains novel but is walking-specific, not universal.
