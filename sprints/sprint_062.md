# Sprint 062 — Compactification Radius & Twisted Boundary Conditions: What CFT is q>4 Potts?

**Status:** Complete (3 experiments).

## Motivation

We know the q→∞ Potts CFT approaches a free boson (Sprint 061: harmonic ratios → k² as O(q⁻³)). But at finite q, c > 1 for q≥5, ruling out a single compact boson (c=1). What is the effective CFT?

Key questions:
1. If we extract a compactification radius R(q) from x₁, does R² ~ ln(q) explain c(q) ~ ln(q)?
2. Can twisted boundary conditions (Z_q flux insertion) measure the spin stiffness directly?
3. Does the stiffness × L converge (Luttinger liquid test)?

## Prior Data Used
- x₁(q) from Sprints 058, 061: q=2→0.125, q=3→0.133, q=5→0.101, q=7→0.086, q=10→0.083, q=15→0.071, q=20→0.055, q=30→0.038
- c(q) from Sprints 054-056: c(2)=0.5, c(3)=0.8, c(5)~1.1, c(10)~1.4
- g_c(q) formula: g_c = (1/5)√(q-1) + 1/20

## Experiments

### Exp 062a: Compactification Radius from x₁

**Method:** Extract effective radius R(q) under two hypotheses:
- H1 (momentum mode): x₁ = 1/(4R²) → R = 1/(2√x₁)
- H2 (winding mode): x₁ = R²/4 → R = 2√x₁

**Results:**

| q | x₁ | R (H1) | R² (H1) | R (H2) | R² (H2) |
|---|-----|--------|---------|--------|---------|
| 2 | 0.125 | 1.414 | 2.000 | 0.707 | 0.500 |
| 3 | 0.134 | 1.367 | 1.869 | 0.731 | 0.535 |
| 5 | 0.102 | 1.569 | 2.463 | 0.637 | 0.406 |
| 10 | 0.083 | 1.734 | 3.009 | 0.577 | 0.332 |
| 20 | 0.055 | 2.132 | 4.545 | 0.469 | 0.220 |
| 30 | 0.038 | 2.565 | 6.579 | 0.390 | 0.152 |

R² ~ ln(q) fit (H1): R² = 1.58·ln(q) + 0.085, RMS=0.591 — **poor fit**, R² grows faster than ln(q).

**x₁(q) fitting:** Power law x₁ ≈ 0.206·q^(-0.449) gives best fit (RMS=0.0103). Neither 1/ln(q) nor (ln(q-1))^(-0.74) fits well.

**Key finding: c·x₁ ≈ 1/9 for q≥3.**

| q | c | x₁ | c·x₁ |
|---|---|-----|------|
| 2 | 0.500 | 0.125 | 0.062 |
| 3 | 0.800 | 0.134 | 0.107 |
| 4 | 1.000 | 0.117 | 0.117 |
| 5 | 1.100 | 0.102 | 0.112 |
| 7 | 1.300 | 0.086 | 0.112 |
| 10 | 1.400 | 0.083 | 0.116 |

For q≥3: c·x₁ = 0.112 ± 0.005, tantalizingly close to 1/9 ≈ 0.1111. q=2 Ising is an outlier (0.062 = 1/16).

**Conclusion:** Single compact boson FAILS. x₁→0 and c>1 are incompatible with c=1 for any single boson. Multiple effective DOF required.

### Exp 062b: Twisted Boundary Conditions

**Method:** Replace periodic BC with Z_q-twisted BC: δ(s_{n-1}, (s_0+k) mod q). Ground state energy shift ΔE(twist k) probes spin stiffness. For a Luttinger liquid: ΔE(k) ∝ k², and ρ_s·L converges to a constant.

**Results — ΔE(twist)/Δ₁ ratio:**

| q | n=4 | n=6 | n=8 |
|---|-----|-----|-----|
| 2 | 1.000 | 1.000 | 1.000 |
| 3 | 1.007 | 1.011 | 1.015 |
| 4 | 1.128 | 1.119 | 1.123 |
| 5 | 1.229 | 1.136 | 1.073 |
| 7 | 1.152 | 0.847 | — |

**q=2 Ising: ΔE(twist) = Δ₁ EXACTLY** at all sizes. The twist field IS the σ field. Perfect Luttinger liquid.

**q=3: ratio ≈ 1.01, slowly growing.** Close to 1, consistent with near-Luttinger behavior. ΔE(twist 2) = ΔE(twist 1) to machine precision — conjugate twists identical.

**q≥4: NOT a Luttinger liquid.** ΔE(2)/ΔE(1) deviates strongly from 4 (free boson prediction):

| q | n | ΔE(2)/ΔE(1) | % from 4 |
|---|---|-------------|----------|
| 3 | 8 | 1.000 | 75.0% |
| 4 | 8 | 1.419 | 64.5% |
| 5 | 8 | 1.785 | 55.4% |
| 7 | 6 | 2.222 | 44.4% |

ΔE(2)/ΔE(1) grows toward 4 (free boson) as q increases, consistent with q→∞ free boson limit, but convergence is very slow.

**ρ_s·L convergence:**

| q | n=4 | n=6 | n=8 | drift |
|---|-----|-----|-----|-------|
| 2 | 0.161 | 0.160 | 0.160 | 1.0% |
| 3 | 0.342 | 0.339 | 0.337 | 1.3% |
| 4 | 0.544 | 0.532 | 0.529 | 2.7% |
| 5 | 0.736 | 0.686 | 0.657 | 10.7% |
| 7 | 0.966 | 0.750 | — | 22.3% |

q=2,3: ρ_s·L nearly constant — Luttinger liquid confirmed. q≥4: significant drift, especially q=7 (22%) — FSS corrections grow with q.

### Exp 062c: Twist-Energy = Scaling Dimension Test

**Method:** In CFT, the ground state energy shift under twisted BC equals the twist field scaling dimension: ΔE(k)·L/(2πv) = x_{twist,k}. Compare twist ratios ΔE(k)/ΔE(1) to spectrum ratios R(σᵏ)/R(σ) from Sprint 057.

**Results — Twist ratios vs spectrum ratios (largest available n):**

| q | k | n | Twist ratio | Spectrum R | Diff |
|---|---|---|-------------|-----------|------|
| 4 | 2 | 8 | 1.419 | 1.68 | -15.5% |
| 5 | 2 | 8 | 1.785 | 2.41 | -26.0% |
| 7 | 2 | 6 | 2.222 | 3.32 | -33.1% |
| 10 | 2 | 6 | 2.917 | 3.66 | -20.3% |

**Twist ratios converge TOWARD spectrum ratios but slowly.** At n=4→6 for q=10: twist ratio grows from 2.20 to 2.92 (target 3.66). The twist field approaches σᵏ in the thermodynamic limit but has large FSS corrections.

**q=10 complete twist sectors at n=4:**

| k | ΔE(k)/ΔE(1) | Target k² |
|---|-------------|-----------|
| 1 | 1.000 | 1 |
| 2 | 2.201 | 4 |
| 3 | 3.131 | 9 |
| 4 | 3.697 | 16 |
| 5 | 3.886 | 25 |

Twist ratios at n=4 saturate rapidly — k=5 barely exceeds k=4. At n=6 they spread out significantly (k=2: 2.9, k=3: 4.8, k=4: 6.2, k=5: 6.7), approaching spectrum values.

**Conjugate twist symmetry:** twist(k) and twist(q−k) have identical energy to machine precision (Δ < 10⁻¹⁴). Exact Z_q conjugation symmetry confirmed.

**x_twist convergence toward x₁:**

| q | n=4 | n=6 | n=8 | x₁ (best) |
|---|-----|-----|-----|-----------|
| 2 | 0.125 | 0.125 | 0.125 | 0.125 |
| 3 | 0.135 | 0.135 | 0.136 | 0.134 |
| 4 | 0.132 | 0.131 | 0.132 | 0.117 |
| 5 | 0.125 | 0.115 | 0.109 | 0.102 |

q=2: x_twist = x₁ exactly. q=3: x_twist ≈ x₁. q≥4: x_twist > x₁ at small sizes, converging from above.

## Key Findings

1. **c·x₁ ≈ 1/9 for q≥3.** The product c(q)·x₁(q) is approximately constant at 0.112 ± 0.005 for all q=3-10. This constrains the CFT: x₁ = (1/9)/c. If c ~ ln(q), then x₁ ~ 1/(9·ln(q)). q=2 Ising is an outlier (c·x₁ = 1/16).

2. **Single compact boson ruled out.** c>1 for q≥5 is incompatible with any single free boson (always c=1). Multiple DOF required.

3. **q=2,3 are Luttinger liquids; q≥4 are not.** ρ_s·L converges for q=2,3 (drift <1.5%) but shows 10-22% drift for q≥5. Twist ratios ΔE(2)/ΔE(1) far from quadratic at finite sizes.

4. **Twist field approaches σ field in thermodynamic limit.** Twist ratios converge toward spectrum harmonic ratios but with large FSS corrections (15-33% deficit at accessible sizes).

5. **Free boson is approached from below as q→∞.** Both twist ratios and spectrum ratios converge to k² prediction, but slowly.

## Surprises
- c·x₁ ≈ 1/9 for q≥3 — no theoretical prediction for this
- q=2 twist field is EXACTLY σ at all sizes (no FSS correction)
- q=3 twist doublet has ΔE(1) = ΔE(2) to machine precision (conjugate pair)
- ρ_s·L drift grows dramatically with q (1% at q=2, 22% at q=7)
- q=10 twist sector ratios at n=4 saturate near k=4 — decompactification signature

**POTENTIALLY NOVEL:** The c·x₁ ≈ 1/9 relation for q≥3 Potts CFT appears previously unmeasured. Combined with the twisted BC spin stiffness analysis showing Luttinger breakdown at q≥4, this constrains the class of CFTs describing q>4 Potts.

[Results: results/sprint_062a_compactification.json, sprint_062b_twisted_bc.json, sprint_062c_twist_dimension.json]
