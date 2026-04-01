# Sprint 063 — Testing c·x₁ ≈ 1/9: Independent Measurement, Clock Model, & Exact Analysis

**Status:** Complete (3 experiments).

## Motivation

Sprint 062 found c(q)·x₁(q) ≈ 0.112 ≈ 1/9 for q=3-10 Potts. Three questions:
1. Does the relation hold at q=15 (independent c measurement)?
2. Is it universal to all Z_q-symmetric models (test on clock model)?
3. Is it an exact identity or a numerical coincidence?

## Experiments

### Exp 063a: c(q=15) via DMRG Entropy Profile

**Method:** DMRG ground state at g_c(15)=0.798 with chi=30, entropy profile fit for c.

**Result:** c(q=15, n=8) = 1.549. With x₁(15) = 0.071 (Sprint 061):
- c·x₁ = 1.549 × 0.071 = **0.110**, within 1% of 1/9 = 0.111.
- Log formula predicted c = 1.606 (+3.7% above measured).

q=20 DMRG exceeded 600s at n=8 — infeasible with current resources.

**The c·x₁ ≈ 0.112 pattern holds at q=15.**

### Exp 063b: Clock Model c·x₁ Test (q=5)

**Method:** Z_5 clock model H = -J·cos(2π(s_i-s_j)/5) - g·(X+X†). Different Hamiltonian, same Z_5 symmetry. Find g_c via energy gap crossing, extract c from DMRG, x₁ from spectrum.

**Clock g_c(q=5) ≈ 0.52 from energy gap method** (n=4,6 crossing at 0.5205). This is much lower than the old MI-CV estimate of 0.67 (Sprint 042), which was a disordered-phase crossover — same error as early Potts sprints.

**Clock q=5 CFT data:**
| Quantity | Clock q=5 | Potts q=5 |
|----------|-----------|-----------|
| g_c | 0.52 | 0.441 |
| c (DMRG) | 1.17 (n=16) | 1.10 |
| \|c/x₁\| | 9.4 | 10.8 |
| x₁ | ~0.124 | 0.102 |
| c·x₁ | **0.146** | **0.112** |

**c·x₁ ≈ 0.146 for clock, NOT 0.112.** The near-constancy of c·x₁ is specific to the Potts Hamiltonian, not universal to Z_q-symmetric models.

Note: c/x₁ appears negative in raw output due to a sign convention in the Casimir formula (inv_diff should use 1/n₂² - 1/n₁² not 1/n₁² - 1/n₂²). The physical value is |c/x₁|.

Also notable: clock q=5 has HIGHER c than Potts q=5 (1.17 vs 1.10) despite same Z_5 symmetry. The cos coupling creates more entanglement than the δ coupling.

### Exp 063c: Exact c·x₁ from Known CFT Data

**Minimal model analysis.** For M(m, m+1): c = 1 - 6/(m(m+1)). The lowest primary is:
- h_{1,2} = (m-2)/(4(m+1)) for m=3 (Ising)
- h_{2,2} = 3/(4m(m+1)) for m≥4 (these are smaller)

**q=3 Potts (M(5,6)) exact values:** c = 4/5, x₁ = 2/15, c·x₁ = **8/75 = 0.10667**. This is NOT 1/9 = 0.1111. The 4% discrepancy confirms c·x₁ is not an exact identity.

**Minimal model c·x₁ formula:** c·x₁ = 3(m+3)(m-2)/(2m²(m+1)²) for m≥4, approaching 0 as m→∞. This is DECREASING — opposite to the Potts pattern which stays approximately constant.

**The Potts c·x₁ ≈ 0.112 relation is approximate, not exact.** Values:
| q | c·x₁ (exact or best) |
|---|---------------------|
| 2 | 1/16 = 0.0625 (outlier) |
| 3 | 8/75 = 0.1067 |
| 4 | 0.117 |
| 5 | 0.112 |
| 7 | 0.112 |
| 10 | 0.116 |
| 15 | 0.110 |

Range for q≥3: 0.107-0.117, spread ±5% around mean 0.112. Weak upward trend for q≥5 (0.112→0.116), slightly below for q=3 (0.107).

## Key Findings

1. **c(q=15) = 1.549 ± ~0.15 (DMRG n=8).** Gives c·x₁ = 0.110, supporting the approximate constancy of c·x₁ up to q=15.

2. **c·x₁ is Potts-specific, not Z_q-universal.** Clock q=5 gives c·x₁ ≈ 0.146, 30% above the Potts value. The near-constancy is a property of the δ coupling, not the symmetry group.

3. **c·x₁ is NOT exactly 1/9.** Exact q=3 value is 8/75 ≈ 0.1067 (4% below 1/9). The approximate constancy is a numerical coincidence of the particular q-dependences of c and x₁ in the Potts model.

4. **Clock g_c(q=5) ≈ 0.52**, correcting old MI-CV estimate of 0.67. Same error pattern as early Potts work.

5. **Clock c(q=5) ≈ 1.17 > Potts c(q=5) ≈ 1.10.** Both models have c > 1 for q=5, confirming the novel CFT regime is not Potts-specific.

## Surprises
- c·x₁ holds to 1% at q=15 (0.110 vs 0.111) despite being approximate
- Clock model has HIGHER c than Potts at same q — cos coupling is "more critical"
- Clock g_c(q=5) = 0.52 << old MI-CV estimate 0.67 — another MI-CV failure
- Minimal model c·x₁ → 0 as m→∞ — opposite to Potts pattern

## Methodological Notes
- q=15 DMRG at chi=30 takes 468s for n=8 — largest q accessible to DMRG entropy profile
- q=20 DMRG infeasible even at n=8
- Clock Hamiltonian sparse construction needs diagonal boundary term (not projectors)

[Results: results/sprint_063a_c_q15_q20.json, sprint_063b_clock_cx1.json, sprint_063c_cx1_analysis.json]
