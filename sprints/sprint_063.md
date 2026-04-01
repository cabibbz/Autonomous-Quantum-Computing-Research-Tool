# Sprint 063 — Testing c·x₁ ≈ 1/9: Independent c Measurement at q=15, Clock Model Test

**Status:** In progress.

## Motivation

Sprint 062 found c(q)·x₁(q) ≈ 1/9 ≈ 0.111 for q=3-10 Potts, potentially novel. But:
1. We only tested q≤10 where both c and x₁ have independent measurements.
2. For q=15,20,30, c values were from the log formula (not measured independently).
3. Is c·x₁ ≈ 1/9 universal to Z_q-symmetric models, or specific to Potts?

Key tests:
- Measure c(q=15) independently via DMRG entropy profile → compute c·x₁ with known x₁≈0.071
- Test clock model (different Hamiltonian, same Z_q symmetry) at q=5
- Check if c·x₁ = 1/9 holds exactly for q=3 (c=4/5, x₁=2/15 → c·x₁ = 8/75 ≈ 0.1067, NOT 1/9)

## Experiments

### Exp 063a: c(q=15) via DMRG entropy profile
(results pending)

### Exp 063b: Clock model c·x₁ test at q=5
(results pending)

### Exp 063c: Exact c·x₁ values from known CFT data
(results pending)
