# Sprint 043 — First-Order Search: Large q Potts (q=10, q=20)

## Question
Does 1D quantum Potts EVER become first-order at large q? The 2D classical Potts model transitions from second-order to first-order at q>4. Sprint 042 showed q=5 is still second-order in 1D quantum. If q=10 or q=20 still show MI-CV crossings, 1D quantum Potts with transverse field may NEVER be first-order — a fundamental difference from the classical case.

## Literature Context
- 2D classical Potts: q≤4 second-order, q>4 first-order (Baxter, exact)
- 1D classical Potts: no phase transition at finite T (Ising argument extends)
- 1D quantum Potts with transverse field: maps to anisotropic 2D classical via Suzuki-Trotter, but anisotropy can change transition order
- Fractal lattice studies (PRE 96, 062105, 2018): transition order depends on geometry and connectivity, not just q

## Prediction
q=10 and q=20 will STILL show crossing curves (second-order). The 1D quantum transverse-field Potts model preserves second-order character for all q because the extreme anisotropy of the quantum-classical mapping suppresses the entropic mechanism that drives first-order transitions in 2D.

## What would be surprising
- Step function (no crossings) at q=10 or q=20 → 1D quantum Potts HAS a first-order boundary
- Crossings that weaken dramatically with q → borderline/marginal behavior

## Experiments

### exp_043a — q=10 Potts timing test (single point)
- **Goal:** Verify q=10 PottsChain runs within 60s at n=8
- **Result:** [pending]

### exp_043b — q=10 Potts n=8 g-sweep
- **Goal:** Find transition region, measure CV profile
- **Result:** [pending]

### exp_043c — q=10 Potts n=12 crossing test
- **Goal:** Check if CV curves cross (second-order) or step (first-order)
- **Result:** [pending]

### exp_043d — q=20 Potts n=8 transition test
- **Goal:** Push to extreme q, find transition, check CV behavior
- **Result:** [pending]

## Findings
[pending]
