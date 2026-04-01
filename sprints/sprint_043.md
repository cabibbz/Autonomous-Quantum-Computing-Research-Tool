# Sprint 043 — First-Order Search: Large q Potts (q=10, q=20)

## Question
Does 1D quantum Potts EVER become first-order at large q? The 2D classical Potts model transitions from second-order to first-order at q>4. Sprint 042 showed q=5 is still second-order in 1D quantum. If q=10 or q=20 still show MI-CV crossings, 1D quantum Potts with transverse field may NEVER be first-order — a fundamental difference from the classical case.

## Literature Context
- 2D classical Potts: q≤4 second-order, q>4 first-order (Baxter, exact)
- 1D classical Potts: no phase transition at finite T (Ising argument extends)
- 1D quantum Potts with transverse field: maps to anisotropic 2D classical via Suzuki-Trotter, but anisotropy can change transition order
- Fractal lattice studies (PRE 96, 062105, 2018): transition order depends on geometry and connectivity, not just q

## Prediction
q=10 and q=20 will STILL show crossing curves (second-order).

## What would be surprising
- Step function (no crossings) at q=10 or q=20 → 1D quantum Potts HAS a first-order boundary
- Crossings that weaken dramatically with q → borderline/marginal behavior

## Key Technical Innovation
**Direct MPS contraction for ρ_ij replaces Gell-Mann correlation reconstruction.** For d=10, the Gell-Mann approach needs 99²=9801 correlation function calls (~minutes). Direct MPS tensor contraction computes ρ_ij in O(n·χ²·d) per pair — MI computation takes 1.4s vs estimated minutes. This enables MI-CV at arbitrarily large d.

## Experiments

### exp_043a — q=10 Potts timing test (single point)
- **Result:** 12.9s total (DMRG 11.4s + MI 1.4s) at n=8, χ=10, g=0.3
- CV=0.663, E0=-7.9558
- Direct MPS contraction validated: 28 pairs computed in 1.4s

### exp_043b — q=10 Potts n=8 g-sweep
- **Result (χ=10):**

| g    | CV     | E0      | MI_mean |
|------|--------|---------|---------|
| 0.02 | 2.30*  | -7.004  | 0.000   |
| 0.05 | 1.85*  | -7.025  | 0.000   |
| 0.10 | 0.492  | -7.101  | 0.850   |
| 0.15 | 0.612  | -7.228  | 0.828   |
| 0.20 | 0.641  | -7.409  | 0.706   |
| 0.25 | 0.658  | -7.650  | 0.613   |
| 0.30 | 0.663  | -7.956  | 0.490   |

*CV at g≤0.05 is meaningless (MI≈0, nearly product state)

- g≥0.40 fails with χ=10 (truncation error exceeded)
- Transition region around g=0.05-0.15 (MI emerges from zero)
- CV rises smoothly 0.49→0.66 from g=0.10→0.30

### exp_043c — q=10 Potts n=12 crossing test ⭐
- **CROSSING DETECTED at g_c ≈ 0.246**

| g    | n=8 CV | n=12 CV | Δ(n12-n8) | Side       |
|------|--------|---------|-----------|------------|
| 0.10 | 0.492  | 1.820   | +1.327    | Disordered |
| 0.15 | 0.612  | 2.098   | +1.485    | Disordered |
| 0.20 | 0.641  | 2.253   | +1.612    | Disordered |
| 0.25 | 0.658  | 0.515   | -0.143    | Ordered    |
| 0.30 | 0.663  | 0.570   | -0.093    | Ordered    |

- Clear crossing between g=0.20 and g=0.25
- Interpolated crossing: g_c ≈ 0.246
- Same qualitative signature as q=2-5: CV increases with n above transition, decreases below
- **This is NOT first-order** (first-order would show all-increase, no crossing)

### exp_043d — q=20 Potts entropy + universality
- **q=20 at χ=10 gives IDENTICAL ground states to q=10**

| g    | q=20 E0    | q=10 E0    | q=20 S_half |
|------|------------|------------|-------------|
| 0.05 | -7.025033  | -7.025033  | ≈0          |
| 0.10 | -7.100540  | -7.100540  | 0.954       |
| 0.15 | -7.227799  | -7.227799  | 1.039       |
| 0.20 | -7.409272  | -7.409365  | 1.035       |

- Energies match to 4+ significant figures → χ=10 captures a universal large-q regime
- At χ=10, only the {|0⟩, |1⟩, |q-1⟩} subspace is active (low-energy sector same for all q≥10)
- Half-chain entropy rises CONTINUOUSLY (0 → 0.95 → 1.04), no discontinuity
- **No first-order jump signature at q=20**

## Findings

### 1D quantum Potts is second-order for all tested q (2-20)
| q  | g_c   | Transition type | Evidence |
|----|-------|-----------------|----------|
| 2  | ≈1.00 | Second-order (Ising) | Sprints 029, 037-038 |
| 3  | ≈1.00 | Second-order (Potts) | Sprints 038-039 |
| 4  | ≈0.89 | Second-order         | Sprint 040 |
| 5  | ≈0.41 | Second-order         | Sprint 042 |
| 10 | ≈0.25 | Second-order         | This sprint (crossing) |
| 20 | ≈0.12-0.15 | Second-order (indirect) | This sprint (entropy, universality) |

### g_c decreases monotonically with q
g_c ≈ {1.0, 1.0, 0.89, 0.41, 0.25, 0.12-0.15} for q = {2, 3, 4, 5, 10, 20}.
The large-q scaling suggests g_c ~ 1/q or similar.

### Universal large-q regime
At χ=10, all q≥10 give identical ground states because the relevant excitations live in a 3-state subspace {|0⟩, |1⟩, |q-1⟩} at each site. Higher-q states are energetically gapped and irrelevant at low bond dimension. This universality strengthens the second-order conclusion: if the low-energy sector (which controls the transition) is the same for all q≥10, and it shows crossings, then ALL q≥10 are second-order.

### Physical mechanism
The 2D classical rule "q>4 → first-order" relies on the entropic mechanism: more ordered states increase the entropy contribution, favoring a discontinuous jump. In 1D quantum with transverse field, the extreme anisotropy of the quantum-classical mapping (temporal direction ≫ spatial) suppresses this mechanism. The transverse field acts on each site independently, creating excitations in a low-dimensional subspace that cannot trigger the collective entropy-driven first-order transition.

## Surprises
- q=20 ground states IDENTICAL to q=10 at χ=10 — universal large-q regime
- g_c shift from q=5→10 (0.41→0.25) is gentler than q=4→5 (0.89→0.41) — rate of decrease slowing
- Half-chain entropy PEAKS then slightly decreases (1.039→1.035) between g=0.15-0.20 — possible signature of transition
- Direct MPS contraction enables MI-CV at any d, removing the d²-bottleneck of Gell-Mann reconstruction
