# Sprint 084 — Entanglement Spectrum at the Walking Boundary

## Hypothesis
Walking breakdown is exclusively an entropy phenomenon (Sprints 082-083). But entropy S = -Σ λ_i ln(λ_i) is a single scalar derived from the entanglement spectrum {λ_i}. The spectrum itself contains much more information. Questions:
1. Does the spectrum **shape** change across the walking boundary, or just the overall scale?
2. Is there a specific part of the spectrum (entanglement gap, tail, degeneracy structure) that carries the walking signature?
3. Can we decompose which eigenvalues contribute the entropy deviation?

## Method
Compute ρ_A for half-chain bipartition at g_c=1/q for q=2-8 on periodic S_q Potts chains. Extract full entanglement spectrum {-ln λ_i} (entanglement energies ξ_i), entanglement gap Δξ, per-eigenvalue entropy contributions s_i = -λ_i ln(λ_i), spectral decomposition by degeneracy level.

Sizes: q=2 n=10-14; q=3 n=8-10; q=5 n=6-8; q=6 n=6-7; q=7 n=5-7; q=8 n=5-6.

## Results

### Experiment 084a — Baseline: q=2,3,5

| q | n | Δξ | λ_max | n_sig | top5 %S | deg_1 |
|---|---|-----|-------|-------|---------|-------|
| 2 | 14 | 1.034 | 0.713 | 5 | 99.5% | 1 |
| 3 | 10 | 1.316 | 0.627 | 6 | 93.6% | 2 |
| 5 | 8 | 1.692 | 0.552 | 10 | 84.8% | 4 |

**First excited entanglement level has (q-1)-fold degeneracy for all q.** q=3: pairs at ξ=1.784. q=5: 4-fold at ξ=2.286. This is the S_q permutation symmetry acting on the entanglement spectrum.

### Experiment 084b — Broken walking: q=6,7,8

| q | n | Δξ | λ_max | n_sig | top5 %S | deg_1 |
|---|---|-----|-------|-------|---------|-------|
| 6 | 7 | 1.870 | 0.544 | 12 | 73.6% | 5 |
| 7 | 7 | 1.982 | 0.526 | 14 | 64.3% | 6 |
| 8 | 6 | 2.128 | 0.527 | 16 | 58.7% | 7 |

Same (q-1)-fold degeneracy pattern continues. Entanglement gap INCREASES with q (opposite to what "breakdown" might suggest).

### Experiment 084c — Spectral decomposition analysis

**Key finding: walking breakdown = entropy concentration in level 1.**

Entropy decomposition by entanglement level:

| q | %S(lev 0) | %S(lev 1) | %S(tail) | c_eff/Re(c) |
|---|-----------|-----------|----------|-------------|
| 2 | 33.1% | 47.8% | 19.1% | 1.000 |
| 3 | 27.3% | 55.9% | 16.8% | 1.116 |
| 5 | 22.1% | 62.7% | 15.2% | 1.012 |
| 6 | 21.0% | 65.8% | 13.2% | 0.890 |
| 7 | 19.8% | 66.8% | 13.5% | 0.784 |
| 8 | 19.2% | 69.1% | 11.7% | 0.739 |

As q increases: Level 1 (the (q-1)-fold degenerate multiplet) absorbs progressively MORE entropy: 48% → 69%. This comes at the expense of both the ground state (33% → 19%) and the tail (19% → 12%).

**Entanglement gap scales as Δξ ≈ 0.47 + 0.78·ln(q)** (R²=0.996). The gap INCREASES with q — walking breakdown is not about closing the entanglement gap.

**Normalized spectrum R₂₁ (second-level gap / first gap) decreases monotonically:** 4.31 (q=2) → 2.42 (q=8). The spectrum compresses above the first excited level.

**Tail entropy correlates with c_eff/Re(c)** (Pearson r = 0.80). The higher entanglement levels carry the walking signature.

## Interpretation

The entanglement spectrum provides the microscopic mechanism for walking breakdown:

1. **Why only entropy breaks:** The energy gap, correlators, and Casimir energy depend on the lowest-lying entanglement levels (ground state + first excited), which are perfectly conformal for ALL q. The entropy deviates because it sums over ALL levels, and the tail distribution changes.

2. **What changes:** The (q-1)-fold first excited multiplet absorbs more entropy weight as q grows. With more degenerate states in level 1, the entropy gets "trapped" there instead of spreading to higher levels. Since c_eff comes from the total entropy scaling, this redistribution breaks the CFT prediction.

3. **Why q=5 is special:** At q=5, level 1 has exactly 4 states absorbing 62.7% of entropy — just enough that the distribution still matches CFT scaling. At q≥6, the 5+ states absorb too much (65.8%+), and the tail deficit drives c_eff below Re(c).

## Surprises
- Entanglement gap INCREASES with q — not a gap-closing phenomenon
- The (q-1)-fold degeneracy of entanglement levels matches the energy spectrum degeneracy exactly
- Level 1 entropy fraction monotonically increases from 48% to 69% across q=2-8
- Normalized R₂₁ ratio is a clean walking discriminator: 2.60 (q=5) → 2.42 (q=8)
- Tail entropy fraction has r=0.80 correlation with c_eff/Re(c)

**POTENTIALLY NOVEL:** First entanglement spectrum decomposition across the walking boundary for S_q Potts chain. First identification of entropy concentration in (q-1)-fold multiplet as the microscopic mechanism for walking breakdown. Literature search found no prior entanglement spectrum analysis at q>4.
