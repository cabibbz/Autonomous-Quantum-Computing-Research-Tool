# Accumulated Knowledge — Edit by topic, not by sprint

## Five Entanglement Archetypes
| Archetype | Example | MI pattern | I3 sign | Negativity | Source |
|-----------|---------|-----------|---------|------------|--------|
| Democratic | GHZ, Ising ordered | Uniform high | +1 (redundant) | Flat | Sprint 005 |
| Distributed | W state | Uniform weak | +0.2 (weak redundant) | Growing with cut | Sprint 005 |
| Geometric | Cluster 1D | Nearest-neighbor only | -1 (irreducible) | Explosive, geometry-dependent | Sprint 007 |
| Topological | Toric code | Non-contractible loops | 0 (zero everywhere) | Uniform floor | Sprint 020 |
| Scale-Free | Critical TFIM | Power-law decay | Mixed | Highest rank, steep decay | Sprint 029 |

## Four Levels of Entanglement Description
Each captures orthogonal information. Different phase transitions are visible at different levels.
1. **Scalar** (entropy) — amount of entanglement
2. **Correlation** (MI, I3) — topology of correlations. I3<0 = irreducible multipartite. I3>0 = redundant/classical-like.
3. **Spectral** (eigenvalue distribution of ρ_A) — symmetry content. U(1) gives doublet degeneracies. Z₃ gives triplets.
4. **Hamiltonian** (H_E = -log ρ_A structure) — locality and entanglement temperature.

## MI-CV as Phase Transition Order Parameter (Sprints 030, 036)
MI uniformity coefficient of variation classifies transition TYPE by curve shape:
- **Discontinuous jump** → first-order transition
- **Monotonic inflection** → Ising / second-order
- **Smooth dome** → BKT / infinite-order

**Confirmed as genuine order parameter (Sprint 036):** Tested at n=8,16,32,50.
- TFIM transition slope diverges as ~n^1.1 (slope 1.7 → 3.9 → 7.8)
- Ordered phase: CV → 0 with n (uniform MI). Disordered: CV → ∞ with n.
- XXZ BKT dome narrows: growth rate ratio 1.34 at Δ=1 vs 1.17 at Δ=1.5 (n=16→32)
- XXZ Néel phase (Δ=2): CV decreases 1.13 → 0.78 (n=8→32)
- Crossing point at h/J ≈ 0.93 (7% finite-size shift from exact h_c=1.0)

**Technique:** All-pairs MI via correlation-function reconstruction of ρ_ij from MPS. Exact for physical states (validated: diff=0 at n=8).

## Archetype Boundaries ≠ Phase Boundaries
I3 sign change occurs at Δ≈0.7 in XXZ, inside the XY phase — not at either thermodynamic transition (Δ=-1 or Δ=1). The entanglement phase diagram has its own topology distinct from thermodynamics.

## Bisognano-Wichmann Locality
H_E ≈ physical Hamiltonian × position-dependent entanglement temperature β(x).
- β(x) follows Unruh-like gradient: hottest at entanglement cut, coldest in bulk.
- BW works in ALL phases, not just at criticality. Peaks in gapped phases.
- Accuracy controlled by H/G-inv ratio (Hamiltonian operator dimension / symmetry-invariant operator dimension).

### H/G-inv Predictor (Sprint 034)
| Model | Symmetry | d | H/G-inv ratio | BW Locality |
|-------|----------|---|---------------|-------------|
| XXZ | U(1) | 2 | 0.143 | 100.0% |
| TFIM | Z₂ | 2 | 0.099 | 91.0% |
| Potts | S₃ | 3 | 0.072 | 76.5% |
| Chiral clock | Z₃ | 3 | 0.040 | 69.7% |
Perfect monotonic correlation. Depends on BOTH symmetry group AND local dimension.

## QEC Arc — Closed (Sprints 014-028)
- [[5,1,3]] is basis-isotropic: confirmed on IBM hardware (Sprint 025, asymmetry 0.040 vs 3-qubit 0.254)
- Small-scale active QEC cannot beat passive encoding (Sprints 026-028)
- Root cause: syndrome extraction (16+ 2Q gates) exceeds code's correction capacity
- Non-FT, flag-FT, repeated measurement — all fail. Threshold theorem requires asymptotic limit.
- First observed QEC advantage: 3-qubit Z-basis Holevo 0.976 vs uncoded 0.959 (Sprint 025)

## Hardware Results (Sprint 025)
Backend: ibm_kingston (Heron, 156 qubits). 20s QPU used.
- Uncoded: avg Holevo 0.967, asymmetry 0.012
- 3-qubit: avg 0.838, asymmetry 0.254
- [[5,1,3]]: avg 0.499, asymmetry 0.040
- Hardware error rates better than conservative noise model (readout 0.35% vs model 1.5%)
- Correlated noise degrades [[5,1,3]] isotropy 4x vs simulator but doesn't destroy it

## Phase Diagram Trajectories
- **TFIM**: Democratic(GHZ) → Scale-Free → Product (one-way)
- **XXZ**: Democratic → Scale-Free → Geometric → Scale-Free → Democratic (loop)
- **Potts**: GHZ-3 → Product (like TFIM but with triplet spectrum from Z₃)
