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

## MI-CV as Phase Transition Order Parameter (Sprints 030, 036, 037)
MI uniformity coefficient of variation classifies transition TYPE by curve shape:

| Transition | MI-CV shape | Crossing? | Gradient scaling |
|-----------|-------------|-----------|------------------|
| Second-order (Ising) | Smooth inflection | Yes (curves cross) | ~n^1.1 |
| First-order | Step function 0→finite | No | ~n^0.9 |
| BKT (infinite-order) | Smooth dome | No | dome narrows |

**Confirmed as genuine order parameter (Sprint 036):** Tested at n=8,16,32,50.
- TFIM transition slope diverges as ~n^1.1 (slope 1.7 → 3.9 → 7.8)
- Ordered phase: CV → 0 with n (uniform MI). Disordered: CV → ∞ with n.
- XXZ BKT dome narrows: growth rate ratio 1.34 at Δ=1 vs 1.17 at Δ=1.5 (n=16→32)
- XXZ Néel phase (Δ=2): CV decreases 1.13 → 0.78 (n=8→32)

**First-order transition (Sprint 037):** FM phase has CV=0.000 exactly. Jump at Δ=-1 grows as ~n^1.0. Presence/absence of curve crossings distinguishes transition orders.

**Critical exponent (Sprint 037→038):** Crossing-point fit gives ν=0.755 at n=8-32, but data collapse (Sprint 038) proves this is a finite-size artifact. Optimal collapse ν converges: 0.80 (all sizes) → 1.04 (n≥16) → 1.12 (n≥24). Ising ν=1 is confirmed. The collapse quality landscape is extremely flat near ν=1 (0.3% difference from optimal), making crossing-point extraction fragile.

**Potts q=3 MI-CV (Sprints 038-039):** 1D 3-state Potts shows same crossing signature as TFIM at n=8,12,16. Crossing point drifts: g_c=0.923 (8,12) → 0.955 (12,16), converging toward self-dual g=1.0.

**MI-CV distinguishes universality classes (Sprint 039).** Joint-optimized data collapse (ν and g_c free) at n=8-16: Potts ν=5/6 gives 14% better collapse than Ising ν=1 (quality ratio 0.86). Optimal ν converges from 0.81 (all sizes) toward 0.87 (large sizes), approaching 5/6≈0.833. CAUTION: Fixed g_c=1.0 gives misleading answer (Ising wins) due to finite-size g_c shift — must jointly optimize.

**Slope exponent as second discriminator:** Potts slope ~ n^1.36, TFIM slope ~ n^1.1. Consistent with 1/ν scaling (Potts 1/ν=1.2 > Ising 1/ν=1.0), with finite-size inflation.

**Potts q=4 MI-CV (Sprint 040):** q=4 (marginal point where 2D Potts transitions from second-order to first-order) shows crossing curves at n=8,12 — same second-order signature as q=3 and TFIM. Crossing at g_c≈0.893 (further below self-dual than q=3's 0.923). Slope at g=1.0 is LOWER than q=3 (1.72 vs 2.27 at n=8). q=4 CV is systematically lower than q=3 above transition (ratio 0.87-0.95) — larger d distributes correlations more evenly. Marginal/BKT character not visible at n≤12.

**q=5 clock MI-CV (Sprint 041):** q=5 clock model STILL shows crossing curves at n=8,12, disproving prediction that crossings would vanish. Crossing at g_c≈0.673 — a 10x larger shift than q=3→4 (0.220 vs 0.030). Slope at g=1.0 halves vs q=4 (0.86 vs 1.72). CV systematically lower than q=4 above transition (0.521 vs 0.727 at g=1.0 n=8).

**g_c scaling law (Sprint 044).** For 1D quantum Potts: g_c = 1.0 for q=2,3 (self-duality protected), and g_c ≈ 0.87*(q-3)^(-0.85) for q≥4. Verified by blind prediction of g_c(7)=0.263 vs measured 0.259 (1.6% error). Full data: q=2→1.0, q=3→1.0, q=4→0.893, q=5→0.41, q=7→0.259, q=10→0.246. The pole at q=3 reflects self-duality breaking. Exponent 0.85 ≈ 5/6 (q=3 Potts ν). Large-q regime flattening begins at q≈7. χ≥20 required for d≥7 (χ=10 gives 25% CV inflation). For clock model: 0.93 (q=2) → 0.923 (q=3) → 0.893 (q=4) → 0.673 (q=5).

**1D quantum Potts is NEVER first-order (Sprint 043).** Tested q=10 (crossing confirmed at g_c≈0.246) and q=20 (identical to q=10 at χ=10, continuous half-chain entropy). All tested q (2, 3, 4, 5, 10, 20) show second-order MI-CV crossings. At q≥10, ground states converge to a universal large-q regime where only the {|0⟩, |1⟩, |q-1⟩} subspace is active. Physical mechanism: the extreme anisotropy of the 1D quantum→2D classical mapping suppresses the entropic mechanism that drives first-order transitions in 2D.

**Clock ≠ Potts for q≥4 (Sprints 041-042).** TeNPy's ClockChain uses cos(2π(s_i-s_j)/q) coupling, which equals Potts δ(s_i,s_j) only for q=2,3. For q≥5, models differ: Clock g_c=0.67 vs Potts g_c=0.41, Potts slope 5.7x steeper. Custom PottsChain model built (Sprint 042) with projector coupling. Both show second-order crossings — the 2D classical "q>4 → first-order" does NOT apply to 1D quantum Potts with transverse field. Anisotropic quantum-classical correspondence preserves second-order character.

**Technique: All-pairs MI.** Two methods:
1. Gell-Mann correlation reconstruction: exact for d≤5 (validated: diff=0 at n=8). Extended to d=3 (Sprint 038), d=4 (Sprint 040), d=5 (Sprint 041, ~16.5s/point). Scales as O(d⁴) correlation calls — impractical for d≥10.
2. **Direct MPS tensor contraction (Sprint 043):** computes ρ_ij by contracting MPS tensors directly. O(n·χ²·d) per pair. 1.4s for 28 pairs at d=10, n=8. Works for ANY d. Replaces Gell-Mann approach for d≥10.

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
- **Potts q=3**: GHZ-3 → Product (like TFIM but with triplet spectrum from Z₃)
- **Potts q=4**: GHZ-4 → Product (crossing signature like q=3, marginal corrections not yet visible)
- **Clock q=5**: GHZ-5 → Product (crossings persist but shifted to g≈0.67, slope halved vs q=4)
- **Potts q=5**: GHZ-5 → Product (crossings at g≈0.41, slope 5.7x steeper than clock, NOT first-order)
- **Potts q=10**: GHZ-10 → Product (crossings at g≈0.25, confirmed second-order)
- **Potts q=20**: GHZ-20 → Product (identical to q=10 at χ=10, continuous entropy, second-order)
