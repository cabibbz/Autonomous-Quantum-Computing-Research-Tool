# Sprint 032 — Entanglement Hamiltonian Locality: Testing Bisognano-Wichmann

**Date:** 2026-04-01
**Status:** Complete (3/3 experiments)

## Motivation

Sprint 031 revealed that the entanglement spectrum encodes symmetry information invisible to scalar/correlation measures. The natural next question: is the entanglement Hamiltonian H_E = -log(ρ_A) *local*? The Bisognano-Wichmann (BW) theorem from QFT predicts that for ground states of local Hamiltonians, H_E is proportional to the physical Hamiltonian restricted to subsystem A, modulated by a position-dependent "entanglement temperature" β(i).

For 1D critical systems described by CFT, the BW envelope should follow:
β(i) ∝ sin(πi/L_A) — parabolic, hottest at the boundary, coolest in the bulk.

This connects our entanglement spectral analysis to one of the deepest results in quantum field theory. If BW holds, it means entanglement is fundamentally *local* despite being a non-local quantity.

**Literature:** Giudici et al. (2018) showed BW holds for 1D critical chains. Dalmonte et al. (2018) extended to lattice models. Key gap: systematic comparison across phases and models at small scale, and what happens at phase transitions vs gapped phases.

## Experiments

### 32a: TFIM Entanglement Hamiltonian at Criticality
- n=8, subsystem A = left 4 qubits
- Compute H_E = -log(ρ_A) from exact ground state
- Decompose H_E in Pauli basis to identify dominant terms
- Compare to BW prediction: H_BW = Σ β(i) h_i where h_i are TFIM terms in subsystem A
- Measure BW fidelity (normalized overlap)

### 32b: BW Fidelity Across TFIM Phase Diagram
- Sweep h/J from 0.1 to 3.0
- Compute BW fidelity at each point
- Prediction: best at criticality (CFT applies), worse in gapped phases

### 32c: BW Fidelity for XXZ Model
- Sweep Δ from -1.5 to 2.5 (covering FM, XY, Néel phases)
- Test BW universality across different symmetry classes
- Compare to TFIM results

## Results

### 32a: TFIM Entanglement Hamiltonian at Criticality

**System:** n=8, subsystem A = left 4 qubits, h/J = 1.0 (critical)

**Key findings:**
- H_E has 136 non-zero Pauli terms, but is **91% TFIM-type** (excluding identity)
  - TFIM terms (ZZ bonds + X sites): 130.1 norm²
  - Non-TFIM corrections: 13.1 norm² (9%)
  - Identity: 522.5 norm² (trivial overall shift)
- **BW locality confirmed** — the entanglement Hamiltonian is almost entirely composed of physical Hamiltonian terms
- Couplings decrease toward entanglement cut:
  - ZZ bonds: 6.09 (far) → 3.86 → 2.13 (near cut)
  - X sites: 6.05 (far) → 5.24 → 2.88 → 1.06 (near cut)
- **Standard sin(πi/L) envelope is WRONG** — BW fidelity only 0.65 with that form
  - Ratio test CV = 1.0 (ratios vary from 2.7 to 79)
  - But bond ratios 6.09:3.86:2.13 ≈ 3:2:1 → linear in distance from cut
  - This matches the half-infinite BW prediction β ∝ x (distance from cut)
- Largest non-TFIM corrections: XX (weight 2), XZZ (weight 3) — next-nearest-neighbor and 3-body terms
- Weight distribution of residual: 55% weight-1, 38% weight-2, 6% weight-3

**Surprise:** BW locality is extremely strong (91%) even at n=8, but the standard finite-size CFT envelope fails. The correct envelope appears to be approximately linear in distance from cut, not sinusoidal. This may be a finite-size effect or an open-boundary subtlety.

### 32b: BW Fidelity Across TFIM Phase Diagram

**Sweep:** h/J = 0.1 to 3.0, 25 points. Four envelope types tested: linear, sin_inv, uniform, log.

**Key findings:**
- **Locality peaks at h/J ≈ 0.80** (92.1%), NOT at criticality (h/J=1.0, 90.9%)
  - Slight shift into ordered phase — may be finite-size effect
  - Gap at h/J=0.80 is 0.13 (small but nonzero)
- **sin_inv envelope wins 21/25 points** — β ∝ sin(π(n_A - i)/L), consistent with distance-from-cut scaling
- Phase averages:
  - Ordered: locality 85.2%, best BW variance 84.7%
  - Critical: locality 90.6%, best BW variance 90.1%
  - Disordered: locality 78.0%, best BW variance 77.6%
- **BW works in ALL phases** — even deep disordered (h/J=3.0) has 69% locality
- Uniform envelope always worst (68% avg) — position dependence IS essential
- Linear vs sin_inv nearly identical (<0.3% difference) — hard to distinguish at n=8

**Surprise:** BW locality does NOT peak at the critical point where CFT applies. It peaks slightly in the ordered phase (h/J≈0.8) and decreases monotonically into the disordered phase. This suggests BW locality is a general property of ground states of local Hamiltonians, not specific to CFT. The ordered phase has *more* local entanglement structure because ZZ correlations are dominant and match the physical Hamiltonian directly.

### 32c: BW Fidelity for XXZ Model

**Sweep:** Δ = -2.0 to 2.5, 26 points. Three envelope types: linear, sin_inv, uniform.

**Key findings:**
- **XXZ achieves 100.0% BW locality** in XY, BKT, and Néel phases!
  - sin_inv envelope captures 99.4-100.0% of H_E variance
  - This is dramatically better than TFIM's maximum 92%
- **FM phase (Δ < -1) has 20% locality** — ground state is product (S≈0), H_E is trivial
- **FM transition (Δ ≈ -1) shows sharp jump**: locality 20% → 78% → 100% in narrow window
- **U(1) symmetry perfectly preserved**: XX/YY coefficient ratio = 1.0000 throughout XY/Néel phases
  - In FM phase (broken U(1)), ratio → ∞ (YY coefficient = 0)
- **Bond coefficient ratios** are stable across phases: (0-1)/(2-3) ≈ 2.4-2.6, (1-2)/(2-3) ≈ 1.9
  - Consistent with sin_inv envelope, not exactly linear (would give 3:2:1)
- Phase averages:
  - FM: 20.0% locality, 10.5% BW variance (product state, meaningless)
  - FM transition: 46.6% locality, 41.1% BW variance
  - XY: 98.6% locality, 98.4% BW variance
  - BKT: 100.0% locality, 99.9% BW variance
  - Néel: 100.0% locality, 99.7% BW variance

**KEY INSIGHT:** Symmetry constrains the entanglement Hamiltonian. The XXZ model's U(1) symmetry (conserving total S_z) restricts H_E to contain only XX+YY+ZZ terms — exactly the physical Hamiltonian terms. The TFIM's Z₂ symmetry allows additional terms (XX, XZZ, etc.) that appear in H_E but not in H. **Higher symmetry → better BW locality**, because symmetry eliminates non-Hamiltonian terms from H_E.

This explains why BW locality is:
- 100% for XXZ (U(1) continuous symmetry — very constraining)
- 92% for TFIM (Z₂ discrete symmetry — weakly constraining)
- The 8% non-BW content in TFIM H_E consists of Z₂-allowed terms not present in the physical Hamiltonian

## Summary and Connections

Sprint 032 establishes a **fourth level** of entanglement description beyond Sprint 031's three levels:
1. Scalar (entropy) = amount
2. Correlation (MI/I3) = topology
3. Spectral (eigenvalue distribution) = symmetry content
4. **Hamiltonian (H_E structure) = locality and temperature** ← NEW

The entanglement Hamiltonian H_E = -log(ρ_A) is NOT just a mathematical object — it's the physical Hamiltonian with a position-dependent temperature (BW envelope). This temperature gradient is Unruh-like: infinite at the entanglement cut, finite in the bulk. The accuracy of this description is controlled by the symmetry group of the model:
- U(1) continuous symmetry → 100% BW accuracy (XXZ)
- Z₂ discrete symmetry → 92% BW accuracy (TFIM)
- Prediction: S₃ (Potts) should be intermediate

**Connection to prior sprints:**
- Sprint 031's spectral analysis found that eigenvalue structure encodes symmetry. Sprint 032 shows the *eigenvectors* encode locality — the BW envelope.
- Sprint 029-030's archetype classification used MI/I3 (level 2). BW locality (level 4) provides complementary information: it tells you HOW MUCH of the entanglement structure is determined by the physical Hamiltonian vs. emergent quantum correlations.
- The 9% non-BW corrections in TFIM are exactly the terms that create its richer entanglement structure (I3 ≠ 0 for non-adjacent triples).

**Open questions for future sprints:**
- Does the TFIM 9% gap shrink with system size n? (finite-size vs fundamental)
- Potts model: S₃ symmetry is between Z₂ and U(1) — intermediate BW accuracy?
- 2D systems: does the entanglement temperature gradient depend on dimensionality?
- Central charge dependence: does c affect the BW envelope shape?
