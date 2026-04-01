# Sprint 032 — Entanglement Hamiltonian Locality: Testing Bisognano-Wichmann

**Date:** 2026-04-01
**Status:** In Progress

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
