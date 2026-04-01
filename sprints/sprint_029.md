# Sprint 029 — Quantum Simulation: Transverse-Field Ising Model Phase Transition

**Date:** 2026-03-31
**Status:** In Progress
**Theme:** New frontier — quantum simulation of condensed matter physics

## Motivation

After comprehensively mapping small-scale QEC (Sprints 014-028), we pivot to quantum simulation. The transverse-field Ising model (TFIM) is the simplest model with a quantum phase transition:

H = -J Σ Z_i Z_{i+1} - h Σ X_i

At h/J = 1 (critical point), the system transitions from ferromagnetic order (h << J) to paramagnetic disorder (h >> J). This is a genuine quantum phase transition driven by quantum fluctuations, not temperature.

**Why this matters:** We've built a rich entanglement analysis toolkit (entropy, MI, I3, entanglement spectrum, archetypes) over 28 sprints. Applying it to a physical model tests whether our measures detect physically meaningful transitions. The TFIM critical point has known exact solutions (Onsager for 2D, Pfaffian for 1D) providing analytic benchmarks.

**Literature check:** arXiv:2602.17662 (Feb 2026) studies VQE+TFIM with entanglement entropy up to 27 qubits. Our gap: nobody has applied the FULL toolkit (MI, I3, entanglement spectrum, archetype classification) to track how the ground state's entanglement *structure* changes across the phase transition. We're not just measuring entropy — we're classifying the entanglement topology.

## Experiments

### 29a: TFIM Ground State Entanglement Across Phase Transition
- Exact diagonalization for n=4,6,8 qubits (open boundary)
- Sweep h/J from 0.0 to 3.0 in 20 steps
- At each point: half-cut entropy, pairwise MI, all-triple I3
- Look for: entropy peak at criticality, MI restructuring, I3 sign changes

### 29b: VQE Optimization Dynamics
- Hardware-efficient ansatz to find TFIM ground states
- Track entanglement measures at each optimization step
- Compare: does VQE discover the correct entanglement structure?
- Focus on h/J = 0.5 (ordered), 1.0 (critical), 2.0 (disordered)

### 29c: Entanglement Spectrum at Criticality
- Full negativity spectrum across all bipartitions
- Compare critical point spectrum to our archetypes (GHZ, W, Cluster)
- Test: does the critical ground state match any known archetype?

---

## Results

### 29a: TFIM Ground State Entanglement — The Ordered Phase IS GHZ

**Setup:** Exact diagonalization, n=4,6,8 qubits, open boundary, h/J swept 0.01 to 3.0 (25 points).

**Key findings:**

1. **Ordered phase (h << J) is literally GHZ:**
   - Half-cut entropy = 1.0 bit (cat state: |↑↑...↑⟩ + |↓↓...↓⟩)/√2
   - Average MI ≈ 1.0 across ALL pairs (uniform maximal correlations)
   - I3 ≈ +1.0 for all triples (redundant/classical-like, exactly like GHZ in Sprint 004)

2. **Disordered phase (h >> J) is product-like:**
   - Half-cut entropy → 0 (approaches |+⟩^n)
   - MI → 0 for all pairs
   - I3 slightly negative (weak irreducible correlations from residual ZZ coupling)

3. **I3 sign change is the sharpest transition marker:**
   - I3 transitions from +1.0 (ordered) through zero to slightly negative (disordered)
   - Sign change occurs near h/J ≈ 1.4-1.5 for all system sizes
   - This is a QUALITATIVE change: redundant → irreducible correlations

4. **For small n, entropy DECREASES through transition:**
   - GHZ-like ordered phase has S=1.0; critical point has S < 1.0
   - In thermodynamic limit, critical entropy S ~ (c/6)log(L) → ∞, would eventually exceed 1.0
   - Crossover to entropy-peak-at-criticality requires n >> 8

5. **Energy gap closes near h/J ≈ 0.5-0.8 (finite-size shifted from h/J=1.0):**
   - Gap minimum shifts toward 1.0 with increasing n (n=8: gap ≈ 0 at h/J ≈ 0.38)
   - Finite-size effects are massive for these small systems

**Surprise:** The TFIM ordered phase is our old friend GHZ — the "classical broadcast in superposition" from Sprint 005. The phase transition is literally GHZ → Product, with I3 sign change as the sharpest qualitative marker.

### 29b: VQE Optimization Dynamics — Ordered Phase is Hardest

**Setup:** 6-qubit HVA (4 layers, 44 params), COBYLA optimizer (200 iters), h/J = 0.5, 1.0, 2.0.

**Key findings:**

1. **VQE difficulty scales with entanglement:**
   - Disordered (h=2.0): 0.28% energy error — product-like state easy to approximate
   - Critical (h=1.0): 2.6% error — harder but manageable
   - Ordered (h=0.5): 10.4% error — GHZ-like long-range correlations are hardest

2. **VQE gets I3 sign WRONG for ordered phase:**
   - VQE avg I3 = -0.04 (negative), exact = +0.73 (positive)
   - VQE has NOT found the GHZ-like structure — it found a locally correlated state instead
   - For critical and disordered phases, VQE gets I3 sign correct

3. **Entanglement trajectory: product → target structure:**
   - Initial |+⟩^n starts at S ≈ 0 (product state)
   - Entanglement builds monotonically during optimization
   - Final entropy: 0.54 (VQE) vs 1.00 (exact) for ordered phase — only half the entanglement built

4. **VQE automatically finds near-correct entanglement for easy phases:**
   - Disordered: S=0.119 (VQE) vs 0.127 (exact) — within 6%
   - Critical: S=0.331 vs 0.467 — 29% deficit
   - Ordered: S=0.539 vs 1.003 — 46% deficit

**Surprise:** VQE's energy-minimization naturally builds the correct entanglement structure for product-like states but FAILS to discover GHZ-like long-range correlations. The HVA ansatz with nearest-neighbor gates builds entanglement locally; GHZ requires global coherence that needs more layers or different topology.
