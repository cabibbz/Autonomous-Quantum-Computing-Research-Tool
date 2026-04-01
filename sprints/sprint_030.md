# Sprint 030 — Heisenberg XXZ Model: Is Scale-Free Universal at All Critical Points?

**Date:** 2026-03-31
**Status:** In Progress

## Motivation

Sprint 029 discovered the Scale-Free entanglement archetype at the TFIM critical point. The TFIM has a single, second-order (Ising universality) phase transition. The Heisenberg XXZ model has a *richer* phase diagram with fundamentally different transition types:

- **Δ < -1:** Ferromagnetic (FM) — doubly degenerate ground state (all-up/all-down superposition = GHZ?)
- **-1 < Δ < 1:** XY/Luttinger liquid — gapless, critical *throughout* the entire phase
- **Δ > 1:** Néel antiferromagnetic — gapped, staggered order
- **Δ = -1:** FM transition (first-order)
- **Δ = 1:** Berezinskii-Kosterlitz-Thouless (BKT) transition (infinite-order, topological)

Key questions:
1. Is the entire XY phase Scale-Free (since it's critical throughout)?
2. Does the BKT transition (topological, infinite-order) look different from the Ising transition in archetype space?
3. What archetype is the Néel phase? It has staggered order — neither GHZ (uniform) nor product.
4. Does MI uniformity CV detect the BKT transition (notoriously hard to detect)?

**Literature check:** Entanglement entropy in XXZ chains is well-studied (concurrence divergence at critical points, entanglement spectrum). However, no prior work maps the full phase diagram onto an entanglement archetype framework or tests MI uniformity CV as a phase detector. The BKT transition is known to be hard to detect with local order parameters — our information-theoretic measures may have an advantage.

## Experiments

### 30a: XXZ Ground State Entanglement Sweep
- Exact diagonalization for n=6 (practical for density matrix ops)
- Sweep Δ from -2.0 to 3.0 in 50 steps
- Compute: half-cut entropy, pairwise MI, I3, negativity spectrum
- Map the entanglement landscape across all three phases

### 30b: Archetype Classification at Each Phase
- Classify each Δ value against the 5 known archetypes
- Compute distances in the (MI_avg, I3_avg, negativity_uniformity) space
- Identify what archetype the Néel phase maps to
- Compare FM phase to TFIM ordered phase (both should be GHZ?)

### 30c: MI Uniformity CV as Universal Phase Detector
- Compute MI uniformity CV across the full Δ sweep
- Compare to half-cut entropy and concurrence as phase detectors
- Test whether CV detects the BKT transition at Δ=1
- Compare TFIM (Sprint 029) and XXZ phase diagrams side by side

---

## Results

### 30a: XXZ Ground State Entanglement Sweep (n=6)

**Sweep:** Δ from -2.0 to 3.0, 67 points (extra near transitions), 3.1s total.

**FM phase (Δ < -1): DEGENERATE — unreliable**
- Ground state degeneracy (gap ≈ 0) causes eigsh to return random superpositions
- All values (MI, I3, entropy) fluctuate wildly depending on phase of superposition
- True ground states are product states |↑↑...↑⟩ and |↓↓...↓⟩ (zero entanglement)
- A 50/50 superposition gives GHZ (S=1, MI=1, I3=1)
- Need Sz-sector decomposition for proper analysis

**FM transition (Δ = -1): ENTROPY SPIKE**
- Half-cut entropy peaks at 1.647 (highest value in entire sweep!)
- This is the first-order transition — discontinuous jump in entanglement
- avg_MI = 0.389, avg_I3 = 0.059 — MI drops while entropy spikes

**XY phase (-1 < Δ < 1): CRITICAL THROUGHOUT**
- Half-cut entropy: 1.113 ± 0.119 (high, slowly decreasing with Δ)
- Avg MI: 0.469 ± 0.087 (moderate, non-uniform)
- Avg I3: 0.118 ± 0.129 (positive → crosses zero at Δ ≈ 0.7)
- Negativity: 0.754 (high)
- MI uniformity CV: 0.791 ± 0.479
- Gap: non-zero at finite n (gapless in thermodynamic limit)

**BKT transition (Δ = 1): INVISIBLE AT THIS SCALE**
- No discontinuity in ANY measure
- All quantities continuous through Δ = 1
- Gap reaches local maximum near Δ ≈ 1 (0.492) — BKT gap is exponentially small in 1/n
- I3 crosses zero near Δ ≈ 0.7 (before the transition!)

**Néel phase (Δ > 1): SLOWLY EVOLVING**
- Half-cut entropy: 1.032 ± 0.006 (nearly constant!)
- Avg MI: 0.417 ± 0.036 (slowly increasing with Δ)
- Avg I3: 0.043 ± 0.057 (positive, growing)
- MI uniformity CV: 1.108 ± 0.197

**Key observations:**
1. I3 sign change occurs at Δ ≈ 0.7 inside the XY phase, NOT at the BKT transition
2. MI uniformity CV: FM ≈ 0 (uniform), XY ≈ 0.79, Néel ≈ 1.11 — monotonic increase
3. The XY phase is NOT uniformly Scale-Free — entropy and MI vary continuously within it
4. The Néel phase has I3 > 0 and growing — resembles GHZ-like redundant correlations (staggered order = "antiferromagnetic broadcast")

### 30b: XXZ Archetype Classification

**Method:** Sz=0 sector projection for FM phase (fixes degeneracy). 10 representative Δ values. Archetype distance computed from (avg_MI, MI_CV, avg_I3, nn_MI_fraction).

**Archetype Map:**

| Phase | Δ | S_half | avg_MI | MI_CV | avg_I3 | Best Archetype |
|-------|---|--------|--------|-------|--------|----------------|
| FM deep | -2.0 | 1.413 | 0.781 | 0.217 | +0.612 | Democratic (GHZ) |
| FM near | -1.2 | 1.923 | 0.500 | 0.142 | +0.153 | Distributed (W) |
| FM transition | -1.0 | 1.469 | 0.629 | 0.000 | +0.356 | Democratic (GHZ) |
| XY near-FM | -0.5 | 1.085 | 0.492 | 0.541 | +0.154 | Scale-Free |
| XY isotropic | 0.0 | 1.039 | 0.424 | 0.960 | +0.053 | Scale-Free |
| XY near-BKT | 0.5 | 1.029 | 0.395 | 1.211 | +0.009 | Scale-Free |
| XY at 0.9 | 0.9 | 1.026 | 0.387 | 1.286 | -0.004 | Geometric (Cluster) |
| BKT transition | 1.0 | 1.026 | 0.387 | 1.289 | -0.005 | Geometric (Cluster) |
| Néel near | 1.5 | 1.029 | 0.395 | 1.225 | +0.008 | Scale-Free |
| Néel deep | 3.0 | 1.045 | 0.501 | 0.703 | +0.175 | Scale-Free |

**Key Discoveries:**

1. **The XY phase is NOT uniformly Scale-Free** — it evolves continuously from Scale-Free (near FM) through to Geometric/Cluster-like (near BKT). MI_CV increases monotonically from 0.54 to 1.29.

2. **The BKT transition is invisible in archetype space** — Δ=0.9 and Δ=1.0 are nearly identical. No sharp transition. The BKT point is just where the Geometric crossover reaches maximum. This contrasts sharply with the TFIM Ising transition, which was clearly visible.

3. **The FM transition (Δ=-1) IS sharp** — clean Democratic with MI_CV=0.000 (perfectly uniform MI). First-order transition produces an abrupt archetype change.

4. **Deep Néel evolves BACK toward Democratic** — staggered order is "antiferromagnetic GHZ" = (|↑↓↑↓...⟩ + |↓↑↓↑...⟩)/√2. I3 grows positive with Δ.

5. **The XXZ traces a LOOP in archetype space:** Democratic → Scale-Free → Geometric → Scale-Free → Democratic. The phase diagram is a continuous journey through archetype space, unlike TFIM's sharp GHZ → Scale-Free → Product path.

6. **The I3 sign change at Δ ≈ 0.7** (inside XY phase, before BKT) marks where correlations transition from redundant (broadcast) to irreducible (geometric). This is NOT at a thermodynamic phase boundary.
