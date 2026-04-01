# Sprint 031 — Entanglement Spectrum: What the Reduced Density Matrix Hides Beyond Entropy

**Date:** 2026-04-01
**Status:** Complete (3/3 experiments)

## Motivation

We've characterized entanglement through scalar measures (entropy, MI, I3, discord, negativity) and discovered 5 archetypes. But the reduced density matrix contains far more information than any single number can capture — its full eigenvalue spectrum (the "entanglement spectrum") encodes edge modes, topological order, and universality classes.

The Li-Haldane conjecture (2008) posits that the entanglement spectrum mirrors the physical edge spectrum of topological phases. A 2025 paper (arXiv 2509.03588) showed this correspondence can *break down* in XXZ chains — entanglement phase transitions can be disconnected from bulk transitions. This directly connects to our Sprint 030 finding that archetype boundaries ≠ phase boundaries (I3 sign change at Δ≈0.7 inside XY phase).

**Key question:** Does the entanglement spectrum reveal structure invisible to scalar entanglement measures? Can it explain the archetype-phase boundary mismatch?

**Literature searched:** arXiv 2509.03588 (entanglement phase transitions in Haldane phase, 2025), Nature Comms 2023 (energy-entanglement spectrum wormhole connection), Quantum 2022 (Li-Haldane conjecture probing).

## Experiments

### 31a: Entanglement Spectrum of 5 Archetypes (n=8)

**Goal:** Compute full eigenvalue spectrum of half-cut RDM for GHZ, W, Cluster 1D, Cluster 2D, and TFIM critical (Scale-Free).

**Results:**

| Archetype | S_vN | Schmidt Rank | NPR | Spectral Gap | Shape |
|-----------|------|-------------|-----|-------------|-------|
| GHZ (Democratic) | 1.000 | 2 | 1.000 | 0.000 | Binary flat (0.5, 0.5) |
| W (Distributed) | 0.967 | 2 | 0.956 | 0.214 | Binary split (0.607, 0.393) |
| Cluster 1D (Geometric) | 1.000 | 2 | 1.000 | 0.000 | Binary flat (0.5, 0.5) |
| Cluster 2D (Topological) | 4.000 | 16 | 1.000 | 0.000 | Maximally flat (all 0.0625) |
| Scale-Free (Critical) | 0.515 | 9 | 0.139 | 0.773 | Steep power-law (α≈12.7) |

**Entanglement energies ξ = -ln(λ):**
- GHZ: [0.693, 0.693] — degenerate
- W: [0.499, 0.934] — small splitting
- Cluster 1D: [0.693, 0.693] — degenerate (identical to GHZ at half-cut!)
- Cluster 2D: [2.773, 2.773, ...×16] — all degenerate at high energy
- Scale-Free: [0.121, 2.179, 7.832, 9.890, 15.831, 17.889, 23.542, 25.600] — structured!

**Surprises:**
- GHZ and Cluster 1D are **spectrally identical** at half-cut — both rank 2 with degenerate eigenvalues. The half-cut can't distinguish Democratic from Geometric topology.
- Scale-Free has the **highest Schmidt rank among 1D states** (9 vs 2) but the **lowest entropy** (0.515 vs 1.0). Most eigenvalues are tiny — the spectrum is steeply decaying.
- Cluster 2D has rank 16 = 2⁴ (maximally mixed RDM) — every degree of freedom across the cut is maximally entangled.
- The entanglement spectrum adds a new discrimination axis: GHZ/Cluster 1D can't be distinguished by entropy, MI, or spectrum alone — you need I3 or discord.

### 31b: Entanglement Spectrum Across TFIM Phase Diagram (n=8)

**Goal:** Track spectrum evolution through ordered → critical → disordered.

**Results:**
- Entropy peak at h/J=0.29 (NOT 1.0 — finite-size shift for n=8)
- Schmidt rank peaks at h/J≈0.36 (rank 10), stays at 10 through critical region, drops in disordered phase
- Entanglement gap minimum near h/J=0.08 (ordered GHZ-like phase)

**CFT Doublet Structure at Criticality (h/J=1.0):**

Entanglement energies (shifted, normalized by first gap):
```
0.000, 1.000, 3.747, 4.747, 7.634, 8.634, 11.381, 12.381
```

Consecutive spacings: **1.0, 2.75, 1.0, 2.88, 1.0, 2.75, 1.0** — alternating!

The spectrum forms **doublets** with internal splitting ~1.0 and inter-doublet spacing ~2.8. This doublet pattern is characteristic of the c=1/2 Ising CFT with two primary fields (identity and energy operator). The alternating spacing encodes the conformal tower structure.

**Phase comparison:**

| Phase | S_vN | Rank | Ent. Gap | NPR |
|-------|------|------|----------|-----|
| Ordered (h<0.5) | 0.857 | 6 | 3.46 | 0.170 |
| Critical (≈1.0) | 0.516 | 10 | 2.06 | 0.125 |
| Disordered (h>2.0) | 0.083 | 6 | 4.59 | 0.130 |

**Surprises:**
- The critical point has the **highest Schmidt rank** despite not having the highest entropy — more eigenvalues active but distributed more unevenly.
- The doublet structure is a spectral fingerprint of the Ising universality class — not visible in any scalar measure.
- The entanglement gap is NON-monotonic through the transition (large → small → large), while entropy is monotonic from ordered to disordered.

### 31c: Entanglement Spectrum Across XXZ Phase Diagram (n=6)

**Goal:** Track spectrum across FM → XY → Néel. Test if I3 sign change at Δ≈0.7 appears in spectrum.

**Results:**

**Universal doublet structure in XY phase:**
Throughout the entire gapless XY phase (Δ=-2 to Δ≈0.95), ALL eigenvalues come in **degenerate pairs** — 4 pairs at rank 8. This is symmetry protection from U(1) conservation of total Sz.

**BKT transition = spectral collapse:**
| Δ | Rank | Degenerate Pairs | S_vN | Ent. Gap |
|---|------|-----------------|------|----------|
| 0.90 | 8 | 4 | 1.313 | 0.000 |
| 0.95 | 6 | 3 | 1.383 | 0.000 |
| 1.00 | 4 | 0 | 1.328 | 0.617 |
| 1.05 | 1 | 0 | 0.000 | ∞ |

The BKT transition causes: (1) degeneracy breaking (4→3→0 pairs), (2) rank collapse (8→6→4→1), (3) gap opening (0→0.617→∞). This happens across a narrow window Δ=0.9-1.05.

**I3 sign change at Δ≈0.7: INVISIBLE in spectrum.**
At Δ=0.6-0.8: rank=8, degeneracies=4, smooth entropy increase (1.11→1.21). The spectrum is completely smooth through the I3 sign change. Archetype boundaries operate on a different level than spectral transitions.

**FM transition (Δ=-1): also invisible.** Rank=8, 4-5 degenerate pairs, entropy≈1.03 throughout. The FM transition at n=6 doesn't register in the half-cut spectrum.

**Entanglement gap anti-correlates with energy gap** (Pearson r=-0.48): states with small energy gaps (critical) tend to have small entanglement gaps (many relevant eigenvalues). But the correlation is moderate — they encode different physics.

## Key Insights

### 1. Three Levels of Entanglement Description

| Level | What it captures | Archetype discrimination |
|-------|-----------------|-------------------------|
| **Scalar** (entropy) | Total entanglement amount | GHZ=Cluster=1.0, can't distinguish |
| **Correlation** (MI, I3, discord) | Entanglement topology/structure | All 5 archetypes distinguished |
| **Spectral** (eigenvalue distribution) | Symmetry structure, CFT content | New: doublet structure, power-law decay |

### 2. Spectral Fingerprints Are Symmetry Fingerprints

The XXZ doublet structure throughout the XY phase reveals that the entanglement spectrum directly reflects the Hamiltonian's symmetries (U(1) → paired eigenvalues). The BKT transition breaks this by opening a gap that destroys the Sz-sector mixing. **The entanglement spectrum is a symmetry probe, while I3/MI are topology probes.**

### 3. Archetype Boundaries ≠ Phase Boundaries ≠ Spectral Transitions

Three distinct "phase diagrams" coexist:
- **Thermodynamic**: FM/XY/Néel boundaries at Δ=-1, Δ=1
- **Archetype**: I3 sign change at Δ≈0.7 (inside XY phase)
- **Spectral**: degeneracy breaking at Δ≈0.95-1.0 (near BKT only)

This confirms and extends the 2025 finding (arXiv 2509.03588) that entanglement transitions can dissociate from bulk transitions. We add a third dissociation: archetype (correlation topology) transitions.

### 4. CFT Content at Criticality

The TFIM critical point's doublet structure (alternating spacings 1.0 and 2.75) encodes the c=1/2 Ising conformal field theory. This is information invisible to entropy, MI, I3, or any scalar measure — it requires the full spectral decomposition. Different universality classes should produce different doublet/multiplet patterns.

## Next Steps

- **Size scaling of doublet structure**: Do TFIM doublet spacings converge to CFT predictions at larger n?
- **Entanglement spectrum of XXZ at isotropic point** (Δ=0, SU(2) symmetry): expect higher degeneracies
- **Alternating (odd/even) bipartition**: should distinguish GHZ from Cluster 1D at spectral level
- **Entanglement Hamiltonian**: H_ent = -log(ρ_A) — does it look like a local Hamiltonian? (Bisognano-Wichmann theorem)
- **Potts model spectrum**: different CFT (c>1/2), should show different multiplet structure
