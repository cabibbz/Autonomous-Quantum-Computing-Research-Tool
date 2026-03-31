# Sprint 006 — Entanglement Monotones: Concurrence, Negativity & Monogamy

**Date:** 2026-03-31
**Goal:** Add proper entanglement monotones (concurrence, negativity) to the archetype classification. These quantify *how much* entanglement exists in each pair, complementing MI (total correlations), I3 (structure), and discord (classical vs quantum nature).

**Background:**
- **Concurrence** (Wootters 1998): defined for 2-qubit states. C=0 means separable, C=1 means maximally entangled. Gold standard for bipartite entanglement.
- **Negativity** (Vidal & Werner 2002): sum of negative eigenvalues of partial transpose. Works for any bipartition. N=0 means PPT (no distillable entanglement), N=0.5 means maximally entangled pair.
- **Residual tangle** (Coffman-Kundu-Wootters): τ = C²(A|BC) - C²(A|B) - C²(A|C). Quantifies genuine tripartite entanglement not accounted for by pairwise entanglement.

**Predictions:**
- GHZ pairwise RDMs are classical mixtures (discord=0), so concurrence should be 0 — all entanglement is multipartite ✓
- W state has nonzero discord, so should have nonzero concurrence — entanglement is distributed in pairs ✓
- Cluster state has zero pairwise discord but strong I3 — unclear! ✓ (turned out: zero pairwise, but massive multipartite)

## Experiment 6a: Pairwise Concurrence & Negativity

**Status:** Complete

| State | n=4 C | n=6 C | n=8 C | Pattern |
|-------|-------|-------|-------|---------|
| GHZ   | 0.000 | 0.000 | 0.000 | Zero everywhere — no pairwise entanglement |
| W     | 0.500 | 0.333 | 0.250 | Uniform C=2/n across ALL pairs |
| Cluster | 0.000 | 0.000 | 0.000 | Zero everywhere — like GHZ |

**Key findings:**
- W state concurrence follows exact formula **C = 2/n** for all pairs, decreasing as entanglement is shared among more qubits
- GHZ and Cluster have ZERO pairwise entanglement — all their entanglement is genuinely multipartite
- Concurrence alone cannot distinguish GHZ from Cluster (both zero)

## Experiment 6b: Negativity Under Noise

**Status:** Complete (but non-discriminating)

Half-cut negativity (3-vs-3 at n=6) under depolarizing noise:
- All three states have identical decay: N(p) = 0.5 × (1-p)^(64/63)
- All die at p ≈ 0.98 — depolarizing noise is too symmetric to discriminate

W pairwise concurrence under noise:
- Dies at p ≈ 0.14 — much more fragile than global entanglement
- This makes sense: pairwise correlations are weak (C=0.33 at n=6)

**Lesson:** Depolarizing noise + half-cut = maximally uninformative combination. Need asymmetric bipartitions or asymmetric noise to see differences.

## Experiment 6c: Entanglement Spectrum (All Bipartitions)

**Status:** Complete

Negativity across all bipartitions of n=6, grouped by subsystem size |A|:

| State | |A|=1 | |A|=2 | |A|=3 |
|-------|-------|-------|-------|
| GHZ   | 0.50 ± 0.00 | 0.50 ± 0.00 | 0.50 ± 0.00 |
| W     | 0.37 ± 0.00 | 0.47 ± 0.00 | 0.50 ± 0.00 |
| Cluster | 0.50 ± 0.00 | 1.37 ± 0.34 | 2.20 ± 1.10 |

**Key findings:**
- **GHZ is flat**: negativity = 0.5 for every bipartition. Perfectly symmetric — doesn't matter how you cut it.
- **W grows with cut size**: less entangled for small subsystems, converges to 0.5 at half-cut. Entanglement is "spread thin."
- **Cluster is explosive**: negativity ranges from 0.5 to **3.5** depending on which qubits you group together! Mean negativity at |A|=3 is 2.2 — over 4x GHZ! And the variance is huge (std=1.1), meaning geometry matters enormously.

The Cluster state's maximum negativity of 3.5 at |A|=3 means certain 3-vs-3 cuts create **7x more entanglement** than GHZ. This connects to Sprint 003's finding that qubit loss can *increase* entropy — the Cluster state has "hidden" entanglement that depends on how you carve up the system.

## Concurrence Monogamy (CKW Inequality)

W state (n=4): For every qubit, the residual tangle τ = 0.000 exactly.

**This is profound:** The W state saturates the monogamy inequality — it has **zero genuine tripartite entanglement**. ALL its entanglement is distributed in pairwise bonds. Compare:
- **GHZ:** Zero pairwise entanglement, all multipartite → τ > 0
- **W:** All pairwise entanglement, zero multipartite → τ = 0

These are opposite extremes of the same monogamy constraint.

## Updated Archetype Table

| Measure | GHZ | W | Cluster |
|---------|-----|---|---------|
| Pairwise MI | 1.0 (uniform) | 0.38 (uniform) | Sparse (NN only) |
| Tripartite info I3 | +1.0 (redundant) | ? | -1.0 (synergistic) |
| Pairwise discord | 0 (classical) | 0.28 (74% quantum) | 0 (classical) |
| Pairwise concurrence | 0 | 2/n (uniform) | 0 |
| Residual tangle | >0 (all multipartite) | 0 (all pairwise) | >0 (all multipartite) |
| Entanglement spectrum | Flat (0.5 everywhere) | Growing (0.37→0.50) | Explosive (0.5→3.5) |

## Synthesis

The entanglement spectrum is the most powerful single discriminator we've found. It reveals three fundamentally different *topologies* of entanglement:

1. **GHZ = Democratic entanglement.** Every cut yields the same amount. It's a single "entanglement bond" shared equally among all qubits. Fragile: lose one qubit and the whole bond breaks (Sprint 003).

2. **W = Distributed entanglement.** Entanglement lives in the pairwise bonds (C=2/n), not in multipartite correlations (τ=0). The more qubits you cut away, the more pairs you include, so negativity grows. Robust: lose a qubit and the remaining pairs survive.

3. **Cluster = Geometric entanglement.** Entanglement depends dramatically on *which* qubits you group together. Some cuts yield 7x more than GHZ. The chain geometry creates "entanglement hotspots" at specific bipartitions. This is why it's useful for measurement-based QC — the right measurements unlock the right entanglement.

## Next Sprint Ideas

- **I3 for W state** to complete the archetype table (currently "?")
- **2D cluster states**: Does the entanglement spectrum change qualitatively with lattice geometry?
- **Entanglement spectrum under qubit loss**: How does losing qubits reshape the spectrum?
- **Robustness under local (non-depolarizing) noise**: Amplitude damping, phase damping — asymmetric noise that might discriminate better
- **Genuine multipartite entanglement witnesses**: Beyond CKW, use GME witnesses for GHZ and Cluster
