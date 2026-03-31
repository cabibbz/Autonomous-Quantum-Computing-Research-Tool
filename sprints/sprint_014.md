# Sprint 014 — Quantum Error Correction: Entanglement as Information Protection

**Date:** 2026-03-31
**Status:** In progress

## Motivation

We've spent 13 sprints building a toolkit of entanglement measures (MI, I3, discord, Phi, OTOC, entanglement spectrum) and studying how information spreads through quantum systems (scrambling, Hayden-Preskill). The natural culmination: **quantum error correction** — the technology that uses entanglement to protect information from noise.

Key connections to prior work:
- **Sprint 004/005:** I3 measures irreducible multipartite correlations — QEC codes need exactly this
- **Sprint 007/008:** 2D cluster states have topological protection — surface codes exploit this
- **Sprint 009:** Structured noise affects states differently — QEC must handle real noise
- **Sprint 010:** Stabilizer measurements detect both GME and errors — same math
- **Sprint 012/013:** Scrambling democratizes information — good codes should scramble

**Central question:** How does a QEC code's entanglement structure relate to its error-correcting ability? Can our information-theoretic toolkit predict code performance?

## Experiments

### 14a: QEC Code States — Entanglement Structure
- Implement 3-qubit bit-flip code and [[5,1,3]] perfect code
- Characterize code states with MI, I3, entanglement spectrum
- Compare to our archetypes (GHZ, W, Cluster)

### 14b: Error Correction in Action
- Apply single-qubit errors at varying rates, measure logical fidelity
- Compare corrected vs uncorrected logical error rates
- Find the "break-even" point where QEC becomes beneficial

### 14c: QEC Meets Scrambling
- Measure how code states spread information (OTOC-like)
- Compare information spreading in code states vs our archetypes
- Test whether better scramblers make better codes

---

## Results

### 14a: QEC Code States — Entanglement Structure

**3-qubit bit-flip code:**
- |0⟩_L = |000⟩ and |1⟩_L = |111⟩ are product states — zero entanglement
- |+⟩_L = (|000⟩ + |111⟩)/√2 is literally a **GHZ state**: MI=1.0 all pairs, I3=0, negativity=0.5
- The simplest QEC code IS the simplest entangled state — they're the same thing

**[[5,1,3]] perfect code — remarkable results:**
- ALL single-qubit entropies = 1.0 (maximally mixed) — every qubit looks random alone
- **ZERO pairwise MI for ALL 10 pairs** — identical to 2D cluster behavior (Sprint 010)!
- **I3 = -1.0 for ALL 10 triples** — uniformly maximal irreducible 3-body correlations
- Negativity spectrum perfectly uniform: |A|=1 → 0.5, |A|=2 → 1.5 for every bipartition
- These properties are **logical-state independent** — |0⟩_L, |1⟩_L, and |+⟩_L all identical

**Comparison table (n=5):**

| Property | GHZ | W | Cluster_1D | [[5,1,3]] |
|---|---|---|---|---|
| S(single) | 1.0 | 0.72 | 1.0 | 1.0 |
| MI (pairs) | 1.0 uniform | 0.47 uniform | sparse (0 or 1) | **0.0 uniform** |
| I3 (triples) | +1.0 uniform | +0.22 uniform | mixed (-1 or 0) | **-1.0 uniform** |
| Neg |A|=1 | 0.5 uniform | 0.4 uniform | 0.5 uniform | 0.5 uniform |
| Neg |A|=2 | 0.5 uniform | 0.49 uniform | mixed (0.5-1.5) | **1.5 uniform** |

**Key insight:** The [[5,1,3]] code is a "perfected cluster state" — it has the irreducible 3-body structure of the cluster but with **full permutation symmetry**. Cluster has I3=-1.0 only for consecutive triples; the QEC code achieves I3=-1.0 for ALL triples. This is the information-theoretic signature of error correction: every pair reveals nothing (zero MI), but every triple has maximal irreducible correlation. Information is perfectly delocalized across 3+ body correlations — no single-qubit error can access or corrupt the encoded information.

The body-order hierarchy from Sprint 010 now has a QEC interpretation:
- 3-qubit code = GHZ-type (2-body correlations → corrects only 1 error type)
- [[5,1,3]] code = symmetric 3-body (→ corrects ALL single-qubit errors)
