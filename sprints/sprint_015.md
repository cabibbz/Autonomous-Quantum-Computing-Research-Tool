# Sprint 015 — Code Diversity: Does Distance Alone Determine Information Structure?

**Date:** 2026-03-31
**Status:** In progress

## Motivation

Sprint 014 revealed that the [[5,1,3]] perfect code is a "perfected cluster state" — zero pairwise MI, I3=-1.0 for ALL triples, perfectly uniform negativity spectrum, and a sharp Page curve at k=3. But is this universal for all distance-3 codes, or specific to [[5,1,3]]?

Three distance-3 codes with very different architectures:
- **[[5,1,3]] perfect code** — non-CSS, maximally efficient (Hamming bound saturated)
- **Steane [[7,1,3]] code** — CSS construction (X and Z errors corrected independently), based on classical Hamming code
- **Shor [[9,1,3]] code** — concatenated (3-qubit phase-flip wrapping 3-qubit bit-flip), the first QEC code ever

If all three share the same MI/I3/Page curve → distance alone determines information structure.
If they differ → code architecture shapes how information is distributed, beyond just the distance.

Key connections to prior work:
- **Sprint 014:** [[5,1,3]] has I3=-1.0 for all triples, MI=0 for all pairs
- **Sprint 013:** Page curve from Hayden-Preskill mirrors QEC recovery threshold
- **Sprint 010:** Body-order hierarchy: code distance = minimum Pauli weight for logical detection

**Central question:** Do all distance-3 codes distribute information identically, or does the code's internal architecture create distinct information topologies?

## Experiments

### 15a: Entanglement Structure of Three Code Families
- Build Steane [[7,1,3]] and Shor [[9,1,3]] codes via stabilizer projection
- Measure MI, I3, negativity spectrum for each
- Compare to [[5,1,3]] results from Sprint 014

### 15b: Error Correction Performance Comparison
- Apply depolarizing noise at varying rates, decode, measure logical fidelity
- Compare break-even thresholds across all three codes
- Connect performance differences to entanglement structure

### 15c: Page Curves — Information Recovery Across Code Families
- Entangle reference qubit with logical qubit, measure MI(ref : k physical qubits)
- Map the recovery transition for each code
- Test: does Singleton bound (k = n-d+1) predict the transition universally?

---

## Results

### 15a: Entanglement Structure of Three Code Families

**Universal properties (all three distance-3 codes share these):**
- All single-qubit entropies = 1.0 (maximally mixed marginals)
- Negativity |A|=1 = 0.5 uniformly (every qubit equally entangled with rest)
- Logical-state independent for MI and negativity

**Where they diverge dramatically:**

| Property | [[5,1,3]] | Steane [[7,1,3]] | Shor [[9,1,3]] |
|---|---|---|---|
| MI range | 0.0 (all pairs) | 0.0 (all pairs) | 0.0–1.0 |
| Total MI | 0.0 | 0.0 | **9.0** |
| I3 range | -1.0 (all triples) | -1.0 to 0.0 | -1.0 to **+1.0** |
| avg I3 | **-1.0** | -0.2 | -0.286 |
| Neg |A|=2 | 1.5 (uniform) | 1.5 (uniform) | 0.5–1.5 |
| Permutation symmetric? | YES | NO | NO |

**Shor code reveals hierarchical structure:**
The MI pattern perfectly exposes the 3-block concatenation:
- Within-block pairs (e.g., 0-1, 3-4): MI = 1.0 (GHZ-like)
- Between-block pairs (e.g., 0-3, 1-6): MI = 0.0

I3 has three categories:
- All 3 qubits from same block (e.g., 0-1-2): I3 = **+1.0** (redundant/GHZ-like) — 3 triples
- One qubit from each block (e.g., 0-3-6): I3 = **-1.0** (irreducible) — 27 triples
- Two from one block, one from another: I3 = 0.0 — 54 triples

The Shor code is a **hierarchical scrambler**: GHZ within blocks (classical broadcast), scrambled across blocks (quantum irreducible). The phase-flip layer creates irreducible 3-body correlations between blocks while preserving within-block redundancy.

**Steane code has selective I3:**
Only 7 of 35 triples have I3 = -1.0; the remaining 28 have I3 = 0.0. Zero MI everywhere (like [[5,1,3]]) but the irreducible correlations are sparse, not uniform. The Hamming code geometry determines which triples carry irreducible information.

**Key insight:** Distance alone does NOT determine information structure. Three codes with identical distance (d=3) have fundamentally different information topologies:
- **[[5,1,3]]**: Perfect information democracy — ALL correlations are irreducible 3-body, fully symmetric
- **Steane [[7,1,3]]**: Selective democracy — zero pairwise leakage but non-uniform 3-body correlations
- **Shor [[9,1,3]]**: Hierarchical — massive within-block pairwise leakage (MI=1.0), irreducible correlations only across blocks

The Shor code's total MI of 9.0 means **each block internally reveals the logical state** — it's three correlated GHZ states. This suggests it should be weaker under noise that can exploit within-block correlations.

