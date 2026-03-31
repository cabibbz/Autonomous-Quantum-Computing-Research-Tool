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

### 15b: Error Correction Performance Comparison

Depolarizing noise at p=0 to 0.3, comparing coded (encode → noise → syndrome correct → decode) vs uncoded fidelity.

**Break-even thresholds (coded fidelity = uncoded fidelity):**

| Code | Break-even p | Fidelity @ p=0.1 | Fidelity @ p=0.2 | Fidelity @ p=0.3 |
|---|---|---|---|---|
| Uncoded | — | 0.933 | 0.867 | 0.800 |
| [[5,1,3]] | **~0.15** | 0.947 | 0.834 | 0.712 |
| Steane [[7,1,3]] | **~0.22** | 0.966 | 0.887 | 0.778 |
| Shor [[9,1,3]] | **>0.30** | 0.978 | 0.921 | 0.831 |

**Key findings:**
- Performance is completely logical-state independent (|0⟩ and |+⟩ give identical curves for all codes)
- **Shor code has the BEST performance** despite having the "worst" information topology (MI=9.0 of pairwise leakage)
- Ranking by break-even: Shor (>0.30) >> Steane (~0.22) >> [[5,1,3]] (~0.15)
- At low noise (p=0.02): Shor fidelity 0.999, Steane 0.998, [[5,1,3]] 0.997 — all excellent

**The surprise:** 15a predicted Shor should be weakest (within-block MI leakage exposes logical state). But Shor is actually **strongest** under depolarizing noise. Why?

**Resolution:** More physical qubits = more redundancy = more error-correcting power, even with information leakage. The Shor code's 9 qubits provide more "room" for noise to act without corrupting the logical state. The relevant comparison isn't break-even threshold but **break-even threshold per physical qubit**:
- [[5,1,3]]: 0.15/5 = 0.030 per qubit — most **efficient**
- Steane: 0.22/7 = 0.031 per qubit
- Shor: 0.30/9 = 0.033 per qubit

Per-qubit efficiency is remarkably similar (~0.03), suggesting that for depolarizing noise, raw qubit count matters more than information topology. The topology differences from 15a would matter for **structured** noise (cf. Sprint 009).

### 15c: Page Curves — Information Recovery Across Code Families

Encoded Bell state: (|0⟩_ref|0_L⟩ + |1⟩_ref|1_L⟩)/√2. Measure MI(reference : k physical qubits) for all subsets of size k.

**Results:**

| k | [[5,1,3]] | Steane [[7,1,3]] | Shor [[9,1,3]] |
|---|---|---|---|
| 1 | 0.000 | 0.000 | 0.000 |
| 2 | 0.000 | 0.000 | 0.000 |
| 3 | **2.000** (all subsets) | 0.400 (range 0–2) | 0.360 (range 0–1) |
| 4 | 2.000 | 1.600 (range 0–2) | 0.800 (range 0–1) |
| 5 | 2.000 | **2.000** (all subsets) | 1.240 (range 1–2) |
| 6 | — | 2.000 | 1.700 (range 1–2) |
| 7 | — | 2.000 | **2.000** (all subsets) |

**Singleton bound k ≥ n-d+1 = n-2 is EXACT for all three codes:**
- [[5,1,3]]: k ≥ 3 → MI jumps from 0 to 2.0. **All** 3-qubit subsets recover fully. Sharp transition.
- Steane [[7,1,3]]: k ≥ 5 → MI reaches 2.0 for all subsets. But at k=3-4, SOME subsets recover (MI=2.0) while others get NOTHING (MI=0.0).
- Shor [[9,1,3]]: k ≥ 7 → MI reaches 2.0 for all subsets. At k=3-6, maximum MI is only 1.0 (not 2.0), showing partial recovery capped by block structure.

**Three qualitatively different Page curves:**

1. **[[5,1,3]] — Perfect knife-edge:** MI = 0 for k < 3, MI = 2.0 for k ≥ 3. Zero spread at all k. Every subset of the same size gives identical recovery. This is the information-theoretic definition of a "perfect" code — the transition is maximally sharp and subset-independent.

2. **Steane [[7,1,3]] — Bimodal transition:** At k=3, some subsets recover everything (MI=2.0) while others recover nothing (MI=0.0). The transition is spread over k=3-4 but with binary outcomes — each subset either works or doesn't. This reflects the Hamming code geometry: only certain 3-qubit subsets contain a full error-correcting set.

3. **Shor [[9,1,3]] — Gradual staircase:** MI increases through intermediate values (0.36, 0.80, 1.24, 1.70) with max MI capped at 1.0 until k=7. This reflects the hierarchical block structure — you can recover block-level information (MI=1.0) before recovering full logical information (MI=2.0). The concatenation creates two recovery thresholds: within-block (partial) and across-block (full).

**Key insight:** The Singleton bound is universal (all three codes obey k ≥ n-2 exactly), but the **shape** of the Page curve is code-specific and directly reflects the information topology discovered in 15a:
- [[5,1,3]]'s perfect symmetry → perfect knife-edge
- Steane's selective I3 → bimodal (some subsets work, others don't)
- Shor's hierarchical MI → gradual staircase with intermediate plateaus

The Page curve shape is a new diagnostic that captures code architecture beyond distance alone. It's the static analog of the Hayden-Preskill protocol from Sprint 013 — but here the "scrambling" is the encoding circuit, and the "recovery transition" is the Singleton bound.

---

## Summary & Key Insights

**Central finding:** Distance alone does NOT determine information structure. Three [[n,1,3]] codes with identical distance have fundamentally different information topologies, noise performance, and recovery characteristics.

**Three information archetypes for QEC codes:**
1. **Democratic** ([[5,1,3]]): All correlations irreducible 3-body, fully symmetric, knife-edge Page curve
2. **Selective** (Steane [[7,1,3]]): Zero pairwise leakage but non-uniform 3-body correlations, bimodal recovery
3. **Hierarchical** (Shor [[9,1,3]]): Strong within-block pairwise correlations, gradual staircase recovery with intermediate plateaus

**The efficiency-topology tradeoff:**
- [[5,1,3]] is the most qubit-efficient (5 qubits) with the sharpest Page curve, but breaks even at the lowest noise (p≈0.15)
- Shor is the least qubit-efficient (9 qubits) with the broadest Page curve, but tolerates the most noise (p>0.30)
- Per-qubit break-even efficiency is nearly identical (~0.03) — symmetric noise can't distinguish topologies

**Connection to prior work:**
- The Page curve shapes parallel Hayden-Preskill recovery (Sprint 013) — encoding IS scrambling
- The topology archetypes (Democratic/Selective/Hierarchical) mirror the entanglement archetypes from Sprint 005-006 (GHZ/W/Cluster) but for codes instead of states
- Shor's hierarchical structure (within-block GHZ) confirms the Sprint 014 finding that GHZ = classical broadcast — the within-block MI=1.0 is GHZ-like redundancy, not quantum irreducibility

