# Sprint 014 — Quantum Error Correction: Entanglement as Information Protection

**Date:** 2026-03-31
**Status:** Complete (3/3 experiments)

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

### 14b: Error Correction Performance

**3-qubit bit-flip code under bit-flip noise:**
- Works perfectly as designed: at p=0.1, fidelity 0.972 (coded) vs 0.900 (uncoded)
- Coded ALWAYS better than uncoded up to p=0.5
- |+⟩ state maintains fidelity 1.0 at all noise levels (|+⟩ is an X eigenstate — bit-flips don't affect it)
- Theoretical: corrects any single bit-flip, fails only on 2+ flips → P_fail = 3p²(1-p) + p³

**[[5,1,3]] code under depolarizing noise:**
- **Break-even at p ≈ 14%** — below this threshold, the code helps
- At p=0.1: fidelity 0.947 (coded) vs 0.933 (uncoded)
- At p=0.3: coded is WORSE (0.712 vs 0.800) — too many multi-qubit errors overwhelm correction
- Performance is **logical-state independent** — |0⟩, |1⟩, |+⟩ all identical (confirms symmetric protection)

**3-qubit code under depolarizing noise (control experiment):**
- **WORSE than uncoded** even at p=0.1 (0.826 vs 0.933)
- The bit-flip correction introduces logical errors when phase errors occur
- Confirms: entanglement structure determines error correction capability

**Key findings:**
1. The break-even threshold (14%) is close to the CHSH violation death point (9.5% from Sprint 002) — both measure when quantum advantage disappears
2. Code performance is state-independent for the [[5,1,3]] code — a direct consequence of its uniform entanglement structure (all-symmetric MI=0, I3=-1.0)
3. The 3-qubit code's failure under depolarizing noise is predicted by its GHZ-like entanglement: 2-body (ZZ) correlations can only detect Z-type errors, leaving X/Y errors invisible
4. Error correction works when noise rate is below the code's "entanglement capacity" — the rate at which multi-qubit errors overwhelm the syndrome space

### 14c: QEC Meets Scrambling — Information Spreading in Code States

**MI recovery pattern (reference qubit entangled with logical, measured against k physical qubits):**

| k qubits | [[5,1,3]] | 3-qubit (GHZ) | GHZ (n=5) |
|---|---|---|---|
| 1 | **0.0** | 1.0 | 1.0 |
| 2 | **0.0** | 1.0 | 1.0 |
| 3 | **2.0** | 2.0 | 1.0 |
| 4 | 2.0 | — | 1.0 |
| 5 | 2.0 | — | 2.0 |

**This is the Page curve from Sprint 013!** The [[5,1,3]] code has a sharp phase transition at k=3 (just over half): below, zero information; above, full recovery. GHZ leaks MI=1.0 to every single qubit but only reaches full recovery (MI=2.0) when ALL qubits are available.

**Operator spreading:**
- [[5,1,3]]: ALL single-qubit ⟨Z⟩ = 0.0. ALL two-qubit ⟨Z_iZ_j⟩ = 0.0. Logical Z (= ZZZZZ) has been spread to weight 5 — invisible at 1-body and 2-body level.
- 3-qubit code: ⟨Z⟩ = ±1.0 on every physical qubit — logical state immediately visible. Zero spreading.

**Quantum Singleton bound verified:**
- [[5,1,3]] code: need n-d+1 = 3 qubits for recovery → confirmed: MI(ref : any 3) = 2.0 exactly
- Any 2 qubits: MI = 0.0 exactly — completely decoupled from logical information
- This is the decoupling theorem in action

**Connection to previous sprints:**

1. **Sprint 013 (Hayden-Preskill):** Scrambling converts position-dependent recovery to position-independent. The [[5,1,3]] code does the same WITHOUT dynamics — it's a "frozen scrambler." Any 3 qubits are equally useful, just as any late-radiation qubit was equally useful post-scrambling.

2. **Sprint 010 (Pauli profiles):** The 2D cluster was "invisible" at the 2-body level. The [[5,1,3]] code is invisible at the 1-body AND 2-body level. The minimum body-order to detect the logical information equals the code distance. **Code distance = body-order of information.**

3. **Sprint 012 (Scrambling):** GHZ had OTOC=1.0 (zero scrambling) despite maximal MI. Here, GHZ/3-qubit code leaks MI=1.0 to every single qubit — this is WHY it can't correct general errors. Information that isn't spread can be locally corrupted.

**The grand unification:**

| Concept | Static (QEC) | Dynamic (Scrambling) |
|---|---|---|
| Information spreading | Encoding circuit | Random circuit evolution |
| Recovery threshold | Quantum Singleton bound (n-d+1) | Decoupling theorem (n/2 + O(1)) |
| Position independence | Any k ≥ n-d+1 qubits work | Any late-radiation qubit works |
| Page curve | MI=0 below threshold, MI=2 above | MI=0 for early radiation, MI=2 combined |
| Information security | d-1 qubits reveal nothing | Pre-scrambling qubits reveal nothing |

**Key insight:** QEC codes are static scramblers. The encoding circuit is a designed, efficient scrambler that achieves the same information-theoretic properties (decoupling, democratization, Page-curve behavior) that random circuits achieve through depth. The code distance is the "scrambling radius" — the minimum number of qubits that can access the logical information. This connects quantum computing (QEC), quantum gravity (Hayden-Preskill, Page curve), and information theory (channel capacity) through a single framework: **how information is distributed across subsystems.**
