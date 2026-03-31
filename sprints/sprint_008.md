# Sprint 008 — Integrated Information Theory: How "Whole" Are Quantum States?

**Date:** 2026-03-31
**Status:** In progress

## Motivation

Integrated Information Theory (IIT) proposes that consciousness arises from systems that are "more than the sum of their parts" — quantified by Phi (Φ), the amount of information lost when the system is split at its weakest point. Setting aside the consciousness angle, Phi is a powerful measure of system integration that should complement our existing toolkit (MI, I3, discord, negativity).

For pure quantum states, the math simplifies beautifully:
- I(A:Ā) = S(A) + S(Ā) - S(AĀ) = 2·S(A) since S(A)=S(Ā) for pure states and S(total)=0
- Φ = min over all bipartitions of I(A:Ā) = 2 · min S(A)

So Phi is determined entirely by the **weakest bipartition** — the cut that loses the least information.

**Predictions:**
- GHZ: Flat entanglement spectrum (S=1.0 everywhere), so Φ = 2.0
- Cluster 1D: Variable spectrum with weak links for distant qubits, so Φ < 2.0
- Cluster 2D: Eliminated weak bipartitions (Sprint 007b), so Φ should be higher than 1D
- W: Entropy grows with cut size, so |A|=1 cut is weakest

**Key question:** Does Phi correlate or anti-correlate with I3? GHZ has positive I3 (redundant) but may have HIGH Phi (uniformly integrated). Cluster has negative I3 (irreducible) but may have LOW Phi (sparse connectivity = weak links). If so, Phi and I3 measure fundamentally different aspects of "wholeness."

---

## Experiment 8a: Phi and MI Spectrum for All Archetypes

**Goal:** Compute Phi (quantum integrated information) for GHZ, W, Cluster-1D, Cluster-2D at n=6. Also compute the full MI spectrum — distribution of I(A:Ā) across all 31 non-trivial bipartitions.

**Results:**

**Phi and MI spectrum (n=6):**

| State | Φ (raw) | Φ (normalized) | MI range | MI mean | MI std |
|-------|---------|----------------|----------|---------|--------|
| GHZ | 2.0 | 0.667 | [2.0, 2.0] | 2.0 | 0.0 |
| W | **1.3** | 0.667 | [1.3, 2.0] | 1.79 | 0.25 |
| Cluster 1D | 2.0 | 0.667 | [2.0, 6.0] | 3.68 | 1.25 |
| Cluster 2D | 2.0 | **1.333** | [2.0, 6.0] | 4.0 | 1.24 |

**Key findings:**

1. **Raw Phi fails to distinguish GHZ from Cluster.** Both have S=1.0 at every single-qubit cut, so raw Phi = 2.0 for both. The "weakest link" is always a single-qubit cut.

2. **W has the lowest Phi (1.3).** Its single-qubit entropy is only 0.65 bits, creating a genuinely weak link. Despite having the only nonzero discord (Sprint 005), W is the *least integrated* system.

3. **Normalized Phi distinguishes 2D Cluster.** Φ_norm accounts for partition size: MI/min(|A|,|Ā|). For all states except 2D cluster, the equal-size partition (|A|=3) has Φ_norm = 0.667. But 2D cluster's minimum |A|=3 MI is 4.0 (not 2.0), giving Φ_norm = 1.333 — **2x higher**. The 2D geometry eliminates weak links at every scale.

4. **MI variance reveals structure.** GHZ: zero variance (perfectly homogeneous — every cut sees the same). Cluster: massive variance (geometry creates rich structure). MI range of [2.0, 6.0] for cluster means some cuts see 3x more information than others.

5. **Total MI (mean) ranks differently from Phi.** Cluster_2D > Cluster_1D > GHZ > W in total integration. But in Phi (weakest link): GHZ = Cluster = 2.0 > W = 1.3. The average vs minimum distinction matters.

**Insight:** Phi measures *uniformity of integration* at its weakest point, not total information content. GHZ is "shallowly but uniformly integrated" — every cut sees exactly 1 bit. Cluster is "deeply but heterogeneously integrated" — some cuts see 3 bits but all see at least 1. The normalized Phi captures 2D's advantage: it raises the floor at EVERY scale, not just single-qubit cuts.

---

## Experiment 8b: Phi Under Depolarizing Noise

**Goal:** How does Phi degrade under noise? Compare to half-cut negativity degradation (Sprint 006b). Does integration die faster or slower than entanglement?

**Results:** *(pending)*

---

## Experiment 8c: Phi Under Qubit Loss

**Goal:** Phi for mixed states after tracing out qubits. For mixed states, I(A:Ā) = S(A) + S(Ā) - S(AĀ), and S(AĀ) ≠ 0. Does 2D cluster maintain higher Phi than 1D under loss?

**Results:** *(pending)*
