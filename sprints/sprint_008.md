# Sprint 008 — Integrated Information Theory: How "Whole" Are Quantum States?

**Date:** 2026-03-31
**Status:** Complete (3/3 experiments)

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

**Results:**

**Phi and negativity under depolarizing noise (n=6):**

| State | Φ(p=0) | Neg(p=0) | Φ(p=0.1) | Neg(p=0.1) | Death(Φ) | Death(Neg) |
|-------|--------|----------|----------|------------|----------|------------|
| GHZ | 2.0 | 0.5 | 1.76 | 0.45 | p=1.0 | p=1.0 |
| W | 1.3 | 0.5 | 1.17 | 0.45 | p=1.0 | p=1.0 |
| Cluster 1D | 2.0 | 0.5 | 1.76 | 0.45 | p=1.0 | p=1.0 |
| Cluster 2D | 2.0 | **3.5** | 1.76 | **3.1** | p=1.0 | p=0.9 |

**Key findings:**

1. **GHZ and Cluster_1D have identical Phi curves.** At every noise level, their Phi values match exactly. This is because both have the same MIP (any single-qubit cut) and S(single qubit) = 1.0 for both. Depolarizing noise affects them identically at this scale.

2. **W always has the lowest Phi** — consistently ~65% of GHZ/Cluster at every noise level.

3. **No interesting phase transitions.** Under uniform depolarizing noise, Phi decays smoothly for all states and only reaches zero at p=1.0 (fully mixed). This confirms that depolarizing noise is "too uniform" to reveal structural differences in integration.

4. **2D Cluster's advantage shows in negativity, not Phi.** Half-cut negativity starts at 3.5 (7x GHZ) and remains much higher throughout. But Phi is identical to GHZ/1D-Cluster because the MIP is still a single-qubit cut.

5. **Depolarizing noise preserves ranking.** Phi ordering (GHZ = Cluster ≥ W) is invariant under noise. The noise doesn't create any crossovers.

**Insight:** Depolarizing noise is too "democratic" — it affects all qubits equally, so it can't distinguish states whose weakest link is the same (single-qubit cut). Local or structured noise (amplitude damping, dephasing on specific qubits) would be more discriminating. This experiment confirms Phi's limitation: it only sees the weakest link, not the overall structure.

---

## Experiment 8c: Phi Under Qubit Loss

**Goal:** Phi for mixed states after tracing out qubits. For mixed states, I(A:Ā) = S(A) + S(Ā) - S(AĀ), and S(AĀ) ≠ 0. Does 2D cluster maintain higher Phi than 1D under loss?

**Results:**

**Single qubit loss (n=6 → n=5):**

| State | Φ(baseline) | Φ(1 lost) | Retention | Position-dependent? |
|-------|-------------|-----------|-----------|---------------------|
| GHZ | 2.0 | 1.0 (uniform) | 50% | No |
| W | 1.3 | 0.918 (uniform) | 70.6% | No |
| Cluster 1D | 2.0 | 1.0–2.0 | 58.3% avg | **YES** |
| Cluster 2D | 2.0 | **2.0 (uniform)** | **100%** | No |

**Two qubit loss (n=6 → n=4):**

| State | Φ(2 lost) | Retention | Position-dependent? |
|-------|-----------|-----------|---------------------|
| GHZ | 1.0 (uniform) | 50% | No |
| W | 0.732 (uniform) | 56.3% | No |
| Cluster 1D | **0.0–1.0** | 33.3% avg | **YES** (Φ=0 for 5 of 15 pairs!) |
| Cluster 2D | 1.0–2.0 | 60% avg | Mild (3 pairs retain Φ=2.0) |

**Key findings:**

1. **2D Cluster is perfectly robust to single qubit loss.** Phi stays at 2.0 regardless of which qubit is lost — the ONLY state with 100% Phi retention. This mirrors the position-independent entropy behavior found in Sprint 007c.

2. **1D Cluster can be completely shattered.** Losing 2 qubits that break the chain (e.g., q0+q2, q2+q3) can drive Phi to **zero** — the system becomes fully decomposable, with zero integration. 5 of 15 possible two-qubit losses completely destroy integration.

3. **GHZ is robustly mediocre.** Loses exactly 50% of Phi on any single qubit loss, uniformly. Never reaches zero. Its rigid structure means uniform degradation — no catastrophic failure, but no resilience either.

4. **W retains the highest FRACTION under single loss (70.6%)** despite starting with the lowest Phi. Its democratic pairwise structure means losing one qubit removes only 1/n of the correlations. But it starts so low that its absolute Phi is still the smallest.

5. **Cluster 1D is the most fragile.** Mean Phi retention of 33.3% under two-qubit loss, with complete disintegration possible. The chain topology has critical links that, when both are broken, fragment the system into uncorrelated subsystems.

**Insight:** Phi under qubit loss reveals a fundamental distinction between **topological** and **geometric** robustness:
- 2D Cluster: Topological — multiple paths connect any pair of qubits, so no single loss can isolate a subsystem. Phi is perfectly maintained.
- 1D Cluster: Geometric — linear chain has critical links. Breaking two links can fragment the chain, driving Phi to zero.
- GHZ: Algebraic — all qubits contribute equally to one global superposition. Losing any qubit degrades the superposition uniformly.
- W: Democratic — correlations are distributed pairwise, so loss is proportional.

This directly connects to quantum error correction: surface codes (built on 2D cluster topology) can tolerate local errors precisely because no local damage can fragment the code's integrated information structure.

---

## Sprint 008 Summary

**Three experiments on Quantum Integrated Information (Phi):**

1. **8a (Phi spectrum):** Raw Phi can't distinguish GHZ from Cluster (both = 2.0). Normalized Phi reveals 2D Cluster is 2x more integrated. MI variance captures structural differences.

2. **8b (Phi under noise):** Depolarizing noise is too uniform to discriminate — GHZ and Cluster_1D have identical Phi curves. No interesting phase transitions.

3. **8c (Phi under qubit loss):** The discriminating experiment. 2D Cluster is **perfectly robust** (100% Phi retention under single loss). 1D Cluster is **catastrophically fragile** (Phi → 0 for certain two-qubit losses). This connects Phi directly to topological protection in quantum error correction.

**The big picture:** Phi measures system *integration* — how much the whole exceeds the sum of its parts at the weakest point. The hierarchy of robustness is:

**2D Cluster > GHZ > W > 1D Cluster** (under qubit loss)

This is different from every other ranking we've found:
- Entanglement (negativity): 2D Cluster >> GHZ = W = 1D Cluster
- Discord: W > GHZ = Cluster = 0
- Tripartite info: Cluster (I3=-1) > W (I3≈0) > GHZ (I3=+1)

Phi under loss is the first measure that makes 1D Cluster look *worse* than GHZ. The chain topology's advantage (rich local entanglement, negative I3) comes at the cost of fragility to targeted damage.

**Next sprint ideas:**
- Local noise models (amplitude damping on specific qubits) — will structured noise discriminate better than depolarizing?
- Genuine multipartite entanglement (GME) witnesses — can we detect entanglement that survives all bipartitions?
- Real IBM hardware comparison — CHSH or small-state tomography on actual QPU
- Quantum error correction codes: compare [[5,1,3]] code integration to these natural states
