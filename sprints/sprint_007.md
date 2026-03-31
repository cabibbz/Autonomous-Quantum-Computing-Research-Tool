# Sprint 007 — 2D Cluster States: Does Geometry Shape Entanglement?

**Date:** 2026-03-31
**Status:** In progress

## Motivation

We've established three entanglement archetypes (GHZ=democratic, W=distributed, Cluster=geometric) using 1D chain cluster states. But real quantum error-correcting codes (surface codes, toric codes) use **2D** cluster states. Does the jump from 1D→2D fundamentally change the entanglement topology?

Key questions:
1. Does 2D geometry create qualitatively different I3 patterns vs 1D chains?
2. Is the entanglement spectrum more explosive in 2D (more bipartition variety)?
3. Does qubit loss in 2D depend on topological position (corner vs edge vs interior)?

Also: filling the W-state I3 gap from Sprint 005 to complete the archetype table.

---

## Experiment 7a: W-State I3 + 2D Cluster MI & I3

**Goal:** Complete the I3 archetype table. Build 2x3 and 2x4 2D cluster states and compute pairwise MI and tripartite information. Compare geometry of correlations to 1D chain.

**Results:**

**W-state I3 (n=6):** Uniform I3 = +0.195 across ALL 20 triples. Positive I3 like GHZ (+1.0) but much weaker. Confirms W-state correlations are redundant/classical-like at the tripartite level — information is duplicated across subsystems, not irreducibly distributed.

**Updated archetype table:**
| State   | MI pattern | Discord | I3 sign | I3 magnitude |
|---------|-----------|---------|---------|-------------|
| GHZ     | uniform max | 0     | +1.0    | strong redundant |
| W       | uniform weak | 0.28  | +0.195  | weak redundant |
| Cluster | sparse NN   | 0     | -1.0    | strong irreducible |

**2D Cluster (2x3) vs 1D Chain (6 qubits):**
- 2D has **8 negative I3 triples** (mean -0.4) vs 1D's **6** (mean -0.3)
- 2D geometry creates MORE irreducible correlations — vertical edges add new pathways for 3-body entanglement
- MI is still nearest-neighbor only in both, but "nearest" now includes vertical neighbors

**2D Cluster (2x4) vs 1D Chain (8 qubits):**
- Both have 8 negative I3 triples, mean -0.143
- At larger sizes, the 2D advantage in I3 count diminishes (perhaps because 2x4 is still very "thin")
- I3 minimum is -1.0 in both — the depth of irreducibility is the same

**Key insight:** 2D geometry amplifies the *density* of irreducible correlations at small scales (2x3), but thin rectangles (2x4) don't show a clear advantage over 1D chains. Need squarer geometries (3x3?) to see full 2D effects — but 9 qubits is near our limit.

---

## Experiment 7b: Entanglement Spectrum — 2D vs 1D Cluster

**Goal:** Compute negativity across all bipartitions for 2D cluster (2x3=6 qubits) and compare to 1D chain (6 qubits). Does 2D geometry create a richer spectrum?

**Results:**

**Entanglement spectrum (n=6, 31 bipartitions):**

| Partition size | 1D chain mean | 1D range | 2D grid mean | 2D range |
|---------------|--------------|----------|-------------|----------|
| |A|=1         | 0.500        | [0.5, 0.5] | 0.500     | [0.5, 0.5] |
| |A|=2         | 1.367 ± 0.34 | [0.5, 1.5] | **1.500 ± 0.00** | [1.5, 1.5] |
| |A|=3         | 2.200 ± 1.10 | [0.5, 3.5] | **2.700 ± 0.98** | [1.5, 3.5] |

**Key findings:**
- **2D eliminates the weakest bipartitions.** 1D chain has |A|=2 partitions with negativity as low as 0.5 (taking distant qubits). 2D grid: ALL |A|=2 partitions have negativity 1.5. The vertical edges ensure every pair of qubits has strong entanglement across the cut.
- **2D raises the floor, not the ceiling.** Max negativity is 3.5 in both. But 2D's minimum at |A|=2 is 1.5 vs 1D's 0.5. This is a **3x improvement** in worst-case entanglement.
- **2D spectrum is more varied overall** (std=0.96 vs 0.90) but the variation is concentrated at |A|=3.
- **Connected vs disconnected partitions in 2D** show similar mean negativity (1.98 vs 2.07) — geometry of the partition matters less than in 1D.

**Insight:** 2D geometry creates "entanglement democracy" at the pairwise level — every 2-qubit cut sees the same entanglement. This is why 2D cluster states are the substrate for surface codes: no weak links. 1D chains have "entanglement deserts" for distant qubits.

---

## Experiment 7c: 2D Cluster Under Qubit Loss

**Goal:** Progressive qubit loss from 2D cluster, tracking entropy. Compare corner vs edge vs interior removal. In 1D, position relative to bipartition mattered — does 2D topology create new behavior?

**Results:** *(pending)*
