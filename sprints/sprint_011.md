# Sprint 011 — Entanglement Dynamics: How Entanglement Builds Gate by Gate

**Date:** 2026-03-31
**Status:** Complete (3/3 experiments)
**Goal:** Track how entanglement measures evolve during circuit construction — do they build gradually or exhibit phase transitions?

## Motivation

Sprints 001–010 characterized entanglement in *final* states. But how does entanglement *emerge*? Is it gradual or sudden? Do different states follow different trajectories through information space? Does the order of gates matter? These questions connect to:
- Quantum circuit complexity (how many gates to reach a given entanglement level?)
- Quantum phase transitions (sudden changes in entanglement structure)
- Circuit optimization (are some construction paths more efficient?)

## Experiments

### 11a: Entanglement Trajectory — Entropy & Negativity Layer by Layer (n=6)

Tracked half-cut entanglement entropy and logarithmic negativity at each gate layer.

| State | Layers | Entropy trajectory | Behavior |
|-------|--------|--------------------|----------|
| GHZ | 7 | 0,0,0,0,**1**,1,1 | Step function at first cross-cut gate |
| W | 7 | 0,0,0,0,**1**,1,1 | Identical step function |
| Cluster 1D | 7 | 0,0,0,0,**1**,1,1 | Identical step function |
| Cluster 2D | 9 | 0,0,0,0,0,0,**1,2,3** | Linear staircase — each vertical edge adds exactly 1.0 |

**Key finding:** Half-cut entropy is a step function: exactly zero until the first gate crosses the bipartition, then jumps to full value. This is trivially determined by the circuit topology relative to the partition. **Exception:** 2D cluster shows linear accumulation — each cross-cut CZ gate adds exactly 1 bit, consistent with the additive entanglement structure found in Sprint 007.

### 11b: MI Dynamics — Pairwise Correlations Layer by Layer

Tracked total pairwise mutual information, MI sparsity, and number of significant pairs at each layer.

| State | MI trajectory | Final total MI | Monotonic? |
|-------|--------------|----------------|------------|
| GHZ | 0→2→3→6→10→15 | 15.0 (all pairs) | **Yes** — each gate adds new correlations |
| W | 0→1.3→2.2→3.3→4.6→5.7 | 5.7 (all pairs) | **Yes** — gradual accumulation |
| Cluster 1D | 0→2→3→**2**→2→2 | 2.0 (2 edge pairs only) | **NO** — MI drops at CZ(2,3) |
| Cluster 2D | 0→2→3→5→6→**2**→**2**→**0** | 0.0 (no pairs!) | **NO** — all MI destroyed |

**Key finding: Cluster states exhibit NON-MONOTONIC mutual information.** Adding entangling gates *destroys* pairwise correlations — MI is consumed to create irreducible multipartite structure. The 2D cluster is most dramatic: builds up to MI=6.0 then collapses to ZERO. All information moves to 3+ body correlations.

**GHZ and W are opposite:** MI grows monotonically, with GHZ accumulating faster (all pairs reach MI=1.0) and W distributing more evenly (all pairs at MI=0.38).

The non-monotonic behavior corresponds exactly to when CZ gates connect previously independent entangled clusters — merging two correlated pairs into a multipartite structure destroys the pairwise correlations.

### 11c: Construction Order Dependence — Does the Path Matter?

Tested 4 different CZ orderings for 1D cluster state and 4 CNOT orderings for GHZ, all producing identical final states.

**GHZ MI trajectories (all orderings identical):**
```
forward:     0 → 2 → 3 → 6 → 10 → 15
reverse:     0 → 2 → 3 → 6 → 10 → 15
outside_in:  0 → 2 → 3 → 6 → 10 → 15
alternating: 0 → 2 → 3 → 6 → 10 → 15
```

**Cluster 1D MI trajectories (wildly different):**
```
forward:     0 → 2 → 3 → 2 → 2 → 2    (peak: 3.0)
reverse:     0 → 2 → 3 → 2 → 2 → 2    (peak: 3.0)
outside_in:  0 → 2 → 4 → 5 → 6 → 2    (peak: 6.0!)
inside_out:  0 → 2 → 3 → 2 → 2 → 2    (peak: 3.0)
```

**Key finding:** GHZ construction is completely **path-independent** — all orderings produce identical MI trajectories. This reflects GHZ's permutation symmetry: the MI at any stage depends only on *how many* qubits are entangled, not *which ones*.

Cluster construction is wildly **path-dependent**. The outside-in ordering builds up **3x the final MI** (peak 6.0 vs final 2.0) by creating isolated entangled pairs that accumulate before the connecting gate destroys them. The last gate CZ(2,3) **destroys 67% of total MI** in a single step.

**Half-cut entropy also path-dependent for cluster:**
- Inside-out: entropy appears at layer 1 (first gate crosses cut)
- Outside-in: entropy appears only at layer 5 (last gate bridges halves)
- All reach the same final entropy = 1.0

## Key Insights

1. **Entangling gates can destroy correlations.** This is counterintuitive — a CNOT or CZ is an entangling gate, yet for cluster states, applying CZ across a boundary *reduces* total pairwise MI. The mechanism: pairwise correlations are "consumed" to create irreducible 3+ body correlations. This is the information-theoretic signature of genuine multipartite entanglement being created.

2. **GHZ is path-independent, cluster is path-dependent.** GHZ has full permutation symmetry — it doesn't matter which qubits you entangle in which order. Cluster states are *geometric* — the spatial structure of the circuit matters enormously. This is the constructional signature of the topological vs democratic distinction found in Sprints 006–007.

3. **2D cluster construction completely annihilates pairwise MI.** The final state has zero MI for ALL 15 qubit pairs. Every bit of pairwise correlation is converted into higher-order structure. This is the most dramatic demonstration of the "body-order hierarchy" from Sprint 010 — 2D cluster lives entirely at the 3+ body level.

4. **Construction-time "phase transitions."** The outside-in cluster ordering shows a dramatic discontinuity: MI drops from 6.0 to 2.0 in a single gate (67% destruction). This is a construction-time analog of a quantum phase transition — a sudden, qualitative change in information structure triggered by a single local operation.

5. **Connection to quantum computation.** In measurement-based quantum computation, single-qubit measurements on cluster states implement gates. Our finding that single CZ gates can destroy 67% of pairwise MI is the time-reverse of this: single measurements on the final state can "release" the locked-up multipartite correlations back into pairwise form. Construction and computation are duals.

## Next Ideas

- **Entanglement speed limits:** What's the maximum rate of entanglement growth (bits per gate)? Is there a Lieb-Robinson-like bound?
- **Circuit complexity from entanglement:** Can we define circuit complexity as the minimum number of gates to reach a target MI/I3 profile?
- **Real hardware comparison:** Run these circuits on IBM hardware — does real noise create non-monotonic MI in states that show monotonic MI in simulation?
- **Scrambling dynamics:** Track operator spreading and out-of-time-ordered correlators (OTOCs) during construction
- **Random circuit dynamics:** What trajectory does a random circuit take through MI space? Is there a universal behavior?
