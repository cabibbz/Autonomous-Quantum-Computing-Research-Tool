# Sprint 003 — Cluster State Entanglement Under Progressive Qubit Loss

**Date:** 2026-03-31
**Status:** Complete

## Idea

Sprint 002c revealed that tracing out one qubit from a cluster state *increases* entanglement entropy (1.0 -> 2.0). This is the most non-classical result so far. Questions:

1. What happens as we trace out qubits one by one? Is there a peak in entanglement?
2. Does it matter WHICH qubit we lose? (Edge vs center of the chain)
3. How does this compare to GHZ and W states under the same progressive loss?
4. Is there a "phase transition" — a critical number of lost qubits where entanglement suddenly collapses?

All experiments at 8 qubits max — well within resource limits.

## Experiments

### Experiment 3a: Progressive Qubit Loss (8 qubits, left-to-right)

Traced out qubits one-by-one from the left edge, measuring half-cut entropy at each step.

| Qubits lost | GHZ | W | Cluster |
|------------|------|-------|---------|
| 0 | 1.000 | 1.000 | 1.000 |
| 1 | 1.000 | 0.954 | 2.000 |
| 2 | 1.000 | 0.954 | 2.000 |
| 3 | 1.000 | 0.811 | 2.000 |
| 4 | 1.000 | 0.811 | 2.000 |
| 5 | 1.000 | 0.544 | 1.000 |
| 6 | 1.000 | 0.544 | 1.000 |

Three completely different regimes:
- **GHZ:** Completely rigid. Entropy = 1.0 no matter how many qubits you lose.
- **W:** Gradual decay. Entanglement bleeds away monotonically.
- **Cluster:** Jumps to 2.0 immediately, holds through 4 lost qubits, then drops. There's a plateau phase.

### Experiment 3b: Position-Dependent Loss (8-qubit cluster)

**Single qubit loss — position matters:**
| Lost qubit | Position | Entropy |
|-----------|----------|---------|
| 0 | edge | 2.000 |
| 1 | near-edge | 2.000 |
| 2 | interior | 2.000 |
| 3 | interior | 1.000 |
| 4 | interior | 1.000 |
| 5 | interior | 1.000 |
| 6 | near-edge | 1.000 |
| 7 | edge | 1.000 |

Sharp transition at the half-cut boundary (between qubits 3 and 4). Losing a qubit from the "left half" (0-2) increases entropy to 2.0, while losing from "right half" (3-7) leaves it at 1.0.

**Two-qubit loss — entropy can reach 3.0:**
Losing pairs (0,2), (0,3), (1,2), or (1,3) pushes entropy to **3 bits** — triple the original. Each lost qubit near the half-cut boundary contributes to this amplification.

Non-adjacent pair loss (mean entropy 1.86) > adjacent pair loss (mean 1.57).

## Analysis

The position dependence reveals something deep: the entropy increase from qubit loss is not about the qubit itself, but about its **relationship to the measurement cut**. Losing a qubit that was "bridging" the bipartition effectively widens the entanglement channel between the two halves.

This is analogous to cutting a knot in a rope: the cut doesn't create tension — it *redistributes* existing connections. In the cluster state, the CZ gates create a web of entanglement. Removing a node near the partition boundary exposes more of that web to the entropy measurement.

The entropy=3.0 cases are particularly interesting: they suggest that the cluster state contains at least 3 ebits (entanglement bits) of hidden entanglement across the bipartition that only become visible when intervening qubits are removed.

## What Surprised Me

1. The sharp boundary at the half-cut: qubits 0-2 give entropy 2.0, qubits 3-7 give 1.0. No gradual transition.
2. Two-qubit loss can give entropy 3.0 — *three times* the original. The entanglement was "trapped" and qubit loss releases it.
3. GHZ's perfect rigidity under loss is itself remarkable — losing 6 of 8 qubits and still having exactly 1 bit of entanglement. The entanglement is stored in a completely delocalized way.

## Next Sprint Ideas

1. Does this position-dependent phenomenon change for 2D cluster states (grid instead of chain)?
2. What about measurement (projective) instead of tracing? Does the basis of measurement matter?
3. Information-theoretic quantities: mutual information, conditional entropy — more nuanced than just entanglement entropy
4. Connect to IIT (Integrated Information Theory) — can we compute Phi for small quantum states?
