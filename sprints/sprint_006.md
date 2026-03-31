# Sprint 006 — Entanglement Monotones: Concurrence & Negativity

**Date:** 2026-03-31
**Goal:** Add proper entanglement monotones (concurrence, negativity) to the archetype classification. These quantify *how much* entanglement exists in each pair, complementing MI (total correlations), I3 (structure), and discord (classical vs quantum nature).

**Background:**
- **Concurrence** (Wootters 1998): defined for 2-qubit states. C=0 means separable, C=1 means maximally entangled. Gold standard for bipartite entanglement.
- **Negativity** (Vidal & Werner 2002): sum of negative eigenvalues of partial transpose. Works for any bipartition. N=0 means PPT (no distillable entanglement), N=0.5 means maximally entangled pair.

**Prediction:** Based on what we know:
- GHZ pairwise RDMs are classical mixtures (discord=0), so concurrence should be 0 — all entanglement is multipartite
- W state has nonzero discord, so should have nonzero concurrence — entanglement is distributed in pairs
- Cluster state has zero pairwise discord but strong I3 — unclear! Nearest-neighbor pairs might have some entanglement

## Experiment 6a: Pairwise Concurrence & Negativity

**Status:** Running...

## Experiment 6b: Negativity Under Noise

**Status:** Pending
