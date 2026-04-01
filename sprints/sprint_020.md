# Sprint 020 — Topological Codes: The Toric Code's Entanglement Structure

**Date:** 2026-03-31
**Status:** Complete (3/3 experiments)

## Motivation

The toric code has been the top "next" item since Sprint 007 when we discovered that 2D cluster states show topological entanglement protection. Sprints 014-019 explored algebraic QEC codes ([[5,1,3]], Steane, Shor) and their information structure. Now we ask: how does a genuinely *topological* code differ?

The toric code on a 2×2 torus encodes 2 logical qubits in 8 physical qubits (distance 2). It's the simplest code with topological order — its logical operators are non-contractible loops on the torus, and no local measurement can distinguish the 4 ground states.

## Experiments

### 20a: Toric Code Entanglement Structure

**Goal:** Prepare [[8,2,2]] toric code ground state, compute MI, I3, negativity. Compare with [[5,1,3]].

**Results:**

| Measure | Toric [[8,2,2]] | [[5,1,3]] |
|---------|-----------------|-----------|
| Total pairwise MI | 4.0 (4 pairs × 1.0) | 0.0 |
| Nonzero MI pairs | 4/28 | 0/10 |
| I3 range | [0.0, 0.0] | [-1.0, -1.0] |
| Negative I3 triples | 0/56 | 10/10 |
| Half-cut entropy | 1.0 | 2.0 |
| Single-qubit entropy | 1.0 (uniform) | 1.0 (uniform) |
| Negativity |A|=1 | 0.5 (uniform) | 0.5 (uniform) |
| Negativity |A|=2 | 0.5–1.5 (non-uniform) | 1.5 (uniform) |

**Key finding:** The toric code has the OPPOSITE information structure from [[5,1,3]]:
- **Toric: all pairwise, no multipartite** — MI=1.0 for 4 specific pairs (non-contractible loop edges), I3=0 everywhere
- **[[5,1,3]]: all multipartite, no pairwise** — MI=0 everywhere, I3=-1.0 for all triples

The 4 MI-correlated pairs are: (0,1), (2,3), (4,6), (5,7) — edges related by translation on the torus, forming the non-contractible loops that define the logical operators.

### 20b: Topological Degeneracy and Local Indistinguishability

**Goal:** Prepare all 4 ground states, verify local indistinguishability — the hallmark of topological order.

**Results:**
- 4 ground states: |++⟩, |+−⟩, |−+⟩, |−−⟩ (eigenvalues of Z_L1 = Z_q0 Z_q1, Z_L2 = Z_q4 Z_q6)
- **Single-qubit level:** ALL 4 states have identical RDMs. Trace distance = 0.000000 for every qubit. Perfect indistinguishability.
- **Two-qubit level:** 4 of 28 pairs can distinguish states (trace distance = 1.0): exactly (0,1), (2,3), (4,6), (5,7) — the non-contractible loop pairs.
- **All 4 states have identical entropy** at every scale: half-cut S = 1.0, single-qubit S = 1.0.

**Key finding:** The toric code's code distance d=2 manifests as perfect 1-qubit indistinguishability with distinguishability at exactly 2 qubits — and ONLY for the 4 specific pairs forming non-contractible loops. This is the operational definition of topological order: logical information is encoded in the global topology, invisible to any local observation. The distinguishing pairs ARE the logical operators.

**Note on TEE:** The Kitaev-Preskill topological entanglement entropy gives TEE=0 for all tripartitions on the 2×2 torus. This is expected — the 2×2 lattice is too small for the area-law/topological separation to work. The KP construction requires regions much larger than the correlation length, which isn't achievable with 8 qubits.

### 20c: Page Curve and Noise Fingerprint

**Goal:** Page curve for logical information recovery; noise response comparison.

**Page Curve Results:**

| Subset size |A| | Toric [[8,2,2]] (mean ± std) | [[5,1,3]] |
|-------------|-------------------------------|-----------|
| 1 | 0.00 ± 0.00 | 0.00 |
| 2 | 0.14 ± 0.35 [0, 1.0] | 0.00 |
| 3 | 0.43 ± 0.49 [0, 1.0] | 1.00 |
| 4 | 0.77 ± 0.42 [0, 1.0] | 1.00 |
| 5 | 1.00 ± 0.00 | — |

**Key finding:** Two fundamentally different Page curve shapes:
- **[[5,1,3]]: knife-edge** — 0 for |A|<3, 1.0 for |A|≥3. All subsets of a given size behave identically (std=0). Matches Singleton bound exactly (n-d+1=3).
- **Toric: gradual + bimodal** — at each size, subsets either recover everything (Holevo=1.0) or nothing (Holevo=0.0), but the *fraction* that recover increases gradually with size. Full recovery guaranteed at |A|=5 (not at n-d+1=7, because we're recovering 1 of 2 logical qubits).

**Noise Fingerprint:**

| Noise type | Toric break-even | [[5,1,3]] break-even | Note |
|------------|-----------------|---------------------|------|
| Depolarizing | p≈0.30 | p≈0.30 | Similar |
| Amplitude damping | γ≈0.50 | γ≈0.50 | Similar |
| Phase damping | **NEVER** | **NEVER** | Both immune! |

**Key finding:** Both codes show PERFECT Holevo preservation (=1.0) under phase damping at ALL noise levels. This is because the logical states are Z-eigenstates, and phase damping (Z-diagonal noise) commutes with the logical Z operator. The fidelity degrades, but the distinguishability of logical states is unchanged. To see phase damping vulnerability, one would need to encode logical X information instead.

Under depolarizing and amplitude damping, the toric code has slightly LOWER fidelity than [[5,1,3]] (8 qubits = 8 noise targets vs 5), but identical Holevo break-even points. More qubits means more noise but also more information spread — these effects cancel for logical information.

## Synthesis: A New Code Archetype

The toric code represents a fourth QEC archetype, distinct from all three algebraic archetypes discovered in Sprint 015:

| Property | Democratic ([[5,1,3]]) | Selective (Steane) | Hierarchical (Shor) | **Topological (Toric)** |
|----------|----------------------|-------------------|---------------------|----------------------|
| Pairwise MI | 0.0 | 0.0 | 9.0 (within blocks) | **4.0 (loop pairs)** |
| I3 | -1.0 (all triples) | -1.0 (7/35) | -1.0 (cross-block) | **0.0 (all triples)** |
| Page curve | Knife-edge | Bimodal | Staircase | **Gradual+bimodal** |
| Negativity | Uniform | Geometry-dep | Block structure | **Loop-dependent** |
| Information topology | Spread everywhere | Hamming geometry | Concatenated layers | **Non-contractible loops** |

**The toric code is the only code with ZERO I3.** Its correlations are entirely pairwise — but specifically along the non-contractible loops of the torus. This is the information-theoretic signature of topological order: information is stored in the *topology* (non-contractible paths), not in local correlations or multipartite structure.

This explains why topological codes are so different from algebraic codes in practice:
- Algebraic codes (esp. [[5,1,3]]) hide information in irreducible multipartite correlations (I3 < 0) — no subset smaller than d reveals anything
- The toric code hides information in topological structure — the pairwise correlations exist and are maximal, but they're along non-local paths that wrap the torus

## Next Sprint Ideas
1. **Toric code under local noise with syndrome extraction** — does topological protection manifest differently from algebraic codes when using actual syndrome-based correction?
2. **Surface code (planar boundary)** — how does removing periodic boundary conditions change the information structure?
3. **Combined T1+T2 realistic noise model** — the real-hardware noise regime
4. **Logical X information under phase damping** — to see the vulnerability we're immune to in Z-basis
5. **Scaling: [[18,2,3]] toric code** (3×3 torus) — if computationally feasible (9 qubits for density matrix)
