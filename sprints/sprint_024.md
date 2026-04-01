# Sprint 024 — Surface Code vs Shor: Same [[9,1,3]], Different Architecture

**Date:** 2026-03-31
**Status:** Complete (3/3 experiments)

## Motivation

Sprint 023 showed that the figure of merit determines the winner: avg Holevo → Shor wins, min-basis Holevo → [[5,1,3]] wins. But we hadn't tested the most practically important code family: the **rotated surface code**.

The rotated surface code at distance d=3 is [[9,1,3]] — the **exact same parameters** as the Shor code. Same n, same k, same d. But completely different architecture:
- **Shor:** Concatenated (phase-flip inside bit-flip)
- **Surface:** Topological (2D lattice with boundary stabilizers)

**Literature search:** No prior work directly compares surface [[9,1,3]] vs Shor [[9,1,3]] using basis-averaged Holevo. Papers on surface codes focus on large-scale logical error rates, not small-scale information-theoretic comparison.

## Results

### 24a: Entanglement Structure — The Boundary Fingerprint

Built the d=3 rotated surface code and compared its information topology:

| Metric | Surface-9 | Shor-9 | [[5,1,3]] |
|--------|-----------|--------|-----------|
| Total pairwise MI | 4.0 | 9.0 | 0.0 |
| Non-zero MI pairs | 4/36 | 9/36 | 0/10 |
| Negative I3 triples | 8/84 | 0/84 | 10/10 |
| Neg |A|=1 (range) | 0.5-0.5 | 0.5-0.5 | 0.5-0.5 |
| Neg |A|=4 (range) | 0.5-7.5 | 0.5-3.5 | — |

**Key finding:** The 4 MI-correlated pairs are exactly the weight-2 boundary stabilizers: (0,1), (2,5), (3,6), (7,8). The surface code is a **hybrid** of the toric code (boundary MI) and [[5,1,3]] (irreducible I3):
- Toric [[8,2,2]]: MI=4.0, I3=0 everywhere (Sprint 020)
- Surface [[9,1,3]]: MI=4.0, I3=-1.0 for 8/84 triples
- [[5,1,3]]: MI=0.0, I3=-1.0 for ALL triples

The surface code stores information at its boundaries — the bulk is "invisible" like 2D cluster states, but boundaries leak correlations.

### 24b: Surface vs Shor Under Combined T1+T2 Noise

Swept 5×5 (γ,λ) grid with basis-averaged and min-basis Holevo:

**Average Holevo winners:** Shor 72%, [[5,1,3]] 28%, Surface 0%
**Min-basis Holevo winners:** [[5,1,3]] 84%, Shor 12%, Uncoded 4%, Surface 0%

Surface **never wins** on either metric! Head-to-head vs Shor:
- Surface wins 6/25 on avg Holevo (only at low γ, high λ)
- Surface wins 10/25 on min-basis Holevo (high dephasing corner)

Basis asymmetry: [[5,1,3]]=0.016, Surface=0.22, Shor=0.31

**The surface code at d=3 is a "jack of all trades, master of none."**

### 24c: Why The Surface Code Loses — Structural Analysis

**Per-channel Holevo at p=0.10:**

| Noise | Code | Z | X | Y | Asymmetry |
|-------|------|---|---|---|-----------|
| Bit-flip | [[5,1,3]] | 0.957 | 1.000 | 0.957 | 0.043 |
| Bit-flip | Shor | 1.000 | 0.655 | 0.655 | 0.345 |
| Bit-flip | Surface | 0.551 | 1.000 | 0.551 | 0.450 |
| Phase-flip | [[5,1,3]] | 1.000 | 0.957 | 0.957 | 0.043 |
| Phase-flip | Shor | 0.464 | 1.000 | 0.464 | 0.536 |
| Phase-flip | Surface | 1.000 | 0.551 | 0.551 | 0.450 |
| Y-flip | ALL codes | ~0.997 | ~0.997 | 1.000 | ~0.005 |
| Depol | [[5,1,3]] | 0.848 | 0.848 | 0.848 | 0.000 |
| Depol | Surface | 0.872 | 0.872 | 0.787 | 0.085 |

**Root cause: Non-uniform stabilizer participation.**

| Qubit | Stabilizers | Role |
|-------|-------------|------|
| 0, 2, 6, 8 | 2 | Corners — least protected |
| 1, 3, 5, 7 | 3 | Edges |
| 4 | 4 | Center — most protected |

[[5,1,3]]: all 5 qubits in exactly 4 stabilizers → **perfect symmetry**.
Surface-9: 2-4 stabilizers per qubit → **non-uniform protection**.

**Surprise: Both surface and Shor are near-perfect under Y-errors (asymmetry 0.005).** Y-flip is the one error type where all codes converge. This is because Y = iXZ, so Y-errors simultaneously flip bit and phase — the concatenation and topological structures both handle this combination efficiently.

## Key Insight

**The surface code's advantage is architectural, not small-scale.** At d=3, the surface code has the same distance as [[5,1,3]] but uses 9 qubits instead of 5 and has non-uniform protection from boundary effects. Its advantage emerges only at larger d:
1. **d scales as O(√n)** for surface codes vs O(1) for fixed codes — the growth rate is the advantage
2. **Local stabilizers** enable practical syndrome extraction (no long-range entanglement needed)
3. **2D layout** matches hardware connectivity (superconducting qubit grids)

At small scale, [[5,1,3]]'s perfect algebraic symmetry (every qubit in exactly 4 stabilizers, isotropic distance) dominates. The surface code pays the "boundary tax" — 4 weight-2 stabilizers create weak points that break isotropy. This tax becomes negligible at large d where boundary qubits are a vanishing fraction.

**The boundary paradox:** The surface code's MI-leaking boundaries are simultaneously its greatest weakness (information exposure) and its greatest engineering advantage (they're where syndrome measurements happen and logical operators terminate).

## Architecture Summary

| Property | [[5,1,3]] | Shor-9 | Surface-9 |
|----------|-----------|--------|-----------|
| Structure | Algebraic | Concatenated | Topological |
| MI pattern | Zero everywhere | Within-block | Boundary pairs |
| I3 pattern | All negative | All zero | Mixed (8/84 neg) |
| Isotropy | Perfect (0.016) | Worst (0.31) | Intermediate (0.22) |
| Avg Holevo winner | 28% | 72% | 0% |
| Min-basis winner | 84% | 12% | 0% |
| Scaling advantage | None (fixed n=5) | Concatenation | O(√n) distance |

## Next Steps
- Real hardware test (QPU): test [[5,1,3]] vs surface code fidelity under real noise
- Surface code at d=5 (25 qubits — beyond our simulator limit, but could test on QPU)
- Syndrome extraction circuits: how does measurement add noise?
- XZZX surface code variant (optimized for biased noise)
