# Sprint 010 — GME Witnesses & Efficient Entanglement Detection

**Date:** 2026-03-31
**Status:** Complete (3/3 experiments)
**Goal:** Can we detect entanglement class from a few Pauli measurements instead of full state tomography?

## Motivation

Sprints 001–009 characterized entanglement using full density matrix access (entropy, MI, I3, discord, concurrence, negativity, Phi, noise fingerprints). But full tomography scales as 3^n measurements — exponentially expensive. Entanglement witnesses and Pauli expectation profiles offer a practical path to entanglement classification.

## Experiments

### 10a: Pauli Expectation Value Profiles (n=6)

Measured all 1-qubit (18) and 2-qubit (135) Pauli expectations for each archetype.

| State | 1-body nonzero | 2-body nonzero | Signature |
|-------|---------------|----------------|-----------|
| GHZ | 0/18 (100% sparse) | 15/135 (ZZ=1 all pairs) | All-pairs ZZ correlation |
| W | 6/18 (Z≠0) | 45/135 (XX=YY=ZZ=1/3) | Uniform weak correlations |
| Cluster 1D | 0/18 | 2/135 (neighbor XZ only) | Stabilizer fingerprint |
| Cluster 2D | 0/18 | **0/135 (100% sparse!)** | Invisible at 2-body level |

**Key finding:** 2D cluster states have ZERO 1-body and 2-body Pauli expectations. All correlations are genuinely 3+ body. This is a qualitative distinction from every other state we've studied.

### 10b: Fidelity-Based GME Witnesses

Constructed entanglement witnesses W = α·I - |ψ⟩⟨ψ| where α = max biseparable fidelity.

| State | Bisep Bound α | GME Margin | GME Death (noise p) |
|-------|--------------|------------|---------------------|
| GHZ | 0.500 | 0.500 | p = 0.55 |
| W | 0.833 | 0.167 | p = 0.20 |
| Cluster 1D | 0.500 | 0.500 | p = 0.55 |
| Cluster 2D | 0.500 | 0.500 | p = 0.55 |

**Key findings:**
- W has the weakest GME witness (margin only 0.167) and GME dies earliest under noise (p=0.20)
- GHZ and both cluster states have identical bisep bound (0.500) and noise thresholds
- Cross-detection is perfectly diagonal: no witness for state A detects GME in state B
- States are nearly mutually orthogonal (cross-fidelities ≈ 0)

### 10c: Minimal Pauli Classifier

Built a decision tree classifier using only 4 Pauli measurement settings:

```
1. <Z> on any qubit     → nonzero? → W state
2. <ZZ> on any pair      → ≈1?     → GHZ state
3. <XZ> on adjacent pair → ≈±1?    → Cluster 1D
4. <XZZ> (2D stabilizer) → ≈1?     → Cluster 2D
5. All ≈ 0                         → Product state
```

**Results:**
- **100% accuracy** on all 5 state classes (pure states)
- Robust to **40% depolarizing noise** (100% accuracy at p=0.4)
- Only fails at p≥0.5 (maximally noisy)
- **182x compression** vs full tomography (4 settings vs 728)

**Stabilizer witness:** 1D cluster stabilizer sum = 6.0 (perfect), exactly 0.0 on all other states.

## Key Insights

1. **The body-order hierarchy is a classifier itself.** GHZ is a 2-body state (all information in ZZ). W is a 1-body state (detectable from single-qubit Z). 1D Cluster is a 2-body stabilizer state (XZ). 2D Cluster is a **3-body** state (XZZ). Each archetype lives at a different level of the correlation hierarchy.

2. **2D cluster's invisibility at the 2-body level** is not a deficiency — it's what makes it computationally powerful. In MBQC, local measurements extract information, but no local or pairwise measurement reveals the computation's content. This is information-theoretic security built into the physics.

3. **4 measurements suffice for classification** — but only if you know which 4. The correct measurements are the **stabilizer generators** of each target class. This connects witness theory to stabilizer formalism: witnessing ≡ stabilizer verification.

4. **W state is the most fragile to detect** despite being the most robust to noise (Sprint 009). Its GME margin is only 0.167 because it has high overlap with biseparable states (α=0.833). Robustness to noise ≠ ease of detection.

## Summary Table (All 10 Sprints)

| Property | GHZ | W | Cluster 1D | Cluster 2D |
|----------|-----|---|------------|------------|
| Body-order | 2 | 1 | 2 | 3+ |
| MI structure | all-pairs | all-pairs (weak) | nearest-neighbor | nearest-neighbor |
| I3 | +1 (redundant) | +0.2 (weak red.) | -1 (irreducible) | more negative |
| Discord | 0 | 0.28 (74% quantum) | 0 | 0 |
| Concurrence | 0 (all multipartite) | 2/n (all pairwise) | 0 | 0 |
| Neg. spectrum | flat | growing | explosive | uniform |
| Phi retention (loss) | 50% | 71% | 0-100% | 100% |
| Best noise channel | dephasing | all similar | amplitude | amplitude |
| Worst noise channel | amplitude | (none specific) | dephasing | dephasing |
| GME margin | 0.500 | 0.167 | 0.500 | 0.500 |
| Min measurements | 1 (ZZ) | 1 (Z) | 1 (XZ) | 1 (XZZ) |
