# Sprint 096 — BW Threshold Mechanism: Non-Potts Operators, Not Spectrum

## Hypothesis
BW breaks down at nA=6 for q=2 (Sprint 095). Is this driven by entanglement spectrum tail weight crossing a critical value, or by operator structure in H_E?

## Experiments

### 096a — Entanglement spectrum vs nA at fixed n (q=2)
q=2 periodic, n=14, nA=3-7.

| nA | nA/n  | tail_weight | %S(tail) | ent_gap | 1-R²     |
|----|-------|-------------|----------|---------|----------|
| 3  | 0.214 | 1.93e-2     | 0.137    | 1.149   | 1.72e-4  |
| 4  | 0.286 | 2.56e-2     | 0.163    | 1.091   | 2.64e-4  |
| 5  | 0.357 | 2.99e-2     | 0.179    | 1.057   | 8.99e-4  |
| 6  | 0.429 | 3.24e-2     | 0.188    | 1.039   | **1.91e-1** |
| 7  | 0.500 | 3.32e-2     | 0.191    | 1.034   | **3.96e-1** |

**Finding:** Entanglement spectrum is SMOOTH across the BW threshold. Tail weight changes only 8% from nA=5→6, but 1-R² jumps 212×. Entanglement gap, participation ratio, number of significant eigenvalues — all smooth. Spectrum does NOT predict BW breakdown.

### 096b — Cross-q confirmation (q=3, q=5)
Same smooth-spectrum-sharp-BW pattern confirmed:

| q | nA | nA/n  | tail_w   | 1-R²     |
|---|-----|-------|----------|----------|
| 3 | 3   | 0.300 | 2.81e-2  | 1.03e-3  |
| 3 | 4   | 0.400 | 3.50e-2  | 9.55e-4  |
| 3 | 5   | 0.500 | 3.72e-2  | **2.11e-1** |
| 5 | 3   | 0.375 | 3.67e-2  | 2.11e-3  |
| 5 | 4   | 0.500 | 4.12e-2  | **4.56e-2** |

Universal: spectrum smooth, BW sharp, for all q tested.

### 096c — Operator range decomposition of BW residual (q=2 n=14)
Decomposed H_res = H_E - α·H_BW into Pauli operators classified by range.

**BW residual dominated by MAXIMUM-RANGE operators:**
| nA | 1-R²    | r=0  | r=1  | r=2  | r≥(nA-2) |
|----|---------|------|------|------|-----------|
| 3  | 1.7e-4  | 0.43 | 0.24 | 0.33 | 0.33      |
| 4  | 2.6e-4  | 0.23 | 0.49 | 0.28 | 0.28      |
| 5  | 8.6e-4  | 0.12 | 0.24 | 0.20 | **0.64**  |
| 6  | 1.7e-1  | 0.02 | 0.03 | 0.06 | **0.96**  |
| 7  | 3.8e-1  | 0.00 | 0.04 | 0.14 | **0.96**  |

At the threshold, 96% of the residual is in operators spanning ≥(nA-2) sites.

Top operators at nA=6: ZIIIXI (r=4, coeff=0.80), ZIIXXI (r=4, coeff=0.69), IZXZII (r=2, coeff=0.60). These are NOT Potts-type operators — they're mixed (clock × shift) combinations.

### 096d — Augmented BW with free Potts operators (q=2 n=14)
Fit H_E using ALL Potts-type operators (delta bonds + shift fields) at all ranges, with free coefficients.

| nA | BW(1-param) | free-NN | free-all-Potts |
|----|-------------|---------|----------------|
| 3  | 1.72e-4     | 9.73e-5 | 9.73e-5        |
| 4  | 2.64e-4     | 1.46e-4 | 1.46e-4        |
| 5  | 9.20e-4     | 7.22e-4 | 6.67e-4        |
| 6  | **1.71e-1** | 1.67e-1 | **1.63e-1**    |
| 7  | **3.61e-1** | 3.59e-1 | **3.47e-1**    |

**Adding longer-range Potts operators barely helps.** At nA=6, going from BW (1 param) to free-fit with 21 Potts operators only improves 1-R² from 0.171 to 0.163 (5% improvement). **The threshold is entirely non-Potts operators, not operator range.**

Non-Potts content (1-R² of best Potts fit) vs nA:
- nA=3: 9.7e-5
- nA=4: 1.5e-4
- nA=5: 6.7e-4
- **nA=6: 0.163** (244× jump from nA=5)
- nA=7: 0.347

The non-Potts fraction itself has a **threshold at nA=6** — identical to the BW threshold.

## Mechanism (Complete Picture)

**BW breaks down when the entanglement Hamiltonian develops significant weight on operators outside the physical Hamiltonian's operator algebra.**

1. The physical Hamiltonian H has only Potts (delta-bond) and shift (field) operators
2. H_E develops mixed operators (XZ-type, "momentum×position" correlators) that grow with nA
3. For nA ≤ 5 (q=2): non-Potts weight < 0.07% — BW works (captures 99.9%+ of H_E)
4. At nA = 6: non-Potts weight jumps to 16.3% — BW becomes inadequate
5. This is a UV lattice effect (Sprint 095): corrections are built from lattice operators absent from the continuum Hamiltonian
6. The entanglement spectrum (eigenvalues of ρ_A) is smooth — the threshold is in the EIGENVECTORS (operator content), not the eigenvalues

This connects Sprint 092 (non-Potts operators = mixed clock-shift type, grow as exp(1.6q)) to Sprint 095 (BW threshold is UV lattice effect). The BW envelope shape is essentially optimal for the Potts subspace — the error comes entirely from the operator algebra mismatch.

## Surprises
- Entanglement spectrum is SMOOTH across BW threshold — spectrum doesn't predict breakdown
- Adding all Potts operators (any range) to BW barely helps: 5% improvement at threshold
- BW envelope is nearly optimal for Potts subspace (free-NN only 2× better than BW at small nA)
- Non-Potts fraction itself has a sharp threshold at nA=6 (244× jump from nA=5)
- The threshold is in EIGENVECTORS of ρ_A, not eigenvalues
- At nA=6, max-range operators carry 96% of residual vs 33% at nA=3

## Implications
- BW accuracy is fundamentally limited by operator algebra completeness, not by entanglement structure
- To improve H_E approximation beyond BW, one needs to add non-Potts (mixed) operators — not longer-range Potts operators
- The non-Potts threshold at nA ≈ 5-6 for q=2 is likely ≈ 2ξ_lattice (twice the lattice correlation length)
- Walking (q=5) shifts threshold down because walking enhances the mixed operators at smaller nA

[Full data: results/sprint_096*.json]
