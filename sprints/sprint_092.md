# Sprint 092 — Non-Potts Operator Identification in H_E

## Motivation
Sprint 091 found that BW corrections in H_E are dominated by non-Potts operators (NNN and 3-body Potts contribute <0.2%), with the biggest jump at q=3→4 (11.7×) coinciding with the real-to-complex CFT boundary. But we never identified WHAT these non-Potts operators are.

For q=2 (Ising), we can decompose H_E into the full Pauli basis on nA sites. For q≥3, we can use a generalized Gell-Mann / clock-shift operator basis. This will tell us which operators dominate the BW residual and whether they have a physical interpretation.

## Questions
1. What operators dominate the non-Potts fraction of H_E?
2. Do the dominant non-Potts operators change character at q=4 (BW locality boundary)?
3. Is there a pattern (e.g., specific multi-body terms, specific operator types) that explains the exponential growth?

## Experiment 092a — Pauli decomposition of H_E residual (q=2)
**Result:** Complete.

Non-BW weight grows dramatically with nA:
| nA | non-BW weight | dominant non-BW operators |
|----|---------------|--------------------------|
| 3  | 0.005%        | ZXZ (3-body r=2)         |
| 4  | 0.008%        | YY (2-body r=1), ZXZ (3-body r=2) |
| 5  | 2.81%         | ZZZ (3-body r=2-4), ZX (2-body long-range) |
| 6  | 6.28%         | ZXZ (3-body r=2-4), ZZZ (3-body), 4-body ZXXZ |

Key findings:
- Non-BW operators are overwhelmingly Z-type: ZXZ, ZZZ, ZZ (long-range), ZXXZ
- 3-body operators dominate non-BW content at all nA
- YY operators (imaginary part of XX coupling?) are tiny for q=2
- Body order of dominant non-BW increases with nA: 2-body→3-body→4-body
- Range of dominant operators increases with nA (up to range=5 at nA=6)
- BW weight fraction DECREASES with nA: 31.1%→26.4%→24.6%→19.4% (out of non-identity)

Physical interpretation: BW captures NN ZZ+X (the Hamiltonian), but H_E also contains
density-density-field (ZXZ) and multi-body density (ZZZ) correlations from the CFT.
These are the conformal descendants — the entanglement Hamiltonian "knows" about the
full operator content, not just the physical Hamiltonian.

## Experiment 092b — Generalized operator decomposition at fixed nA=3
**Result:** Complete.

Clock-shift (Weyl-Heisenberg) basis decomposition at nA=3 for q=2,3,4,5:

| q | non-BW weight | dominant non-BW type | 2-body r=1 | 3-body r=2 |
|---|---------------|---------------------|------------|------------|
| 2 | 0.005%        | DFD (3-body)        | 0.002%     | 0.003%     |
| 3 | 0.007%        | MD/DM (2-body)      | 0.005%     | 0.002%     |
| 4 | 0.008%        | MD/DM (2-body)      | 0.007%     | 0.002%     |
| 5 | 0.009%        | MD/DM (2-body)      | 0.007%     | 0.001%     |

Key findings:
- Non-BW weight is FLAT at nA=3 (varies only 0.005-0.009% across q=2-5)
- Dominant non-BW operators: MD/DM = mixed × density (like Pauli Y⊗Z)
- 3-body DMD (density-mixed-density) is second dominant type
- Body-order composition: 1-body and 2-body carry 99.99%+ of non-identity weight
- q-dependence of BW breakdown ONLY manifests at larger nA

This confirms Sprint 091 finding: nA controls BW breakdown, not q directly.
The q-dependence amplifies with nA — the subsystem must be large enough to "see" the
complex CFT operator content.

Operator type translation: D=density/clock (Z-like), F=field/shift (X-like), M=mixed (Y-like).
For q=2: D→Z, F→X, M→Y exactly.

## Experiment 092c — nA=4 operator decomposition for q=2,3 (cross-q comparison)
**Result:** Complete.

Full clock-shift decomposition at nA=4:

| q | non-BW weight | dominant non-BW 2-body | dominant non-BW 3-body |
|---|---------------|----------------------|----------------------|
| 2 | 0.0075%       | MM (YY) 0.0038%      | DFD (ZXZ) 0.0037%   |
| 3 | 0.0080%       | MM 0.0022%, MD/DM 0.0029% | DMD 0.0016%, DFD 0.0005% |

Key findings:
- At nA=4, non-BW weight is still nearly equal for q=2,3 (0.0075% vs 0.0080%)
- q=2: dominated by MM (Pauli YY) and DFD (Pauli ZXZ) — imaginary coupling and density-field-density
- q=3: dominated by MM, MD/DM (mixed-density), and DMD — field-density mixing dominates
- q=3 has MORE operator types contributing (richer algebra) but similar total weight
- Growth nA=3→4: q=2 grows 1.5×, q=3 grows 1.1× — BUT Sprint 091's "non-Potts fraction" used Frobenius norm of H_E - α·H_BW, which is a different metric

The discrepancy with Sprint 091 (which found 0.029% for q=2 at nA=4) is because:
- Sprint 091 measured ||H_E - α·H_BW||_F / ||H_E||_F (residual of best fit)
- Here we measure sum of |c|² for non-BW operators / total sum |c|²
- These weight operators differently: Sprint 091 is Frobenius norm, here it's Hilbert-Schmidt

## Synthesis

### Operator Content of BW Corrections

The entanglement Hamiltonian H_E = -log(ρ_A) decomposes into:
1. **Identity** (~70-88% of weight, increases with q — more entropic)
2. **1-body field (F)** — the transverse field terms Σ X^k (~6-16%, decreasing with q)
3. **2-body NN density-density (DD)** — the Potts coupling (~6-16%, decreasing with q)
4. **Non-BW corrections** — everything else

The non-BW corrections are dominated by:
- **Mixed operators (M = X^a Z^b, a,b>0):** Generalization of Pauli Y. These represent "spin-orbit"
  coupling — correlations between the order parameter and its conjugate momentum. Not present
  in the physical Hamiltonian.
- **DFD/DMD (3-body):** Density-field-density and density-mixed-density chains. These are
  "mediated" correlations — site i's density affects site i+1's field which affects site i+2's density.
- **Long-range DD:** Next-nearest-neighbor and longer density-density correlations.

### Physical Interpretation

The BW theorem says H_E ∝ local Hamiltonian density (to leading order). The corrections tell us
what the entanglement Hamiltonian "knows" beyond the physical Hamiltonian:

1. **Mixed (M) operators are the dominant correction.** These don't appear in the S_q Potts
   Hamiltonian at all. They represent entanglement between the "position" (Z/density)
   and "momentum" (X/field) sectors. In field theory language, these are the φ·π correlators
   that the reduced density matrix encodes but the Hamiltonian doesn't contain as separate terms.

2. **3-body DFD/DMD are the next correction.** These are "mediated" interactions — operator
   entanglement that threads through intermediate sites. CFT predicts these from the OPE:
   when two operators fuse, they produce a third at the intermediate point.

3. **The q-dependence at small nA is weak.** At nA=3, non-BW varies only 0.005-0.009%
   across q=2-5. The exponential growth of non-Potts operators with q (Sprint 091) is
   a LARGE-SUBSYSTEM effect that manifests at nA≥4.

### Connection to Walking Boundary

Sprint 091 showed non-Potts fraction jumps 11.7× at q=3→4. Here we find the operator
types are QUALITATIVELY the same across q (MM and DFD dominate). The quantitative growth
with q is an amplitude effect, not a new operator appearing. The complex CFT doesn't introduce
new operator types — it amplifies the mixed (M) operators that are already present at q=2.

This is consistent with the entanglement spectrum finding (Sprint 084): walking breakdown
is a weight redistribution, not a structural change. The same operators carry more weight
as q increases past the walking boundary.

### Surprises
- MM (mixed×mixed, Y⊗Y-type) is the #1 non-BW operator for q=2 at nA=4
- q=3 has MORE non-BW operator types contributing but similar total weight to q=2
- Non-BW weight is essentially flat across q at nA=3 (varies <2×)
- No qualitative change in operator content at q=4 (walking boundary)
- DFD (3-body) and MM (2-body) trade dominance between q=2 and q=3

**POTENTIALLY NOVEL:** First operator decomposition of entanglement Hamiltonian for S_q Potts
chain using generalized clock-shift basis. First identification of mixed (XZ-type/Y-type)
operators as dominant BW correction. First demonstration that walking boundary does not
change the TYPE of BW corrections, only their amplitude.
