# Sprint 034 — BW Operator Algebra: Predicting Locality from Symmetry Constraints

**Date:** 2026-04-01
**Status:** In Progress

## Idea

Sprint 033 overturned the prediction that larger symmetry groups give better BW locality: S₃ (order 6, d=3) gave 76.5%, LOWER than Z₂ (order 2, d=2) at 91%. The proposed correction: BW accuracy = dim(Hamiltonian terms) / dim(symmetry-invariant terms), depending on BOTH group G AND local dimension d.

This sprint tests this quantitatively:
1. **Operator counting**: Explicitly compute dimensions of G-invariant operator spaces for different (G, d) combinations
2. **Z₃ clock model**: Same d=3 as Potts but SMALLER group (order 3 vs 6) — isolates group size effect at fixed d
3. **Prediction verification**: Does the ratio quantitatively predict all BW data from Sprints 032-033?

**Literature:** Giudici et al. (2018) confirmed BW for Potts qualitatively. No prior work systematically counts symmetry-invariant operator dimensions to predict BW accuracy. The operator counting approach is new.

**Predictions:**
- Z₃(d=3) should have MORE invariant operators than S₃(d=3) since smaller group constrains less → LOWER BW locality than S₃'s 76.5%
- The ratio dim(H_terms)/dim(G-inv) should quantitatively predict measured BW locality across all models

## Experiment 34a: Operator Counting

**Result:** Explicit G-invariant operator dimensions computed for all four symmetry groups.

### Key data (2-body, n_A=2-4):

| Group | d | |G| | n_A | Total tl ops | G-inv tl ops | G-inv/total | H terms | H/G-inv |
|-------|---|-----|-----|-------------|-------------|-------------|---------|---------|
| Z₂ | 2 | 2 | 4 | 255 | 71 | 27.8% | 7 | 0.099 |
| U(1) | 2 | ∞ | 4 | 255 | 42 | 16.5% | 6 | 0.143 |
| S₃ | 3 | 6 | 3 | 728 | 69 | 9.5% | 5 | 0.072 |
| Z₃ | 3 | 3 | 3 | 728 | 125 | 17.2% | 5 | 0.040 |

### Analysis:
1. **H/G-inv ratio is monotonically correlated with measured BW locality:** U(1) (0.143) → 100%, Z₂ (0.099) → 91%, S₃ (0.072) → 76.5%
2. **At fixed d=3:** S₃ (order 6) constrains MORE than Z₃ (order 3): 69 vs 125 invariant operators. So S₃ should have BETTER BW → confirmed by data
3. **Across d:** d=3 models have far more invariant operators than d=2 models at same n_A, overwhelming the group-size advantage
4. **Z₃ prediction:** H/G-inv = 0.040 (lowest of all four) → BW locality should be LOWEST, below S₃'s 76.5%

### Surprises:
- At single-site level, all groups within same d give identical counts (1 invariant operator each). Discrimination appears only at 2+ sites.
- U(1) has FEWER invariant operators than Z₂ at 2-site level (4 vs 5) despite being an infinite group — continuous symmetry constrains more effectively than any finite subgroup

## Experiment 34b: Z₃ Clock Model BW Locality

**Result:** Peak BW locality = **76.5%** at h/J=0.75 — identical to S₃ Potts!

### Critical discovery: Clock ≡ Potts at q=3

The Z₃ clock interaction Z_i Z_{i+1}† + h.c. = diag(2,-1,-1,-1,2,-1,-1,-1,2) = 3δ(σ_i,σ_{i+1}) - I.
The transverse field X + X† = P + P† is also S₃-invariant (P → P† under transposition).

**The Z₃ clock Hamiltonian secretly has full S₃ symmetry!** This explains the exact BW match.

### Data:
- Ordered phase: S = log₂(3) = 1.585 bits (constant), locality 43-63%
- Peak: 76.5% at h/J = 0.75 (near critical region)
- BW envelope: linear wins over sin_inv throughout (unlike TFIM where sin_inv wins)
- Coupling profile: Unruh-like gradient (couplings decrease from far-from-cut to near-cut)

### Implication:
- Cannot test Z₃ vs S₃ with standard clock model — need chiral perturbation (φ≠0)
- The "prediction overturned" is actually "experiment didn't test what we thought"
- Redirect 34c: chiral clock model with genuine Z₃ (not S₃) symmetry

## Experiment 34c: Prediction Verification

*[Results to be filled]*
