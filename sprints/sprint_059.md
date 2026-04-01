# Sprint 059 — Conformal Tower Analysis: Descendant Structure in q>4 Potts CFT

**Status:** In progress

## Motivation

We have measured five CFT quantities for q=2-10 Potts: g_c, c, ν, x₁, and operator content (Sprints 051-058). The key open question: is the q>4 Potts critical point a **genuine CFT** with Virasoro descendant structure, or only approximately conformal?

In a true CFT on a periodic chain of length N, the energy spectrum organizes into **conformal towers**:
- Each primary operator with scaling dimension x has descendants at x+1, x+2, ...
- Descendants from L₋₁ (or L̄₋₁) carry momentum k=±1 relative to the primary
- The descendant degeneracy at level n follows the Virasoro partition function

**CFT tower prediction for primary at x₁ (spin field σ):**
- Level 0: x₁ at k=0, degeneracy = (q-1) [spin field primaries]
- Level 1: x₁+1 at k=±1, degeneracy = 2(q-1) [L₋₁σ and L̄₋₁σ]
- Level 2: x₁+2 at k=0,±2, degeneracy = more complex

In gap ratios R = (E-E₀)/Δ₁: descendants at R = 1 + n/x₁.

**Literature context:**
- Lao et al. (PRB 2019, arXiv:1811.11189): Approximate towers for q=5,6 but conformal data drift — they assumed weak first-order
- Tang et al. (PRL 2024, arXiv:2403.00852): Clean towers for q=5 but only with non-Hermitian model
- **No clean Hermitian q>4 tower with momentum resolution exists in the literature**
- Our chain shows genuinely second-order transitions (Sprint 043), so towers should converge

## Plan
1. **Exp 059a:** Validate momentum-resolved tower on q=2 Ising (known: σ at x=1/8, L₋₁σ at x=9/8)
2. **Exp 059b:** q=3 Potts towers (known: σ at x=2/15, descendants at 2/15+1)
3. **Exp 059c:** q=4,5 towers — first momentum-resolved measurement in novel regime
4. **Exp 059d:** q=7,10 towers — deep novel regime

## Experiments

