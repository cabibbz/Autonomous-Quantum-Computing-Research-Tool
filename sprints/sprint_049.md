# Sprint 049 — Correlation Length & Entropy FSS: q=3 Potts Critical Point WRONG

**Status:** Complete (6 experiments)

## Motivation

MI-CV has systematic problems at large d: dead-pair bias (Sprint 048), χ non-convergence of all-pairs MI, and Gell-Mann failure for d≥4. This sprint develops MI-independent methods for ν(q) extraction: correlation length from correlator decay, and half-chain entropy finite-size scaling.

## Key Questions
1. Can we extract ν=1 for TFIM from correlation length and entropy FSS?
2. Can we measure the central charge c directly from S(g_c, n) ~ (c/6)ln(n)?
3. Do these methods work for q≥3 Potts?

## Experiments

### 049a: TFIM Correlation Length Validation (χ=60, n=40,80)

Extracted ξ from exponential decay of ⟨σx_i σx_j⟩_conn using bulk sites (n/4 to 3n/4).

| g | ξ(n=40) | ξ(n=80) | ξ_exact=1/ln(g) | R²(80) |
|---|---------|---------|-----------------|--------|
| 1.02 | 13.4 | 20.9 | 50.5 | 0.967 |
| 1.05 | 9.7 | 12.9 | 20.5 | 0.983 |
| 1.10 | 6.7 | 8.0 | 10.5 | 0.992 |
| 1.20 | 4.2 | 4.7 | 5.5 | 0.997 |
| 1.50 | 2.1 | 2.3 | 2.5 | 0.999 |
| 2.00 | 1.3 | 1.4 | 1.4 | 1.000 |

**Findings:**
1. Method works: R² > 0.95 in disordered phase, clean exponential decay
2. Finite-size saturation: ξ < ξ_exact when ξ_exact > n/4
3. ξ(80)/ξ(40) = 1.93 ≈ 2.0 at g_c → confirms criticality
4. Naive ν extraction: 0.62 (n=40), 0.72 (n=80). Both far below ν=1 due to finite-size saturation. The method underestimates ν because ξ is capped at ~n/4.

### 049b: TFIM + Potts Correlation Length (q=2,3)

q=2 at n=20,40,80 with χ=60: same ξ trends as 049a. Consistent results.

**q=3 Potts correlation length: CATASTROPHICALLY WRONG.** ξ ≈ 0.7 at ALL tested g values (0.80-1.02). Entropy S = 0.047 constant across all sizes. Initially suspected DMRG convergence failure, but exact diag at n=8 CONFIRMS S=0.047 and E=-18.608.

The q=3 Potts model at g=1.0 is deep in the disordered phase with ξ < 1. The critical point is NOT at g=1.0.

### 049c: TFIM Entropy FSS (n=16,24,32,48,64)

S(n, g) at g_c=1.0:

| n | S(g_c=1.0) |
|---|-----------|
| 16 | 0.4234 |
| 24 | 0.4602 |
| 32 | 0.4857 |
| 48 | 0.5211 |
| 64 | 0.5458 |

**Central charge c = 0.529 from full fit, converging to exact 0.500:**
- Pairwise c: 0.544 (16,24) → 0.532 (24,32) → 0.523 (32,48) → 0.516 (48,64)
- Clean convergence toward c=1/2 from above

**Pseudo-critical point from dS/dg peak shifts toward g=1.0:**
- n=16: g_pseudo = 0.930
- n=32: g_pseudo = 0.970
- n=64: g_pseudo = 0.990

**Entropy FSS collapse: ν extraction fails.** Best collapse at ν=3.0, not ν=1.0. This is expected — entropy has logarithmic singularity S ~ (c/6)ln(ξ) at criticality, not a power law. Standard FSS collapse form doesn't apply to entropy directly.

### 049d-e: q=3 Potts Wide Scan — True Critical Point Found

**Exact diag (n=6,8) reveals the full entropy landscape:**

| g | S(n=6) | S(n=8) | Phase |
|---|--------|--------|-------|
| 0.00 | 0.000 | 0.000 | Ground state degenerate |
| 0.02 | 1.099 | 1.099 | Ordered (GHZ, S=ln3) |
| 0.15 | 1.092 | 1.099 | Ordered |
| 0.20 | 1.042 | 1.089 | Ordered → critical |
| 0.25 | 0.864 | 0.988 | **S grows with n** |
| 0.30 | 0.607 | 0.698 | **S grows with n** |
| 0.35 | 0.416 | 0.487 | Transitioning |
| 0.50 | 0.183 | 0.184 | Disordered (S saturated) |
| 1.00 | 0.047 | 0.047 | Deep disordered |

**DMRG n=12 confirms (chi=40):**

| g | S(n=8) | S(n=12) | ΔS | c estimate |
|---|--------|---------|-----|-----------|
| 0.24 | 1.024 | 1.099 | +0.075 | 1.10 |
| 0.28 | 0.826 | 0.999 | +0.173 | 2.56 |
| 0.30 | 0.698 | 0.851 | +0.153 | 2.26 |
| 0.32 | 0.582 | 0.678 | +0.096 | 1.42 |
| 0.36 | 0.413 | 0.438 | +0.025 | 0.37 |

Central charge c crosses the exact q=3 Potts value c=4/5=0.800 between g=0.32 (c=1.42) and g=0.36 (c=0.37), suggesting **g_c ≈ 0.33-0.35 for q=3 Potts** with our Hamiltonian normalization.

### Why g_c ≠ 1.0

TFIChain uses Z₂ parity conservation, which the Potts Hamiltonian H = -J Σ_a P_a⊗P_a - g(X+X†) does NOT have (X+X† breaks any Z_q charge). The self-dual point depends on the normalization of the kinetic term relative to the coupling.

For the standard Potts model on a square lattice, self-duality gives exp(K_c) = 1+√q, i.e., K_c = ln(1+√q). The quantum-classical mapping relates J/g to the anisotropy of the 2D lattice. At the isotropic point:

g_c/J = 1/(2(q-1)) × ... (model-specific)

For q=2: the kinetic term is σz (eigenvalues ±1) and coupling is σx⊗σx. The self-dual condition gives g_c = J = 1.0. ✓

For q=3: the kinetic term X+X† has eigenvalues {2, -1, -1} (range = 3) while the coupling has range 1 (same/different). The asymmetry shifts g_c to ~1/3.

**Key insight: Sprints 038-048 were studying a CROSSOVER in the disordered phase, not the actual phase transition.** The MI-CV crossings near g≈0.9-1.0 are real phenomena but occur deep in the gapped disordered phase of the q=3 Potts model.

---

## Summary of Findings

1. **TFIM central charge c ≈ 0.52, converging to exact 0.50.** Clean logarithmic scaling S = (c/6)ln(n) + const at g_c=1.0 across 5 sizes (n=16-64). Pairwise c converges: 0.544 → 0.516.

2. **Correlation length method works but is finite-size limited.** ξ from correlator decay agrees with exact ξ=1/ln(g) far from g_c (within 5% at g=2.0) but saturates at ξ ~ n/4 near g_c. Naive ν extraction gives 0.62-0.72, well below ν=1.

3. **Entropy FSS is NOT suitable for ν extraction.** The logarithmic singularity of entropy at criticality (S ~ ln ξ ~ ν ln|g-g_c|) means standard power-law collapse fails. Central charge extraction works; ν does not.

4. **q=3 Potts critical point is at g_c ≈ 0.33, NOT 1.0.** Exact diag confirmed. The entropy landscape shows ordered (GHZ-like, S=ln3) for g<0.25, critical transition at g≈0.33, and disordered (product-like, S→0) for g>0.40. All prior q≥3 Potts results (Sprints 038-048) assumed g_c≈1.0 and were studying a disordered-phase crossover.

## Surprises

- **q=3 Potts at g=1.0 has S=0.047** — confirmed by exact diag. Not a DMRG convergence failure.
- **GHZ entropy S=ln(3)** persists deep into the ordered phase (up to g≈0.20 for n=8)
- **Entropy FSS collapse gives ν=3 for TFIM** — the method doesn't work for ν (logarithmic vs power-law singularity)
- **10 sprints of Potts results were at the wrong g** — the most expensive mistake of the project so far
- **TFIChain Z₂ conservation is essential** — without it, q=2 would also fail at large n

## Impact on Prior Results

| Sprint | Claim | Status |
|--------|-------|--------|
| 038 | q=3 Potts MI-CV crossings at g≈0.9 | Crossings real, but in disordered phase |
| 039 | Potts data collapse, ν=5/6 | Measured at wrong g_c |
| 040-042 | q=4,5 Potts MI-CV | All at wrong g values |
| 043 | q=10,20 "second-order for all q" | Phase location uncertain |
| 044 | g_c scaling law | Fitted to wrong g_c values |
| 045-048 | ν(q) extraction | All at wrong critical points |

The qualitative MI-CV signatures (crossings, domes, steps) remain valid as phenomena in the disordered phase. Whether they correlate with the ACTUAL phase transition behavior is an open question.

## Next Steps

1. **Find true g_c for q=3-10 Potts** — use entropy peak / central charge method validated here
2. **Re-examine MI-CV at true g_c** — do the crossing signatures persist near the real critical point?
3. **Implement iDMRG** for correlation length — avoids finite-size saturation problem
4. **Investigate the disordered-phase crossover** — what causes MI-CV crossings at g>>g_c?
