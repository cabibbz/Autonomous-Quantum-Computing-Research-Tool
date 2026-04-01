# Sprint 048 — ν(q=10) Extraction: Dead-Pair Bias Invalidates Raw MI-CV at Large d

**Status:** Complete (5 experiments)

## Motivation
ν(q) has a sharp peak at q=4 (≥2.2, BKT-like), then drops: q=5 (2.0), q=7 (0.5). Does ν(q=10) confirm the large-q → mean-field limit? Existing q=10 data (Sprint 043) used χ=10, which was known to be insufficient for d≥7. This sprint reruns with χ=20.

## Key Questions
1. Does ν(q=10) ≈ 0.5 (mean-field), confirming trend from q=7?
2. Are MI-CV crossings present at q=10 with χ=20?
3. Does the crossing g_c match the scaling law prediction g_c(10)=0.246?

## Experiments

### 048a: q=10 MI-CV at n=8, χ=20
Dense sweep over g ∈ [0.15, 0.35] (9 points). Total time: 810s.

CV is remarkably flat across the entire range:
| g | CV | E₀ |
|---|---|---|
| 0.150 | 0.597 | -7.228 |
| 0.200 | 0.601 | -7.409 |
| 0.240 | 0.611 | -7.597 |
| 0.250 | 0.613 | -7.650 |
| 0.260 | 0.614 | -7.706 |
| 0.280 | 0.616 | -7.802 |
| 0.300 | 0.617 | -7.956 |
| 0.350 | 0.611 | -8.341 |

CV varies only 3% (0.597–0.617) across the entire transition. Compare q=7 n=8 where CV varies from 0.597 to 0.438 (26% range). The d=10 local Hilbert space smears out the CV signal.

### 048b: q=10 MI-CV at n=12, χ=20
Focused sweep (8 points). Total time: 1094s.

**NO CROSSINGS.** n=12 CV is BELOW n=8 CV at ALL g values:
| g | n=8 CV | n=12 CV | diff |
|---|---|---|---|
| 0.150 | 0.597 | 0.422 | -0.175 |
| 0.200 | 0.601 | 0.481 | -0.120 |
| 0.240 | 0.611 | 0.503 | -0.108 |
| 0.250 | 0.613 | 0.508 | -0.105 |
| 0.260 | 0.614 | 0.512 | -0.102 |
| 0.280 | 0.616 | 0.533 | -0.084 |
| 0.300 | 0.617 | 0.537 | -0.080 |
| 0.350 | 0.611 | 0.508 | -0.103 |

Gap narrows near g_c≈0.25 (minimum -0.080 at g=0.30) then widens at g=0.35.

**Sprint 043 χ=10 "crossing confirmed at g_c≈0.246" is INVALIDATED.** At χ=10, n=12 in the ordered phase (g=0.1) had mean MI ≈ 3×10⁻⁵ (essentially zero) → DMRG failed to break symmetry → artificial high CV. At χ=20, n=12 g=0.15 has mean MI = 0.42 (well-converged).

### 048c: ν extraction from raw MI-CV
Slope ratio dCV/dg:
- At [0.20, 0.26]: slope ratio = 2.46, ν = 0.45
- At [0.24, 0.30]: slope ratio = 5.88, ν = 0.23
- Global [0.20, 0.30]: ratio = 3.61, ν = 0.32

Data collapse: ν = 0.24, g_c = 0.339 (far from known g_c ≈ 0.246). Collapse is overfitting with only 2 sizes. ν extraction from raw MI-CV is unreliable.

### 048d: χ convergence check — ENERGY CONVERGES BUT MI DOES NOT
At n=8, g=0.20:
| χ | CV | E₀ | mean MI | time |
|---|---|---|---|---|
| 20 | 0.601 | -7.409391 | 1.264 | 128s |
| 30 | 0.588 | -7.409415 | 1.512 | 191s |
| 40 | 0.370 | -7.409442 | 1.824 | 263s |

**CRITICAL: Energy converges (5th decimal) but MI does NOT.** Mean MI jumps 44% from χ=20→40. The "dead pairs" at χ=20 are partially artifacts — they acquire non-zero MI at higher χ. CV drops 38% from χ=30→40. At d=10, χ=40 captures only 40/10000 possible Schmidt values.

**Energy convergence ≠ MI convergence for large d.** This is because MI requires accurate local 2-site density matrices, not just the global energy. The bond dimension needed for MI convergence is much larger than for energy convergence when d is large.

### 048e: Filtered MI-CV — Dead-Pair Bias Mechanism
**CRITICAL METHODOLOGICAL FINDING.**

At n=8 (28 total pairs), 7 pairs (25%) have MI < 0.01. At n=12 (66 total pairs), 11 pairs (17%) have MI < 0.01. This asymmetric dead-pair fraction creates a systematic bias in CV.

When near-zero MI pairs are filtered out (threshold = 0.01):
| g | n=8 filt | n=12 filt | diff |
|---|---|---|---|
| 0.150 | 0.234 | 0.249 | +0.015 |
| 0.200 | 0.242 | 0.273 | +0.031 |
| 0.250 | 0.264 | 0.313 | +0.049 |
| 0.260 | 0.266 | 0.320 | +0.053 |
| 0.280 | 0.270 | 0.294 | +0.023 |
| 0.300 | 0.271 | 0.300 | +0.029 |
| 0.350 | 0.262 | 0.285 | +0.024 |

**Sign flips!** n=12 > n=8 everywhere when dead pairs removed. But still no crossings — n=12 > n=8 at ALL g. This is the same q=4 pattern (monotonically ordered n=12 > n=8 everywhere). Slope ratio < 1 in filtered data → ν → ∞ (BKT-like).

## Key Findings

1. **q=10 raw MI-CV shows NO crossings at χ=20** — partially a dead-pair artifact AND partially χ non-convergence
2. **Sprint 043's χ=10 crossing was wrong**: insufficient χ for d=10
3. **Dead-pair bias mechanism**: fraction of near-zero MI pairs decreases with n at large d, creating systematic CV deflation at larger n. BUT dead pairs are themselves partially χ artifacts at d=10.
4. **Energy convergence ≠ MI convergence**: E₀ converges at χ=20 (5th decimal) but MI structure changes dramatically through χ=40. At d=10, MI-based analysis requires much higher χ than energy-based analysis.
5. **MI-CV as FSS tool breaks down for d ≥ 10**: neither raw nor filtered MI-CV is reliable at χ=20-40
6. **ν extraction not possible** at q=10 with current DMRG capability

## Surprises
- Raw MI-CV gives OPPOSITE size ordering from filtered truth (n=12 < n=8 raw vs n=12 > n=8 filtered)
- n=8 CV is remarkably flat (3% variation) — d=10 smears the transition signal
- **χ=20 → χ=40 changes CV by 38%** while energy changes by 0.0005% — MI and energy converge at vastly different rates
- Mean MI jumps 44% from χ=20→40 — "dead pairs" are partially χ artifacts
- DMRG energy convergence is NOT sufficient to trust MI-derived quantities at large d

## Implications
1. **MI-CV FSS requires χ > d² for reliable results.** At d=10, need χ > 100 (infeasible for current setup).
2. **Alternative FSS probes needed for large d:** half-chain entropy (S vs g peak), correlation length from transfer matrix, or entanglement spectrum gaps. These depend on fewer density matrix elements than all-pairs MI.
3. **Previous q=10 and q=20 results (Sprint 043) should be treated as qualitative only.**
4. **q=2,3 MI-CV results (d=2,3) are reliable** — χ=20 is well above d²=4,9.
5. **q=4,5,7 results (d=4,5,7) need scrutiny** — χ=20 gives d²=16,25,49, marginal to insufficient.
