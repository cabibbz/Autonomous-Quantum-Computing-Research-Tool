# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 062 — Compactification Radius & Twisted BC: c·x₁ ≈ 1/9 Discovery

## Active Research Thread
**Complete characterization of 1D quantum Potts CFT, now with c·x₁ constraint.** Nine quantities measured:

| q | g_c | c | ν | x₁ | c/x₁ | c·x₁ | C_sse | R(σ²)/4 |
|---|-----|---|---|-----|------|------|-------|---------|
| 2 | 0.250 | 0.500 | 1.00 | 0.125 | 4.0 | 0.062 | 0.50 | — |
| 3 | 0.333 | 0.800 | 0.86 | 0.133 | 6.0 | 0.107 | 0.54 | 0.25 |
| 5 | 0.441 | ~1.10 | 0.85 | ~0.10 | 10.8 | 0.112 | 0.38 | 0.60 |
| 10 | 0.684 | ~1.40 | 1.12 | ~0.083 | 16.9 | 0.116 | ~0.21 | 0.92 |
| 30 | 1.127* | ~1.90* | — | ~0.038† | — | — | 0.04† | 1.00 |

*formula prediction, †n=4 only

**Sprint 062 key findings:** c·x₁ ≈ 1/9 for q≥3 (POTENTIALLY NOVEL). Single compact boson ruled out. q=2,3 are Luttinger liquids, q≥4 are not. Twisted BC spin stiffness confirms breakdown.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **c·x₁ = 1/9 test at larger q** — Measure c and x₁ independently at q=15,20 to check if c·x₁ = 1/9 holds. Requires entropy profile at larger n for c, or use c from log formula + x₁ from desc gap.
2. **QPU validation of Potts CFT** — 580s budget unspent. Encode q=3 Potts in qutrit-from-qubit mapping, measure energy gap ratio on hardware.
3. **Effective number of bosons** — If c·x₁=const, the system has n_eff = c/1 = c effective modes each with x₁_eff = 1/(9n_eff) = 1/(9c). Test: does spectrum factorize into c independent sectors?

## What's Been Ruled Out
- Single compact boson for q≥5 (Sprint 062: c>1)
- Luttinger liquid for q≥4 (Sprint 062: ρ_s·L drifts 10-22%)
- c/x₁ = 2q for q≥4 (Sprint 058)
- Free boson at finite q (Sprint 057, approaches as q→∞)
- x₁ saturation at large q (Sprint 061)
- Simple analytic formula for C_sse(q) or c/x₁(q)

## Key Tools Available
- Exact diag: n=4 for q≤30 (810k dim), n=6 for q≤10 (1M dim)
- DMRG: open BC, any q at moderate n (chi=20-80)
- Z_q charge resolution via G = X₁⊗...⊗Xₙ
- IBM QPU: 580s remaining
