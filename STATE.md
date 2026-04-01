# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 063 — Testing c·x₁: c(q=15) confirmed, clock model differs, not exactly 1/9

## Active Research Thread
**1D quantum Potts CFT characterization complete through q=30.** Key data table:

| q | g_c | c | ν | x₁ | c·x₁ | C_sse |
|---|-----|---|---|-----|------|-------|
| 2 | 0.250 | 0.500 | 1.00 | 0.125 | 0.062 | 0.50 |
| 3 | 0.333 | 0.800 | 0.86 | 0.133 | 0.107 | 0.54 |
| 5 | 0.441 | ~1.10 | 0.85 | ~0.10 | 0.112 | 0.38 |
| 10 | 0.684 | ~1.40 | 1.12 | ~0.083 | 0.116 | ~0.21 |
| 15 | 0.798* | 1.549 | — | ~0.071† | 0.110 | 0.09† |

*formula, †n=4 only

**Sprint 063 findings:** c·x₁ ≈ 0.112 is Potts-specific (clock gives 0.146), approximate (not 1/9 exactly), but robust (holds to 1% at q=15). Clock g_c corrected.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **QPU validation** — 580s unspent. Strongest prediction: q=3 Potts energy gap ratio at g_c=1/3, or Ising scaling dimension x₁=1/8 on hardware. Worth burning budget.
2. **c·x₁ mechanism** — WHY is c·x₁ ≈ constant for Potts? The (q-1) spin primaries have total "weight" proportional to c. Each has dimension ~x₁. If weight per primary ≈ constant, c·x₁ ~ (q-1)·const·const = const. Test degeneracy-weighted sum rule.
3. **Entanglement Hamiltonian at q>4** — H_E = -log(ρ_A) structure at the novel q>4 critical points. Does BW locality hold? H/G-inv ratio predictions from Sprint 034.

## What's Been Ruled Out
- c·x₁ = 1/9 exact (Sprint 063c: q=3 exact is 8/75)
- c·x₁ universal to Z_q models (Sprint 063b: clock ≠ Potts)
- Single compact boson for q≥5 (Sprint 062)
- Luttinger liquid for q≥4 Potts (Sprint 062)
- Old clock g_c ≈ 0.67 (Sprint 063b: true value 0.52)

## Key Tools Available
- Exact diag: n=4 for q≤30, n=6 for q≤10, n=8 for q≤5
- DMRG: q≤15 at n≤8 (chi=30), q≤5 at n≤24
- IBM QPU: 580s remaining
