# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 067 — Two-transition scan: hybrid has NO floating phase, unlike clock model.

## Active Research Thread
**Potts-clock hybrid universality class — FULLY characterized, phase structure confirmed.**

Four differences from clock model (q≥5):

| Property | Hybrid | Clock |
|----------|--------|-------|
| Transition type | Power-law 2nd order | BKT |
| ν | ~0.83 (finite) | →∞ |
| Floating phase | **None** (Sprint 067) | **Yes** (g=[0.30,0.92] at q=5) |
| c·x₁ | 0.112 | 0.146 |

Key CFT data for hybrid:

| q | g_c | c | ν | x₁ | c·x₁ | C_sse |
|---|-----|---|---|-----|------|-------|
| 2 | 0.250 | 0.500 | 1.00 | 0.125 | 0.062 | 0.50 |
| 3 | 0.333 | 0.800 | 0.86 | 0.133 | 0.107 | 0.54 |
| 5 | 0.441 | ~1.10 | 0.83 | ~0.10 | 0.112 | 0.38 |
| 10 | 0.684 | ~1.40 | 1.12? | ~0.083 | 0.116 | ~0.21 |

Sprint 067 findings:
- No second gap minimum at q=5,7 (scanned g=0.02-3.0)
- DMRG entropy: c_eff = 1.22 at g_c, drops to 0.003 at g=0.55 (sharp transition)
- Clock q=5 floating phase Δg=0.62 validated with same method
- Large-g gap*N → 2·sin(2π/q): 1.90 (q=5), 1.56 (q=7)

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — QPU budget 580s unspent. Strongest prediction: energy gap scaling at q=3 (exact ν=5/6) or q=2 conformal tower on real hardware. Measure Δ·N at multiple n to test FSS.
2. **2D generalization** — Does the hybrid universality class persist in 2D? Small plaquette exact diag (e.g., 2×3 at q=3) to check if 1D results extend.
3. **Entanglement Hamiltonian in floating phase** — The clock model's floating phase has c=1 (free boson). Measure BW locality in the floating phase vs at g_c — does the entanglement Hamiltonian become exactly local (BW=100%)?

## What's Been Ruled Out
- Floating phase in hybrid model (Sprint 067: 3 independent probes)
- Weakly first-order at q=10 (Sprint 066: no signal up to n=6)
- Hybrid = clock universality (Sprint 065+067: 4 differences)
- c·x₁ = 1/9 exact (Sprint 063c: q=3 exact is 8/75)
- c·x₁ universal to Z_q models (Sprint 063b: clock ≠ hybrid)
- Single compact boson for q≥5 (Sprint 062)

## Key Tools Available
- Exact diag CPU: n=4 for q≤30, n=6 for q≤10, n=8 for q≤5
- DMRG: q≤15 at n≤8 (chi=30), q≤5 at n≤24. Clock DMRG ~10x slower.
- IBM QPU: 580s remaining
- results.db: SQLite with all key measurements from sprints 50-67
