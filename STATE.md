# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 066 — Weakly first-order test at q=10: NO first-order signal. Transition is continuous.

## Active Research Thread
**Potts-clock hybrid is a genuinely new universality class — now FULLY characterized.**

Three models, three universality classes for q≥5:

| Model | q>4 transition | ν | Key evidence |
|-------|---------------|---|-------------|
| S_q Potts | First-order | — | Literature (Gorbenko et al.) |
| Z_q clock | BKT (ν→∞) | ∞ | Sprint 065c: slopes flat, ν diverges |
| **Our hybrid** | **Continuous 2nd order** | **~0.83** | Sprint 065c + 066: no first-order signal |

Key CFT data for hybrid:

| q | g_c | c | ν | x₁ | c·x₁ | C_sse | BW peak |
|---|-----|---|---|-----|------|-------|---------|
| 2 | 0.250 | 0.500 | 1.00 | 0.125 | 0.062 | 0.50 | 91% |
| 3 | 0.333 | 0.800 | 0.86 | 0.133 | 0.107 | 0.54 | 76.5% |
| 5 | 0.441 | ~1.10 | 0.83 | ~0.10 | 0.112 | 0.38 | 42% |
| 10 | 0.684 | ~1.40 | 1.12? | ~0.083 | 0.116 | ~0.21 | 39% |

Sprint 066 findings:
- Δ·N at q=10 increases but DECELERATES (continuous, not first-order)
- Gap minimum ~0.030 size-independent (no exponential closing)
- q=3 converges from above, q≥5 from below (sign-flipped corrections)
- Anomaly saturates ~5.5% for q≥7

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Two-transition scan** — Does hybrid have a SECOND transition like clock? Scan g over wide range at q=5,7 looking for a second critical point or intermediate floating phase.
2. **Hardware validation** — QPU budget 580s unspent. Test gap scaling at q=3 (known ν=5/6) or conformal tower structure on real hardware.
3. **Entanglement entropy on hardware** — Can we measure S(A) scaling on QPU via randomized measurements? Test c=0.5 (Ising) as proof of concept.

## What's Been Ruled Out
- Weakly first-order at q=10 (Sprint 066: no signal up to n=6)
- Hybrid = clock universality (Sprint 065: ν qualitatively different)
- c·x₁ = 1/9 exact (Sprint 063c: q=3 exact is 8/75)
- c·x₁ universal to Z_q models (Sprint 063b: clock ≠ hybrid)
- Single compact boson for q≥5 (Sprint 062)
- Luttinger liquid for q≥4 (Sprint 062)

## Key Tools Available
- Exact diag CPU: n=4 for q≤30, n=6 for q≤10, n=8 for q≤5
- DMRG: q≤15 at n≤8 (chi=30), q≤5 at n≤24. Clock DMRG ~10x slower.
- IBM QPU: 580s remaining
- results.db: SQLite with all key measurements from sprints 50-66
