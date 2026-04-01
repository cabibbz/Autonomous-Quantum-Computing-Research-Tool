# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 061 — Large-q Limit: CFT Data at q=15-30 and Free Boson Convergence

## Active Research Thread
**Complete characterization of 1D quantum Potts CFT, now extended to q=30.** Eight quantities measured:

| q | g_c | c | ν | x₁ | c/x₁ | C_sse | R(σ²)/4 |
|---|-----|---|---|-----|------|-------|---------|
| 2 | 0.250 | 0.500 | 1.00 | 0.125 | 4.0 | 0.50 | — |
| 3 | 0.333 | 0.800 | 0.86 | 0.133 | 6.0 | 0.54 | 0.25 |
| 5 | 0.441 | ~1.10 | 0.85 | ~0.10 | 10.8 | 0.38 | 0.60 |
| 10 | 0.684 | ~1.40 | 1.12 | ~0.083 | 16.9 | ~0.21 | 0.92 |
| 15 | 0.798* | ~1.61* | — | ~0.071† | — | 0.09† | 0.98 |
| 20 | 0.922* | ~1.73* | — | ~0.055† | — | 0.06† | 0.99 |
| 30 | 1.127* | ~1.90* | — | ~0.038† | — | 0.04† | 1.00 |

*formula prediction, †n=4 only

**Sprint 061 key findings:** Free boson is the q→∞ limit. Harmonic ratios converge to k² as O(q^{-3}). R_epsilon grows linearly with q at n=4. C_sse ~ q^{-0.9}. x₁ → 0 (not saturating).

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **QPU validation of Potts CFT** — 580s budget unspent. Encode q=3 Potts in qutrit-from-qubit mapping, measure energy gap ratio on hardware. First hardware test of Potts critical physics.
2. **R_epsilon at larger n for q=15** — Is R_eps=13.7 at n=4 a finite-size artifact or physical? DMRG at g_c=0.798 for n=16-32, extract entropy profile for c(q=15).
3. **Free boson radius from scaling dimension** — If x₁ = 1/(4πR²), extract effective radius R(q) from measured x₁. Does R ~ √(ln q) as predicted?

## What's Been Ruled Out
- c/x₁ = 2q for q≥4 (Sprint 058)
- Free boson at finite q (Sprint 057, but approaches rapidly as q→∞)
- x₁ saturation at large q (Sprint 061: x₁ continues to 0.038 at q=30)
- Epsilon merging with sigma² (Sprint 061: different charge sectors, diverge)
- n=4,5 energy gap crossing for g_c at large q (21% FSS, unreliable)
- Simple analytic formula for C_sse(q) (Sprint 060-061)

## Key Tools Available
- Exact diag: n=4 for q≤30 (810k dim), n=6 for q≤10 (1M dim)
- DMRG: open BC, any q at moderate n (chi=20-80)
- Z_q charge resolution via G = X₁⊗...⊗Xₙ
- IBM QPU: 580s remaining
