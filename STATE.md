# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 071 — Ly=2 cylinder geometry for q=2,5. g_c(cyl) measured: q=2 = 0.451, q=5 = 0.714. No first-order signal for q=5 (smooth order parameter, convergent gap crossing). DMRG impractical for q≥5 cylinder; exact diag gap×Lx crossing is reliable.

## Active Research Thread
**2D extension of Potts-clock hybrid — still inconclusive for full 2D at q=5.**

Cylinder (Ly=2 ladder) geometry results:

| q | g_c(1D) | g_c(cyl) | g_c(2D) | cyl/1D | Lx pairs |
|---|---------|----------|---------|--------|----------|
| 2 | 0.250 | 0.451 | 0.771 | 1.80 | (4,5),(5,6),(6,7) |
| 5 | 0.441 | 0.714 | 1.588 | 1.62 | (3,4) |

Full 2D (torus) results from Sprint 068:

| q | g_c(2D) | Sizes | gap×L | Verdict |
|---|---------|-------|-------|---------|
| 2 | 0.771 | L=2,3,4 | 2.363 | Continuous confirmed |
| 3 | 1.267 | L=2,3 | 6.10 | Consistent with continuous |
| 5 | 1.588 | L=2,3 | 3.50 | Inconclusive (2 sizes only) |

Sprint 071 cylinder findings:
- Order parameter ⟨δ(s_i,s_j)⟩ smooth for q=5 on cylinder — no first-order signal
- Gap×Lx crossing exists and converges for q=2 (3 pairs)
- q=5 has single crossing pair (Lx=3,4) — needs more sizes
- DMRG impractical for q≥5 cylinder (d=5 per site too large for modest chi)

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **q=3 cylinder comparison** — q=3 is exactly solvable (like q=2). Measure g_c(cyl, q=3) from gap crossings Lx=3-6 (dim up to 3^12=531k). Compare cyl/1D ratio to q=2,5. Continuous behavior would further support universality of the hybrid's continuous transition.
2. **Hardware validation** — 580s QPU remaining. Best candidate: 1D q=2 conformal tower (n=6-8 qubits). Clearest numerical prediction to test.
3. **Cylinder Ly=3** — For q=2 only (dim=2^{3·Lx}=8^Lx). Lx=4: dim=4096, Lx=5: dim=32k. Compare g_c convergence as Ly→∞ toward 2D g_c=0.771.

## What's Been Ruled Out
- Floating phase in 1D hybrid model (Sprint 067)
- Weakly first-order in 1D at q=10 (Sprint 066)
- Hybrid = clock universality (Sprint 065+067)
- c·x₁ = 1/9 exact (Sprint 063c)
- L=2 in scaling regime for 2D torus (Sprints 068-070)
- DMRG for q≥5 cylinders at chi<50 (Sprint 071b — too slow/inaccurate)
- Entropy peak as g_c proxy on finite-width cylinders (Sprint 071a — drifts with Lx)
- First-order transition on Ly=2 cylinder for q=5 (Sprint 071c — smooth order parameter)

## Key Tools Available
- 2D torus exact diag: L×L up to ~10 sites total
- Cylinder exact diag: Ly=2, Lx up to 4-5 for q=5 (dim≤10M)
- 1D exact diag CPU: n≤8 for q≤5. GPU extends to n≤10.
- DMRG: q≤5 1D chains at n≤24. q=2 cylinder at Lx≤20. q≥5 cylinder IMPRACTICAL.
- IBM QPU: 580s remaining
