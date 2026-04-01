# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 057 — Operator Content & Scaling Dimensions: What CFT Describes q>4 Potts?

## Active Research Thread
**CFT operator content of 1D quantum Potts.** Energy spectrum on periodic chains reveals the full operator content. The Z_q spin field sector has exactly (q-1) primaries organized in conjugate pairs. Harmonic ratios x(σᵏ)/x(σ) are sub-quadratic and approach k² (free boson) as q→∞.

| q | c | x₁ | Spin fields | R_ε | R(σ²/σ) |
|---|---|-----|------------|-----|---------|
| 2 | 0.500 | 0.125 | 1 | 8.0 | — |
| 3 | 0.800 | 0.133 | 2 | 6.0 | — |
| 4 | 1.000 | ~0.12 | 3 | 6.6 | 1.68 |
| 5 | ~1.10 | ~0.10 | 4 | 7.1 | 2.41 |
| 7 | ~1.30 | ~0.09 | 6 | 8.1 | 3.32 |
| 10 | ~1.40 | — | 9 | 8.3 | 3.66 |

**Key finding:** The growing number of spin field primaries (q-1) below the energy scale directly explains c(q) ~ ln(q). The q→∞ limit approaches a free boson.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Larger sizes for q=4,5 absolute x₁** — q=4 at n=10 (dim=10⁶, need sparse or DMRG excited states) or use Lanczos at n=8 with more eigenvalues to identify conformal towers more precisely.
2. **x₁(q) formula** — With x₁ reliably extracted for q=2,3, see if x₁ ∝ 1/q or follows a parafermion pattern. Need larger sizes for q≥4 (currently too small).
3. **QPU validation of g_c(q=3)=1/3** — 580s budget. Prepare ground state approximation near g_c, measure order parameter.

## What's Been Ruled Out
- Free boson k² for harmonic ratios at any finite q (Sprint 057)
- sin²(πk/q)/sin²(π/q) formula — 2% for q=7 k=2 but 15-19% elsewhere (Sprint 057)
- iDMRG S-vs-ln(xi) for c extraction at criticality (Sprint 055)
- Small-scale QEC active correction (Sprints 026-028)
- Pauli fraction as BW metric (Sprint 035)
- Gell-Mann MI for q≥4 (Sprints 045-046)
- Entropy FSS for ν extraction (Sprint 049)
- All Sprints 038-048 Potts g_c values: WRONG (Sprint 049)
- Data collapse for ν at n≤10: UNRELIABLE (Sprint 053)
- q=7+ DMRG profile method: too slow at d≥7 (Sprint 055)
- Analytic continuation of Potts CFT for q>4: WRONG (Sprint 056)

## Key Tools Available
- Exact diag with periodic BC: n≤12 (q=2), n≤10 (q=3), n≤8 (q=4), n≤6 (q=5,7,10)
- DMRG (TeNPy): finite chain, open BC only
- Energy gap Δ·N crossing: best method for g_c
- CFT spectrum on periodic chain: scaling dimension ratios from gap ratios
- IBM QPU: 580s remaining
