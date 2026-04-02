# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 078 — S_q Potts g_c = 1/q exact (self-duality). DMRG walking: gap×N increases n=8→12 at q=5 (no breakdown). c_eff(q=5) = 1.15 matches complex CFT Re(c) ≈ 1.138 to 1%. S_q ≠ hybrid at q=2 (field factor 2).

## Active Research Thread
**S_q Potts complete characterization.** Self-duality gives exact g_c = 1/q. Walking regime extends beyond n=12. Complex CFT predictions (Gorbenko et al.) confirmed quantitatively.

S_q Potts data:

| q | g_c (exact) | c_eff (DMRG) | c_predicted | 1st degen | R_ε |
|---|------------|-------------|------------|-----------|-----|
| 3 | 1/3 | 0.89→0.80 | 4/5 (exact) | identical to hybrid | — |
| 5 | 1/5 | 1.15 | Re(c)≈1.138 | 4-fold (S₅) | ~3.9 |
| 7 | 1/7 | — | — | 6-fold (S₇) | 3.48 |
| 10 | 1/10 | — | — | (q-1)-fold | — |

Hybrid vs S_q: different g_c, different c, different degeneracies. Same model only at q=3.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **S_q Potts walking length at larger n** — DMRG at q=5 with higher chi (120-200) for n=24-50. Find where gap×N starts decreasing. Need ~10 min per size.
2. **Hardware validation** — 580s QPU remaining (52+ sprints since last use). Best candidate: conformal tower degeneracy (4-fold S₅ vs 2-fold Z₅) on n=6-8 qubits.
3. **S_q Potts c(q=7,10)** — DMRG entropy at exact g_c = 1/7 and 1/10. Compare to complex CFT predictions at higher q. Test if c_eff continues to match Re(c).

## What's Been Ruled Out
- g_c(S_q) ≠ 1/q: CONFIRMED g_c = 1/q exact (Sprint 078a)
- Walking breakdown at n≤12 for q=5 (Sprint 078b)
- S_q ≡ hybrid at q=2: WRONG, S_q field = X, hybrid = 2X (Sprint 078a)
- S_q = hybrid universality at q≥4 (Sprints 076-078)
- First-order S_q at q=5,7,10 at n≤8 (Sprints 076-077)
- All items from previous STATE.md

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤5. GPU: n≤10 for q≤5, n≤8 for q=7.
- S_q Potts DMRG: working for q=3,5 (TeNPy, orthogonal_to as kwarg)
- Exact g_c = 1/q for S_q Potts — no need for gap crossing scans
- IBM QPU: 580s remaining
