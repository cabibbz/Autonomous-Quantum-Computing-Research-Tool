# Current State — Rewrite this completely each sprint

## Last Sprint
Sprint 088 — Tail weight scaling at q=2,3 (real CFT baseline). Power-law tail growth w_tail ~ n^2.0 is UNIVERSAL across all critical q=2,3,5 — not walking-specific. Sprint 087 exponents (1.72, 2.52) were real-space fitting artifacts; log-log gives b≈2.0 for all q. q=7 appears higher (b≈3.0) but only 4 data points. Entanglement gap closure Δξ ~ -0.43·ln(n) also universal.

## Active Research Thread
**Entanglement spectrum universality vs walking-specific effects.**

Tail weight scaling (DMRG midchain, open BC, log-log fit):

| q | exponent b | prefactor A | R² | %S(tail) at n=24 | CFT type |
|---|-----------|-------------|------|------------------|----------|
| 2 | 1.98 | 6.7e-6 | 0.97 | 4.4% | real (Ising) |
| 3 | 2.01 | 7.9e-6 | 0.97 | 4.2% | real (3-Potts) |
| 5 | 2.01 | 8.1e-6 | 0.97 | 3.8% | walking |
| 7 | 2.97 | 7.5e-7 | 0.96 | 1.4% (n=12) | broken |

Key reframing: walking breakdown manifests through level REDISTRIBUTION (how (q-1)-fold multiplet vs ground state share weight), NOT through differential tail growth. The tail grows the same for all q.

## QPU Budget
- Used: 20s of 600s (Sprint 025: ibm_kingston)
- Remaining: 580s

## Top 3 Next Experiments
1. **Hardware validation** — 580s QPU unused for 63 sprints. q=2 Ising at n=6-8 on QPU.
2. **q=7 tail weight at larger n** — Only 4 data points (n=6-12). Extend to n=16,20 via DMRG to test whether b≈3 is genuine or pre-asymptotic.
3. **Walking discriminator in level redistribution** — Compare %S(lev0) and %S(lev1) trajectories across q. The q-dependent part is HOW weight redistributes between ground+multiplet, not the tail.

## What's Been Ruled Out
- Tail weight growth as walking-specific (universal n^2 for all q, Sprint 088)
- Sprint 087 exponents 1.72, 2.52 (real-space fitting artifact)
- Tail weight saturation at accessible sizes (power law with R²>0.97)
- α=3 as Re(c) recovery probe (single-size extraction artifact, Sprint 086)
- c_∞ (min-entropy) as Re(c) recovery probe (degrades like c_1)
- Open-BC CC profile fit for Rényi c_α extraction (R²<0.9 at n≤12)
- Oscillatory correlator corrections at r ≤ 7 (Im(x_σ) too small)
- Open-BC raw power law for x extraction (inflates η by 5×)
- x_σ as walking discriminator (nearly constant across q)
- Casimir energy as walking discriminator (follows Re(c) for all q)
- Entanglement gap as walking discriminator (INCREASES with q at fixed n)
- q=6 as second walking case (c_eff drops 2.9% n=8→12)

## Key Tools Available
- 1D exact diag CPU: n≤8 for q≤6, GPU: n≤10 for q≤5, n≤8 for q=6
- S_q Potts DMRG: q=2,3 (very fast, n≤24+), q=5 (fast, n≤24), q=6 (slow, n≤12), q=7 (n≤12)
- Periodic exact diag correlator: q=2 n≤14, q=3 n≤10, q=5 n≤8, q=7 n≤7
- Exact g_c = 1/q for S_q Potts
- IBM QPU: 580s remaining
