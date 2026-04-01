# Sprint 053 — ν(q) at True Critical Points: Method Validated, Old ν Values WRONG

**Date:** 2026-04-01
**Status:** Complete (6 experiments)

## Goal
Extract ν(q=3) at g_c=1/3 to validate methodology (exact ν=5/6), then apply to q=2,4,5,7,10.

## Key Results

**Validated method: corrected energy gap slope.**
d(Δ·N)/dg at g_c scales as N^{1/ν}·(1 + b/N), with b=0.86 from q=3 calibration.
Works to <1% for q=2 (ν=1.00) and ~3% for q=3 (ν=0.86 vs 5/6).

**ν(q) table — corrected:**

| q | g_c | Sizes | Raw ν | Corrected ν | Old ν (Sprints 045-047) | Exact (2D classical) |
|---|-----|-------|-------|-------------|------------------------|---------------------|
| 2 | 0.250 | 4,6,8 | 1.17/1.15/1.13 | 0.99/1.00/1.01 | 1.0 | 1.0 |
| 3 | 0.333 | 4,6,8,10 | 0.99→0.93 | 0.86 (1/N extrap) | 5/6 ✓ | 5/6 |
| 4 | 0.392 | 4,6,8 | 0.94/0.92/0.91 | 0.82/0.82/0.82 | ≥2.2 (WRONG) | 2/3 + logs |
| 5 | 0.441 | 4,6 | 0.97 | 0.85 | ~2.0 (WRONG) | — (1st order in 2D) |
| 7 | 0.535 | 4,6 | 1.14 | 0.97 | ~0.5 (WRONG) | — (1st order in 2D) |
| 10 | 0.684 | 4,5 | 1.37 | 1.12 | unreliable | — (1st order in 2D) |

**ALL old ν values for q≥4 were WRONG.** Sprint 045-047 used MI-CV data collapse which overestimates ν by ~40% at small sizes. The previous "BKT peak" at q=4 (ν≥2.2) and "large ν" at q=5 (ν≈2.0) are artifacts of bad methodology at wrong critical points.

**ν(q) for q=3,4,5 is nearly constant at ~0.82-0.86.** This contrasts dramatically with the old non-monotonic picture. The 2D classical exact values (1 → 5/6 → 2/3 for q=2,3,4) predict a monotone decrease; we see approximate constancy. q=7,10 values are from the smallest size pairs and may be overestimated.

**Data collapse is UNRELIABLE.** Gives ν=1.19 for q=3 (exact 5/6) and ν=1.16 excluding n=4. Systematically overestimates by ~40%.

## Experiments

### Exp A: Energy gap slopes q=3, n=4,6,8 (exact diag)
- 57 g-points per size near g_c=1/3. Fast: 0.6s (n=4), 2.1s (n=6), 6.9s (n=8).
- d(Δ·N)/dg at g_c: 11.25 (n=4), 16.97 (n=6), 22.96 (n=8)
- Raw pairwise ν: 0.986 (4,6), 0.951 (6,8) — converging downward
- Δ·N crossing points: 0.318, 0.322, 0.325 — converging to 1/3 from below

### Exp B: Push to n=10 (dim=59049, exact diag)
- 35 g-points in 39s. n=10 slope: 29.20
- Pairwise ν(8,10) = 0.928, continuing downward trend
- DMRG excited states FAILED (orthogonal_to gives gap=0 for q=3 Potts)

### Exp C: ν extrapolation and method comparison
- **Direct power-law** (slope ~ N^{1/ν}): ν = 0.961 (15% error)
- **Corrected power-law** (+ b/N term): ν = 0.859 (3.1% error), b = 0.86
- **1/N extrapolation** of pairwise ν: ν(∞) = 0.860 (3.2% error)
- **Data collapse** (fixed g_c): ν = 1.19 (43% error) — WORST method
- Crossing drift exponent: -1.87, consistent with ν=5/6, ω=4/5 prediction (-2.0)

### Exp D: ν(q=5) at g_c=0.441
- n=4,6 (dim 625, 15625): 1.7s, 12.7s
- Crossing at g=0.420 (FSS correction → g_c≈0.44)
- At g_c=0.441: raw ν=0.974, corrected ν=0.850

### Exp E: ν(q=2,4,7) at true g_c values
- **q=2 (sanity check):** 3 size pairs give corrected ν = 0.994, 1.000, 1.008. **Validates calibration.**
- **q=4:** 3 size pairs give corrected ν = 0.822, 0.823, 0.823. Remarkably consistent.
  Old ν≥2.2 → 0.82. NOT BKT. Not even anomalously large.
- **q=7:** n=4,6 only (dim 2401, 117649). At g_c=0.535: corrected ν=0.974.

### Exp F: ν(q=10) at g_c=0.684
- n=4 (10000), n=5 (100000) — pushed limits
- Crossing at g=0.646
- At g_c=0.684: corrected ν=1.12 (suspect — only (4,5) pair, b calibration may not transfer)

## Surprises
1. **ν(q=4) = 0.82, NOT ≥2.2** — the "BKT peak" was an artifact of MI-CV data collapse at wrong g_c
2. **ν(q=5) = 0.85, NOT 2.0** — same artifact
3. **Data collapse is the worst ν extraction method** at these sizes (43% error for q=3)
4. **ν(q=3,4,5) ≈ 0.82-0.86** — nearly constant, contrasting with old non-monotonic picture
5. **Corrected power-law and 1/N extrapolation agree** to 0.001 for q=3

## Methodology Assessment
| Method | Error (q=3) | Recommendation |
|--------|-------------|----------------|
| Data collapse (Δ·N) | 43% | **DO NOT USE** at n≤10 |
| Direct power-law fit | 15% | Rough estimate only |
| Pairwise slope ratio | 5-18% | Good with multiple sizes |
| Corrected power-law | 3% | **BEST** with ≥3 sizes |
| 1/N extrapolation | 3% | **BEST** with ≥3 consecutive pairs |
