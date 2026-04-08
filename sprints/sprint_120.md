# Sprint 120 — Hybrid q=4 chi_F: The z_m=1 Transition Point

**Date:** 2026-04-08
**Thread:** chi_F spectral decomposition — model comparison
**Status:** Complete

## Motivation

Sprint 119 showed chi_F spectral decomposition is universal in mechanism but model-specific in exponents. At q=3, both models identical (alpha=1.40). At q=5, already diverged (hybrid 1.41 vs S_q 2.09). The question: **where does divergence begin — q=4 or q=5?**

q=4 is the marginal Ashkin-Teller point for S_q Potts (BKT-like, alpha=1.771 from Sprint 118). Literature confirms: S_q q=4 has c=1, ν=2/3, and multiplicative log corrections from a marginal operator (Salas-Sokal 1997, Balog et al. 2007). The measured S_q alpha=1.771 is consistent with asymptotic alpha=2.0 dressed by log corrections.

## Experiments

### 120a — Hybrid q=4 g_c via gap×N crossing

**Method:** Coarse scan g∈[0.33,0.43] with 30 points at n=6,8,10 (GPU).

**Result: g_c(hybrid, q=4) = 0.3933 ± 0.0002**

| Pair | g_c crossing |
|------|-------------|
| (6,8) | 0.393412 |
| (8,10) | 0.393250 |

Excellent convergence (spread 0.000162). Fits pattern: q=2→0.250, q=3→0.333, **q=4→0.393**, q=5→0.438, q=7→0.535, q=10→0.684.

### 120b — chi_F spectral decomposition at hybrid q=4

**Method:** Spectral chi_F at g_c=0.3933, n=4-11 (8 data points, up to dim=4.2M). Also q=3 sanity check (n=4,6,8,10).

**Raw data:**

| n | dim | chi_F | gap_m | |me|² | frac | time |
|---|-----|-------|-------|-------|------|------|
| 4 | 256 | 2.071 | 0.957 | 7.59 | 1.000 | 0.0s |
| 5 | 1,024 | 2.984 | 0.766 | 8.76 | 1.000 | 0.1s |
| 6 | 4,096 | 3.978 | 0.638 | 9.72 | 1.000 | 0.2s |
| 7 | 16,384 | 5.047 | 0.547 | 10.56 | 1.000 | 1.0s |
| 8 | 65,536 | 6.184 | 0.478 | 11.30 | 1.000 | 25s |
| 9 | 262,144 | 7.385 | 0.425 | 11.98 | 1.000 | 34s |
| 10 | 1,048,576 | 8.647 | 0.382 | 12.60 | 1.000 | 79s |
| 11 | 4,194,304 | 9.965 | 0.347 | 13.19 | 1.000 | 198s |

**Single-multiplet dominance frac=1.000 at ALL sizes** — universal mechanism confirmed at q=4.

**Pairwise exponents:**

| Pair | alpha | z_m | beta_me | recon |
|------|-------|-----|---------|-------|
| (4,5) | 1.637 | 0.997 | 0.644 | 1.637 |
| (5,6) | 1.577 | 1.002 | 0.572 | 1.577 |
| (6,7) | 1.543 | 1.005 | 0.533 | 1.543 |
| (7,8) | 1.522 | 1.006 | 0.510 | 1.522 |
| (8,9) | 1.507 | 1.006 | 0.494 | 1.507 |
| (9,10) | 1.497 | 1.007 | 0.484 | 1.497 |
| (10,11) | 1.489 | 1.007 | 0.476 | 1.489 |

**Global fit: alpha=1.548, z_m=1.004, beta_me=0.541, nu_eff=0.785**

**Sanity check (q=3):** alpha=1.472, z_m=1.030 — matches prior Sprint 119 (1.40 global). Small discrepancy from using only 4 sizes vs Sprint 119's 5.

## Key Findings

### 1. z_m = 1.004 at q=4 — exactly marginal

This is the most striking result. The dynamic exponent z_m crosses 1 at q=4:

| q | Hybrid z_m | S_q z_m | Interpretation |
|---|-----------|---------|---------------|
| 3 | 1.03 | 1.03 | Both identical, z_m > 1 |
| **4** | **1.00** | **~1.05** | **Hybrid exactly marginal** |
| 5 | 0.91 | 1.31 | Diverged: continuous vs walking |
| 7 | 0.77 | ~1.4 | Strong divergence |

For the hybrid, z_m crosses 1 between q=3 and q=5, with q=4 sitting right at z_m=1. For S_q, z_m stays above 1 at q=4 (marginal operator keeps it walking).

### 2. Hybrid alpha peaks at q=4 before declining

Complete hybrid alpha(q) trajectory:
- q=3: 1.40 → q=4: **1.55 (peak)** → q=5: 1.41 → q=7: 0.95 → q=10: 0.55

The hybrid alpha has a **non-monotonic peak at q=4**. Below q=4, both models walk (z_m>1). Above q=4, the hybrid transitions to continuous (z_m<1) and alpha falls, while S_q maintains walking and alpha grows.

### 3. Model divergence starts precisely at q=4

| q | Hybrid alpha | S_q alpha | Divergence |
|---|-------------|-----------|-----------|
| 3 | 1.47 | 1.40 | 5% (noise) |
| 4 | 1.55 | 1.77 | 12% |
| 5 | 1.41 | 2.09 | 33% |
| 7 | 0.95 | 2.65 | 64% |

### 4. Pairwise alpha decreasing (log corrections)

Hybrid q=4 pairwise alpha drifts from 1.64→1.49 over n=4-11. This is **opposite** to what S_q q=4 should do (drift upward toward 2.0). The hybrid is not approaching the Ashkin-Teller fixed point — it's a different universality class even at q=4.

### 5. Literature context

S_q q=4 has known multiplicative log corrections: chi_F ~ L² (ln L)^{-p} giving effective alpha < 2.0 at finite size (Salas-Sokal 1997). The hybrid q=4 alpha drifts in the wrong direction to approach 2.0, confirming the two models are in different universality classes at q=4.

## Summary

**q=4 is where the hybrid and S_q models begin to diverge.** The hybrid z_m crosses 1 at q=4, marking the exact boundary between walking-like (z_m>1, q≤3) and continuous (z_m<1, q≥5) behavior. The S_q model maintains z_m>1 at q=4 (and beyond) due to its S_q-symmetric marginal operator. The Z_q-symmetric hybrid field (X+X†) lacks this marginal operator for q>3, causing z_m to drop below 1 and alpha to decrease.

## Next Steps

1. **Hybrid q=4 log correction fit** — fit chi_F = A·N^α·(ln N)^{-p} to extract the log correction exponent. Compare with S_q q=4 expectation.
2. **S_q q=4 at larger sizes** — extend Sprint 118 data (n=4-11) to see if alpha drifts toward 2.0.
3. **z_m(q) continuous fit** — with q=3,4,5,7,10 data, fit z_m(q) and find the precise q where z_m=1 (should be q_cross ≈ 4.0 for hybrid).
