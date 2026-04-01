# Sprint 051 — Energy Gap Method to Locate g_c for q=4,5 Potts

**Date:** 2026-04-01
**Status:** Complete (5 experiments)

## Motivation

After Sprint 049 revealed all Potts g_c values for q≥3 were wrong, and Sprint 050 used self-duality to get exact g_c for q=2 (0.25) and q=3 (1/3), we still lack g_c for q≥4 where self-duality is broken. The energy gap Δ = E₁ - E₀ closes as 1/N at criticality (CFT), so Δ·N is scale-invariant at g_c. Curves of Δ·N vs g for different N cross at g_c. This is cleaner than entropy (log singularity) or MI-CV (dead pairs, χ issues).

**Literature:** q=3 Potts has c=4/5, q=4 is related to Ashkin-Teller with c=1. The energy gap method is standard in numerical studies of quantum phase transitions.

## Experiments

### Exp A: Validate energy gap on q=2,3 (known g_c)
**Prediction:** Δ·N curves for n=4,6,8 cross at g_c=0.25 (q=2) and g_c=1/3 (q=3).

### Exp B: q=4 energy gap → g_c
**Prediction:** Crossing somewhere in [0.30, 0.40] based on Sprint 050 pseudo-critical estimate.

### Exp C: q=5 energy gap → g_c
**Prediction:** Unknown — previous g_c≈0.41 was at wrong point.

---

## Results

### Exp A Results: Energy gap validated on q=2,3

**Method works.** Δ·N curves for different n cross near the known g_c, with systematic finite-size convergence.

**q=2 (g_c = 0.25 exact):**
| n₁, n₂ | g_cross |
|---------|---------|
| 4, 6   | 0.2386  |
| 4, 8   | 0.2413  |
| 6, 8   | 0.2439  |
| 6, 10  | 0.2450  |
| 8, 10  | 0.2462  |

Largest-size crossing (8,10) at 0.246 — 1.5% below exact. Converging monotonically upward.

**q=3 (g_c = 1/3 ≈ 0.333 exact):**
| n₁, n₂ | g_cross |
|---------|---------|
| 4, 6   | 0.3178  |
| 4, 8   | 0.3215  |
| 6, 8   | 0.3251  |

Largest-size crossing (6,8) at 0.325 — 2.5% below exact. Same monotonic convergence.

**Finite-size correction pattern:** crossings approach g_c from below. Rate: ~2-3% error for largest available pair. For q=4,5 where g_c is unknown, the (6,8) crossing should be ~2-3% below true g_c.

### Exp B Results: q=4 energy gap → g_c ≈ 0.39

**Crossings converge monotonically:**
| n₁, n₂ | g_cross |
|---------|---------|
| 4, 6   | 0.3728  |
| 4, 8   | 0.3777  |
| 6, 8   | 0.3823  |

Best estimate: g_c(q=4) ≈ 0.382 raw, ~0.392 corrected. Sprint 050 pseudo-critical at g≈0.34 (n=8) was 13% too low.

### Exp C Results: q=5 energy gap → g_c ≈ 0.44

| n₁, n₂ | g_cross |
|---------|---------|
| 4, 6   | 0.4199  |
| 6, 8   | 0.4303  |

Best estimate: g_c(q=5) ≈ 0.430 raw, ~0.441 corrected. Previous Sprint 042 value g_c≈0.41 was close but ~7% low.

q=5 n=8 (dim=390,625) took ~7s/point. Full scan feasible but timed out at g=0.48; targeted scan near crossing worked.

### Exp D Results: q=7 energy gap → g_c ≈ 0.52

No crossings in [0.05, 0.35]! Extended range found:

| n₁, n₂ | g_cross |
|---------|---------|
| 4, 6   | 0.5106  |

Best estimate: g_c(q=7) ≈ 0.511 raw, ~0.524 corrected. The old value g_c≈0.26 was off by 2×!

### Exp E: q=10 feasibility

q=10 n=4 (dim=10,000): 0.1s/point — full scan done. q=10 n=6 (dim=1,000,000): 23.5s/point — too slow for full scan. Would need DMRG for n=6+.

## Key Discovery: g_c INCREASES with q

The old (invalidated) trend was g_c decreasing: 1.0 → 1.0 → 0.89 → 0.41 → 0.25. The TRUE trend is g_c **increasing**:

| q | g_c (exact/best estimate) | Method |
|---|--------------------------|--------|
| 2 | 0.250 | Self-duality (exact) |
| 3 | 0.333 | Self-duality (exact) |
| 4 | ~0.39 | Gap crossing (n=6,8) + 2.5% correction |
| 5 | ~0.44 | Gap crossing (n=6,8) + 2.5% correction |
| 7 | ~0.52 | Gap crossing (n=4,6) + 2.5% correction |

**Physical interpretation:** More ground-state degeneracy (q-fold) requires a STRONGER transverse field to destroy order. The field term X + X† creates transitions between adjacent states only (|s⟩ → |s±1⟩), so it takes longer to "mix" a q-fold degenerate ground state as q increases. g_c appears to grow sublinearly with q.

**This COMPLETELY reverses the old picture.** The invalidated scaling law g_c ∝ (q-3)^{-0.85} predicted g_c→0 for large q. The true trend is g_c→∞ (or saturation).

## What went wrong with previous g_c estimates

The old approach (MI-CV crossings from DMRG, Sprints 038-048) measured disordered-phase crossovers, not the true critical point. The energy gap method directly probes the vanishing of the excitation gap — the defining feature of a quantum phase transition. It has:
- No dead-pair bias
- No χ convergence issues (exact diag)
- Clean convergence with system size
- Validation against known exact results

## Surprises

1. **g_c INCREASES with q** — opposite to all 10 sprints of (wrong) Potts results
2. **q=7 g_c ≈ 0.52, not 0.26** — off by factor of 2 from old estimate
3. **Energy gap method converges faster** than entropy or MI-CV methods (~2% error at n=6,8 vs ~30% for MI-CV)
4. **q=5 old estimate g_c≈0.41 was closest** to truth (~0.44) — perhaps MI-CV was least wrong at q=5
5. **Physical mechanism clear**: q-fold degeneracy is harder to break with nearest-neighbor X+X†
