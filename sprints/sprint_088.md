# Sprint 088 — Tail Weight Scaling at q=2,3: UNIVERSAL Power Law, Not Walking-Specific

**Status:** Complete (3 experiments).

**Question:** Sprint 087 found unbounded power-law tail weight growth for q=5 (walking) and q=7 (broken). Is this unique to walking/complex CFT (q>4), or does it occur for real CFT too (q=2,3)?

**Answer: UNIVERSAL.** All q=2,3,5 have identical exponent b ≈ 2.0 in log-log space. Power-law tail growth is a generic feature of critical entanglement spectra, not a walking signature.

---

## Experiment 088a — q=3 entanglement spectrum n=8-24

DMRG at g_c=1/3 (c=4/5 real CFT). Open BC, midchain Schmidt spectrum.

| n | chi | λ_max | w_mult | w_tail | %S(l0) | %S(l1) | %S(tail) | Δξ | S |
|---|-----|-------|--------|--------|--------|--------|----------|------|------|
| 8 | 40 | 0.860 | 0.139 | 0.00047 | 0.257 | 0.735 | 0.008 | 2.515 | 0.504 |
| 12 | 60 | 0.833 | 0.165 | 0.00132 | 0.264 | 0.718 | 0.018 | 2.311 | 0.574 |
| 16 | 80 | 0.815 | 0.182 | 0.00230 | 0.269 | 0.705 | 0.027 | 2.190 | 0.620 |
| 20 | 100 | 0.802 | 0.195 | 0.00330 | 0.271 | 0.694 | 0.035 | 2.107 | 0.654 |
| 24 | 120 | 0.791 | 0.205 | 0.00429 | 0.273 | 0.686 | 0.042 | 2.043 | 0.681 |

(q-1)=2 multiplet perfectly degenerate (spread = 0.000000 at all sizes).

## Experiment 088b — q=2 entanglement spectrum n=8-24

DMRG at g_c=1/2 (c=1/2 Ising). Open BC, midchain Schmidt spectrum.

| n | chi | λ_max | w_mult | w_tail | %S(l0) | %S(l1) | %S(tail) | Δξ | S |
|---|-----|-------|--------|--------|--------|--------|----------|------|------|
| 8 | 16 | 0.890 | 0.110 | 0.00037 | 0.297 | 0.694 | 0.009 | 2.091 | 0.350 |
| 12 | 30 | 0.870 | 0.129 | 0.00102 | 0.309 | 0.672 | 0.019 | 1.906 | 0.393 |
| 16 | 37 | 0.856 | 0.142 | 0.00176 | 0.315 | 0.657 | 0.028 | 1.798 | 0.422 |
| 20 | 45 | 0.846 | 0.151 | 0.00250 | 0.319 | 0.645 | 0.036 | 1.724 | 0.443 |
| 24 | 54 | 0.839 | 0.158 | 0.00324 | 0.321 | 0.635 | 0.044 | 1.668 | 0.459 |

## Experiment 088c — Cross-q comparison and power-law fits

### Power-law exponents (log-log fit: w_tail ~ n^b)

| q | b (log-log) | b (real-space) | A (log-log) | R² (real) | n range | CFT type |
|---|------------|----------------|-------------|-----------|---------|----------|
| 2 | **1.98** | 1.69 | 6.7e-6 | 0.992 | 8-24 | real (Ising, c=1/2) |
| 3 | **2.01** | 1.72 | 7.9e-6 | 0.993 | 8-24 | real (3-Potts, c=4/5) |
| 5 | **2.01** | 1.72 | 8.1e-6 | 0.993 | 8-24 | walking (complex, Re(c)=1.138) |
| 7 | **2.97** | 2.52 | 7.5e-7 | 0.983 | 6-12 | broken (complex, Re(c)=1.351) |

**q=2, 3, and 5 have IDENTICAL exponent b ≈ 2.0** in log-log space. q=7 appears higher (b≈3.0) but has only 4 points over a narrow range — likely contaminated by pre-asymptotic corrections.

### Methodological correction to Sprint 087

Sprint 087 reported b=1.72 (q=5) and b=2.52 (q=7) using **real-space** curve_fit, which overweights large-n points. The standard approach for power-law scaling is log-log regression (equal relative weighting). Log-log gives b≈2.0 for q=2,3,5 — identical within error. The Sprint 087 conclusion of "unbounded power law" remains correct; only the exponent values need correction.

### Entanglement gap: Δξ = a + b·ln(n), all q

| q | a | b | R² | Δξ(8) | Δξ(24) |
|---|---|---|---|-------|--------|
| 2 | 2.88 | -0.39 | 0.993 | 2.09 | 1.67 |
| 3 | 3.39 | -0.43 | 0.994 | 2.51 | 2.04 |
| 5 | 4.07 | -0.43 | 0.994 | 3.19 | 2.71 |
| 7 | 4.78 | -0.50 | 0.994 | 3.90 | 3.55 |

Logarithmic closure is also universal. Intercept a grows with q (higher initial gap), slope b ≈ -0.43 for q=2-5, slightly steeper at q=7.

### Entropy fractions at n=24

| q | %S(lev0) | %S(lev1) | %S(tail) | λ_max | w_tail |
|---|----------|----------|----------|-------|--------|
| 2 | 0.321 | 0.635 | 0.044 | 0.839 | 0.00324 |
| 3 | 0.273 | 0.686 | 0.042 | 0.791 | 0.00429 |
| 5 | 0.225 | 0.737 | 0.038 | 0.787 | 0.00437 |

At matched sizes (n=24), %S(tail) is nearly identical across q=2,3,5 (3.8-4.4%). The walking regime is invisible at this level.

### Extrapolated n where w_tail = 10%

| q | n(10%) | method |
|---|--------|--------|
| 2 | ~130 | log-log extrapolation |
| 3 | ~109 | log-log extrapolation |
| 5 | ~108 | log-log extrapolation |
| 7 | ~53 | log-log extrapolation (4 pts, uncertain) |

---

## Key Findings

1. **Tail weight power law is UNIVERSAL** — q=2,3,5 all show w_tail ~ n^2.0 with identical exponent and similar prefactors. This is a generic property of critical entanglement spectra, not a signature of walking or complex CFT.

2. **Sprint 087 exponents were fitting artifacts.** Real-space curve_fit gives biased exponents that vary with q. Log-log fit (standard for power laws) reveals universal b≈2.0 for q=2,3,5.

3. **q=7 exponent appears larger (b≈3.0)** but only 4 points over n=6-12. This could be pre-asymptotic corrections or genuinely faster growth in the walking-broken regime. Needs n>12 data to resolve.

4. **Entropy fractions at matched sizes are q-independent.** %S(tail) ≈ 4% at n=24 for all q=2,3,5. The walking regime does NOT manifest in tail weight at accessible sizes.

5. **Entanglement gap closure is universal.** Δξ ~ -0.43·ln(n) for q=2-5. All critical spectra have logarithmically closing gaps.

6. **What IS q-dependent:** the initial gap intercept (grows with q as ~0.47+0.78·ln(q)), and the multiplet weight (grows with q due to (q-1) degeneracy). The RATE of tail growth is universal.

## Surprises

- q=2,3,5 tail exponents identical to <1% — far more universal than expected
- The "energy-entropy decoupling" hierarchy from Sprint 087 is reframed: ALL critical points have the same tail growth. Walking breakdown must operate through the entanglement TEMPERATURE (spacing between levels), not the tail weight.
- Sprint 087's n≈147 (q=5) vs n≈72 (q=7) crossover scales are correct in spirit but the mechanism is different: same tail weight growth rate, different sensitivity of observables to that tail.

## What This Changes

The narrative from Sprints 082-087 was: "walking breakdown = entropy concentration = tail weight growth." But tail weight grows identically for real CFT (q=2,3) where there IS no walking breakdown. **The walking-specific phenomenon must be the REDISTRIBUTION between levels** (how multiplet and ground state share weight), not the absolute tail growth rate.

Specifically: at q=5-7, %S(lev1) is higher (74% vs 64% for q=2) because the (q-1)-fold multiplet has more states. This leaves less entropy "budget" for the ground state level. The tail draws from the same universal pool, but the consequences differ because the starting distribution is q-dependent.

[Full data: results/sprint_088a_entspec_dmrg_q3.json, results/sprint_088b_entspec_dmrg_q2.json, results/sprint_088c_tail_comparison.json]
