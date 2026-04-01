# Sprint 045 — ν(q) Extraction: Does Universality Class Change with q?

## Motivation
We know ν=1 for q=2 (Ising) and ν≈5/6 for q=3 (3-state Potts). For q>4, 1D quantum Potts is always second-order (Sprint 043) but the universality class is unknown. The quantum-classical correspondence maps 1D quantum to an ANISOTROPIC 2D classical model — not the standard isotropic 2D Potts where q>4 is first-order. Extracting ν(q) at q=5 will reveal whether there's a q-dependent family of universality classes.

## Predictions
1. ν should decrease with q (larger q → more mean-field-like → ν closer to 1/2)
2. If the exponent 0.85 in g_c(q) is exactly 5/6 = ν(q=3), there may be a deep connection
3. q=7 and q=10 may share the same ν (large-q universality class)

## Experiments

### 045a — q=7 Data Collapse with Existing n=8,12 Data
**Inconclusive.** Only 3 overlapping g-values (0.25, 0.30, 0.35), all on the ordered/near-critical side. The collapse quality landscape is flat and the free fit gives physically unreasonable ν=0.44. Need more data (but n=12 DMRG at d=7 takes >120s/point, too slow for dense sweeps).

### 045b — q=5 Potts n=8 Dense Sweep (chi=20)
8 g-values: 0.30-0.55, ~5s/point. CV monotonically decreasing from 0.373 (g=0.30) to 0.107 (g=0.50), then rising to 0.132 (g=0.55). The minimum around g=0.50 is the MI-CV minimum between ordered and deep-disordered phases.

### 045c — q=5 Potts n=12 Dense Sweep (chi=20)
4 new g-values (0.35, 0.38, 0.42, 0.45) filling gaps from Sprint 042d.

### 045d — Verify Sprint 042d Data + n=16
**CRITICAL BUG FOUND IN SPRINT 042.** Sprint 042d used Gell-Mann correlation reconstruction for MI, which gives WRONG results at q=5:
- g=0.30: 042d=0.134, actual=0.353 (62% error)
- g=0.40: 042d=0.516, actual=0.299 (42% error)
- g=0.50: 042d=1.360, actual=0.121 (11x error!)

The Gell-Mann reconstruction is fundamentally unreliable for q≥5 — numerical errors in the 24×24 correlation matrix accumulate. Direct MPS contraction (Sprint 043) is the only reliable method for d≥5.

**Implication:** Sprint 042's crossing analysis mixed χ=15 (n=8) and χ=20 (n=12) AND used bad Gell-Mann MI. The reported g_c≈0.41 for q=5 was coincidentally close to the true value but for wrong reasons.

**n=16 data obtained** at g=0.30, 0.35, 0.38, 0.40, 0.42, 0.45, 0.50 (~20-60s/point).

### 045f — q=5 Data Collapse (n=8, 12, 16)
**ν(q=5) ≈ 2.0, g_c ≈ 0.45.** This is the headline result.

Data collapse with joint (ν, g_c) optimization across n=8,12,16:

| Model | ν | g_c | Quality | Ratio |
|-------|---|-----|---------|-------|
| **Free fit** | **2.24** | **0.446** | **0.000265** | **1.00** |
| ν=2 | 2.00 | 0.442 | 0.000276 | 1.04 |
| ν=3/2 | 1.50 | 0.435 | 0.000421 | 1.59 |
| ν=1 (Ising) | 1.00 | 0.429 | 0.000980 | 3.69 |
| ν=5/6 (Potts q=3) | 0.833 | 0.427 | 0.001528 | 5.76 |
| ν=2/3 | 0.667 | 0.424 | 0.002754 | 10.38 |
| ν=1/2 | 0.500 | 0.421 | 0.005631 | 21.22 |

**Independent confirmation from slope scaling:** MI-CV slope at transition scales as n^0.487, giving ν = 1/0.487 = 2.05.

Crossing points: n=(8,12) g_c=0.332, n=(12,16) g_c=0.346. The large finite-size shift (crossing 0.33-0.35 vs thermodynamic g_c≈0.45) is EXPECTED for large ν — correlation length diverges slowly, so finite-size effects are amplified.

## Key Results

### Prediction 1 was WRONG
ν INCREASES with q, opposite to mean-field prediction:
| q | ν | 1/ν | Trend |
|---|---|-----|-------|
| 2 | 1.0 | 1.0 | Ising |
| 3 | 5/6 ≈ 0.83 | 1.2 | 3-state Potts |
| 5 | ~2.0 | ~0.5 | **NEW** |

The universality class changes dramatically between q=3 and q=5. ν goes DOWN from q=2 to q=3, then UP sharply at q=5. This non-monotonic behavior suggests the self-duality breaking at q=3→4 fundamentally alters the transition character.

### Physical interpretation
Large ν means the transition is very gentle — the correlation length diverges slowly as ξ ~ |g-g_c|^(-2). This is consistent with:
1. The MI-CV curves changing gradually (not sharply) for q=5
2. The enormous finite-size shift in crossing points
3. The suppression of the first-order mechanism in 1D quantum

In 2D classical Potts, q>4 is first-order (ν→∞). In 1D quantum, the extreme anisotropy of the quantum-classical mapping prevents the transition from going first-order, but the TENDENCY toward first-order manifests as anomalously large ν. The transition is "almost first-order" — second-order but with a diverging correlation length exponent.

### Gell-Mann MI reconstruction: BROKEN for q≥5
Sprint 043's direct MPS contraction method is essential for d≥5. The Gell-Mann method has uncontrolled numerical errors. This invalidates Sprint 042's quantitative results for q=5 (but the qualitative conclusion — second-order, not first-order — remains valid since it was confirmed with direct MPS in Sprint 043).

## Surprises
1. **ν=2 at q=5** — 2x the Ising value, opposite to mean-field prediction
2. **Non-monotonic ν(q)**: decreases q=2→3, then increases sharply q=3→5
3. **Sprint 042d data was wrong by up to 11x** — Gell-Mann reconstruction fails for q=5
4. **Finite-size crossing points 30% below thermodynamic g_c** — large ν amplifies finite-size effects
5. **g_c ≈ 0.45** for q=5 — higher than Sprint 042's estimate of 0.41

## Files
- exp_045a_q7_collapse.py — q=7 collapse (inconclusive)
- exp_045b_potts_q5_n8_dense.py — n=8 dense sweep
- exp_045c_potts_q5_n12_dense.py — n=12 dense sweep
- exp_045d_verify_and_n16.py — verify 042d + n=16 data
- exp_045e_verify_n12_g030.py — verify n=12 g=0.30 + more n=16
- exp_045f_q5_collapse.py — data collapse analysis
- results/sprint_045*.json — all raw data
