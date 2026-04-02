# Sprint 100 — DMRG Casimir Energy: Open vs Periodic BC

**Status:** Complete (3 experiments).

## Motivation

Our confirmed novel finding (Sprint 098): Casimir energy tracks Re(c) of complex CFT across the walking boundary, 16x more consistent than entropy. Established with periodic-BC exact diag only (N=4-14 for q=2, N=4-10 for q=5, N=4-8 for q=7).

**Key question:** Can DMRG extend this to larger N via open-BC Casimir extraction?

**Literature:** Carlon & Igloi (PRB 57, 1998) extracted Casimir for q=2,3,4 Potts via DMRG. No paper reports successful Im(c) oscillation detection from real Hamiltonian Casimir energy — the oscillation period is ~e^300 in system size (Sprint 099 finding confirmed by literature).

## Experiments

### 100a: DMRG Casimir precision test — q=2 (known c=0.5), open BC

**7 DMRG data points** (N=8,10,12,14,16,20,24). DMRG agrees with exact diag to ~10⁻¹⁴ at N≤14. Open BC formula: E₀ = ε_∞·N + e_s - πvc/(24N) + B/N².

| Method | c | c/0.5 | Notes |
|--------|---|-------|-------|
| 4-param fit | 0.485 ± 0.003 | 0.971 | All 7 sizes, R²=1.0 |
| Richardson extrap | 0.486 | 0.972 | Pairwise extrapolation |
| Pairwise (20,24) | 0.467 | 0.934 | Converging monotonically |

**Boundary corrections ~3% at N~24.** Pairwise c/0.5 increases monotonically: 0.88 → 0.93, converging slowly.

### 100b: Open-BC Casimir for q=5,7 — NEGATIVE

Applied same method to q=5 (exact diag N=6-10 + DMRG N=12) and q=7 (exact diag N=4-8).

| q | Points | N range | c(4p) | c/Re(c) | c/c_eff | Verdict |
|---|--------|---------|-------|---------|---------|---------|
| 2 | 7 | 8-24 | 0.485 | 0.971 | 0.971 | OK (3% off) |
| 5 | 6 | 6-12 | 0.729 | 0.640 | 0.633 | **FAILS** (36% off) |
| 7 | 5 | 4-8 | 0.483 | 0.357 | 0.456 | **FAILS** (64% off) |

**Root cause:** Open-BC Casimir formula has 1/N boundary corrections that scale with q. At q=5 (d=5 local dim), boundary corrections are ~5× larger than q=2. At q=7 (~7×). The small accessible sizes (N≤12 for q=5, N≤8 for q=7) mean boundary terms completely dominate the Casimir signal.

**Fundamental limitation:** Open-BC Casimir extraction requires N ≫ 1/Δ_boundary where Δ_boundary depends on the boundary operator content. For S_q Potts with free boundaries, Δ_boundary shrinks with q, requiring exponentially larger N.

### 100c: Periodic-BC Casimir reanalysis — confirms Re(c) tracking

Reanalyzed Sprint 099a dense periodic-BC data (all integer N). The periodic formula E₀/N = ε_∞ - πvc/(6N²) has 1/N² corrections (not 1/N), making it far more robust.

**Pairwise c/Re(c) from consecutive periodic-BC sizes:**

| q | (N₁,N₂) | c/Re(c) | Best pair | Notes |
|---|---------|---------|-----------|-------|
| 2 | (5,6) | 1.000 | **(13,14): 0.984** | 11 pts, convergent |
| 5 | (7,8) | 1.002 | **(9,10): 0.994** | 7 pts, convergent |
| 7 | (5,6) | 1.000 | **(7,8): 0.974** | 5 pts |

**2-param fit c/Re(c):** q=2: 0.999, q=5: 1.019, q=7: 1.006.

**Periodic BC gives c/Re(c) ≈ 1.00 for ALL q at modest sizes.** The 1/N² scaling of periodic corrections vs 1/N for open BC makes periodic extraction 10-100× more accurate at same N.

## Key Findings

1. **Open-BC Casimir extraction FAILS for q≥5 S_q Potts** at accessible sizes. Boundary corrections scale with q and overwhelm the Casimir signal. Would need N≥50 for q=5, N≥100 for q=7 — impractical for DMRG with d≥5.

2. **Periodic-BC Casimir extraction SUCCEEDS at modest sizes.** Sprint 098/099 periodic data gives c/Re(c) = 1.00 ± 0.02 for all q=2-7. The confirmed novel finding stands.

3. **The Casimir-Re(c) result CANNOT be extended to large N via DMRG.** Periodic DMRG is too expensive for q≥5, and open-BC extraction fails. The periodic exact diag sizes (N≤14 for q=2, N≤10 for q=5) are the practical limit.

4. **Literature confirms:** No paper has detected Im(c) oscillations from Casimir energy at any accessible size. The period is astronomically long (~e^300 for q=5). The non-Hermitian deformation approach (Tang et al. 2024, Shimizu & Kawabata 2025) is the only viable route to access Im(c).

## Surprises

- Open-BC boundary corrections scale dramatically with q (3% at q=2, 64% at q=7 for same N range)
- Periodic BC is 10-100× more accurate than open BC at same N — the 1/N vs 1/N² scaling difference is enormous
- The confirmed novel finding (Casimir tracks Re(c)) is robust but size-limited — it rests on periodic exact diag at N≤10-14

## What This Changes

The Casimir-Re(c) result is confirmed but cannot be extended to DMRG sizes. The next frontier requires either: (1) non-Hermitian formulation, (2) transfer matrix approach for periodic BC at large sizes, or (3) a completely different observable.

[Data: results/sprint_100a_dmrg_casimir_q2.json, results/sprint_100b_dmrg_casimir_q5q7.json, results/sprint_100c_casimir_analysis.json]
