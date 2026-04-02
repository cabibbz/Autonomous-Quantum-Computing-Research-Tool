# Sprint 086 — Rényi Entropy Scaling via DMRG: α=3 Was an Extraction Artifact

**Date:** 2026-04-02
**Status:** Complete (3 experiments)

## Motivation

Sprint 085 discovered that α=3 Rényi entropy uniquely recovers Re(c) in the walking-broken regime (q=7-8), while α=1 (von Neumann) is most accurate for walking q≤6. Those measurements used **single-size extraction**: c_α = 12·S_α / ((1+1/α)·ln(N/π)), which includes the non-universal constant c'_α. This sprint tests whether that finding survives at larger sizes via DMRG, and whether **size-pair extraction** (which cancels c'_α) gives the same answer.

## Experiments

### 086a — DMRG Rényi entropies at q=5, n=8,10,12

Open-BC DMRG at g_c=1/5. Extracted Schmidt spectrum at every bond, computed S_α for α=0.5,1,2,3,5,10,∞. Profile fits to CC formula gave poor R² (0.79-0.88) due to open-BC boundary corrections. Key data is midchain S_α for size-pair extraction.

**Midchain S_α (q=5, DMRG open BC):**

| n | chi | S₀.₅ | S₁ | S₂ | S₃ | S_∞ |
|---|-----|------|-----|-----|-----|------|
| 8 | 60 | 1.106 | 0.610 | 0.300 | 0.230 | 0.153 |
| 10 | 80 | 1.167 | 0.660 | 0.335 | 0.257 | 0.172 |
| 12 | 60 | 1.216 | 0.699 | 0.362 | 0.279 | 0.186 |

### 086b — DMRG Rényi entropies at q=7, n=6,8,10,12

Open-BC DMRG at g_c=1/7. q=7 has d=7 per site, so DMRG is slower (226s at n=10 chi=40).

**Midchain S_α (q=7, DMRG open BC):**

| n | chi | S₀.₅ | S₁ | S₂ | S₃ | S_∞ |
|---|-----|------|-----|-----|-----|------|
| 6 | 40 | 1.162 | 0.540 | 0.228 | 0.173 | 0.115 |
| 8 | 50 | 1.251 | 0.609 | 0.267 | 0.203 | 0.135 |
| 10 | 40 | 1.307 | 0.655 | 0.294 | 0.224 | 0.149 |
| 12 | 30 | 1.346 | 0.690 | 0.314 | 0.239 | 0.160 |

**Entanglement spectrum confirms (q-1)=6 fold degeneracy at all sizes:**
- n=6: λ_max=0.891, 6×λ₁=0.109, tail=0.01%
- n=12: λ_max=0.852, 6×λ₁=0.147, tail=0.11%

Tail weight grows 8× from n=6→12, directly tracking walking breakdown.

### 086c — Three-way comparison of c_α extraction methods

Compared:
1. **Single-size periodic** (Sprint 085): c_α = 12·S / ((1+1/α)·ln(N/π))
2. **Size-pair periodic** (Sprint 085 data): c_α = 6·ΔS / ((1+1/α)·ln(N₂/N₁))
3. **Size-pair open DMRG** (Sprint 086 data): c_α = 12·ΔS / ((1+1/α)·ln(L₂/L₁))

**Best α by method and q:**

| q | single(periodic) | pair(periodic) | pair(open DMRG) |
|---|-----------------|----------------|-----------------|
| 2 | α=0.5 | α=0.5 | — |
| 3 | α=1 | α=0.5 | — |
| 5 | α=1 | α=1 | α=0.5 |
| 6 | α=1 | α=1 | — |
| 7 | **α=3** | **α=2** | **α=1** |
| 8 | **α=3** | **α=2** | — |

**The α=3 result is specific to single-size extraction.** Size-pair extraction (both periodic and open BC) gives α=1 or α=2 as best for q=7-8.

**Periodic size-pair c_α/Re(c) at best pair:**

| q | pair | c₁/Rec | c₂/Rec | c₃/Rec | c_∞/Rec | best |
|---|------|--------|--------|--------|---------|------|
| 2 | (10,14) | 0.996 | 1.055 | 1.114 | 1.101 | α=0.5 (1.000) |
| 3 | (8,10) | 0.994 | 1.081 | 1.131 | 1.080 | α=0.5 (1.001) |
| 5 | (6,8) | 1.004 | 1.114 | 1.127 | 1.040 | α=1 (1.004) |
| 6 | (6,8) | 0.997 | 1.102 | 1.098 | 1.006 | α=1 (0.997) |
| 7 | (6,7) | 0.817 | 0.904 | 0.889 | 0.809 | α=2 (0.904) |
| 8 | (6,7) | 0.799 | 0.877 | 0.852 | 0.772 | α=2 (0.877) |

**Walking breakdown in DMRG c₁ (q=7, adjacent pairs):**

| pair | c₁ | c₁/Re(c) |
|------|-----|----------|
| (6,8) | 1.425 | 1.055 |
| (8,10) | 1.252 | 0.927 |
| (10,12) | 1.127 | 0.834 |

c₁ drops 22% from (6,8) to (10,12) — clear walking breakdown. For q=5, c₁ drops only 5% over the same range.

## Key Findings

1. **Sprint 085's α=3 finding was a SINGLE-SIZE EXTRACTION ARTIFACT.** The non-universal constant c'_α varies with α. At n=6-7 for q=7, c'₃ happens to partially compensate the walking deviation, making c₃ look closest to Re(c). Size-pair extraction (which cancels c'_α) shows α=1-2 is consistently best.

2. **Optimal α shifts gradually: 0.5 (real CFT) → 1 (walking) → 2 (broken).** This is a smoother progression than Sprint 085 suggested. The shift reflects increasing weight redistribution in the entanglement spectrum.

3. **No Rényi index recovers Re(c) for broken walking.** At q=7, even the best α gives only 90% of Re(c) from periodic pairs and degrades with n. Walking breakdown is a genuine multi-eigenvalue phenomenon — no single α can compensate.

4. **Open-BC DMRG profile fits unreliable for c_α extraction (R² < 0.9).** Boundary corrections dominate at n≤12. Size-pair extraction from midchain S_mid is more robust.

5. **Walking breakdown clearly visible in DMRG size-pair c₁.** q=7 c₁ drops 22% over n=6-12, while q=5 drops only 5%. This is the most direct DMRG evidence of walking breakdown in Rényi entropies.

## Corrections to Sprint 085

- ~~α=3 uniquely recovers Re(c) in walking-broken regime~~ → α=3 appeared optimal only due to single-size extraction contamination by c'_α
- ~~Rényi spread (c₂-c_∞) is a new monotonic walking discriminator~~ → This was computed from single-size extraction and is contaminated by c'_α differences. Size-pair Rényi spread may still be a valid discriminator but needs re-measurement.
- Sprint 085's α=1 finding for walking q≤6 IS confirmed by size-pair extraction

## Surprises

- Open-BC profile fit R² < 0.9 even at n=12 — boundary corrections much larger than expected
- q=7 DMRG entanglement tail grows 8× from n=6→12 while λ_max drops only 4%
- Periodic and open BC give DIFFERENT optimal α (periodic: α=2, open: α=1 for q=7) — BC matters for Rényi structure
- Single-size vs size-pair extraction gives qualitatively different conclusions about optimal α

[Full data: results/sprint_086a_dmrg_renyi_q5.json, results/sprint_086b_dmrg_renyi_q7.json, results/sprint_086c_renyi_analysis.json]
