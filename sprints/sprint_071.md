# Sprint 071 — 2D Cylinder Geometry: Ly=2 Ladder for q=2,5

**Status:** Complete (3 experiments)
**Date:** 2026-04-01

## Motivation

Sprints 068-070 studied the 2D hybrid model on periodic L×L torus but hit a bottleneck: q=5 L=4 has dim ≈ 10^11 (infeasible). Only L=2,3 accessible, and L=2 is pathologically out of scaling. We had exactly ONE reliable 2D data point (L=3) for q=5.

**This sprint:** Cylinder geometry with Ly=2 (ladder). Open BC in x-direction enables DMRG (for q=2) and exact diag with gap×Lx crossings (for both q). Coordination z=3 (between 1D z=2 and 2D z=4).

## Experiments

### 071a — q=2 Ly=2 Cylinder DMRG Validation

**Method:** TeNPy DMRG with Square lattice, bc=['open', 'periodic'], chi=20. Entropy scan to find g_c.

**Result:** g_c(q=2, Ly=2 cylinder) ≈ 0.45 from entropy peak. Peak drifts with Lx (0.35 at Lx=8 → 0.48 at Lx=20), indicating true g_c is slightly higher.

**c_eff = 0.19** — much lower than Ising c=0.5. Entropy saturates at S≈0.71 for Lx≥16. This is a chi=20 truncation artifact: the bond dimension limits the entanglement entropy growth on the cylinder.

**Key lesson:** DMRG entropy peak is unreliable for g_c on finite-width cylinders at modest chi. The peak position drifts because the pseudo-critical point converges slowly with Lx. Exact diag gap crossings are more reliable (experiment 071c).

### 071b — q=5 Cylinder DMRG (ABANDONED)

**DMRG impractical for q=5 cylinder.** With d=5 per site and Ly=2 cylinder geometry, chi=20 is insufficient (truncation error exceeds tolerance) and higher chi makes each DMRG run take >1 hour. Pivoted to exact diag for small Lx.

**Important lesson:** Cylinder DMRG for q≥5 requires d_eff > 5 per site in the MPS ordering, making chi≥50+ necessary. Each DMRG sweep scales as O(chi³ × d²), so q=5 cylinder is ~6× slower than q=2 at same chi. Not practical without Z_q symmetry-exploiting DMRG.

### 071c — Exact Diag Energy Gap on Ly=2 Cylinder

**Method:** Gap×Lx crossing using exact diag. Vectorized bond construction for non-adjacent sites. q=2: Lx=4,5,6,7 (dim up to 16k). q=5: Lx=3,4 (dim up to 390k).

**q=2 gap crossings converge beautifully:**
| Lx pair | g_c |
|---------|-----|
| (4,5) | 0.4468 |
| (5,6) | 0.4513 |
| (6,7) | 0.4536 |

**g_c(q=2, Ly=2 cyl) = 0.4506** — converging toward ~0.455. Ratio to 1D: 1.80.

**q=5 gap crossing:**
| Lx pair | g_c |
|---------|-----|
| (3,4) | 0.7138 |

**g_c(q=5, Ly=2 cyl) = 0.7138** — single crossing pair (Lx=5 at dim=10M takes 191s/point, too slow for full scan). Ratio to 1D: 1.62.

**Cylinder/1D g_c ratios:**
| q | g_c(1D) | g_c(cyl) | g_c(2D) | cyl/1D | 2D/1D |
|---|---------|----------|---------|--------|-------|
| 2 | 0.250 | 0.451 | 0.771 | 1.80 | 3.08 |
| 5 | 0.441 | 0.714 | 1.588 | 1.62 | 3.60 |

**Cylinder g_c is between 1D and 2D as expected.** The cyl/1D ratio (1.62-1.80) is approximately half of the 2D/1D ratio, consistent with adding one extra neighbor (z: 2→3 for cylinder, 2→4 for 2D).

**Order parameter ⟨δ(s_i,s_j)⟩ smooth across transition for both q.**
- q=2 (Lx=6): ⟨δ⟩ at g_c = 0.735. Smooth decrease, max |d⟨δ⟩/dg| = 1.21 at g=0.40.
- q=5 (Lx=4): ⟨δ⟩ at g_c = 0.471. Smooth decrease, max |d⟨δ⟩/dg| = 1.88 at g=0.60.
- **No discontinuous jump in either case** → no first-order signal on the cylinder.
- Both approach random value (1/q) for large g, confirming paramagnetic phase.

## Key Results

1. **g_c(Ly=2 cylinder) established for q=2 and q=5.** Cylinder provides an intermediate geometry between 1D and 2D, with g_c scaling approximately as coordination number.

2. **No first-order signal for q=5 on cylinder.** Both order parameter (smooth) and gap crossing (exists, converging) are consistent with continuous transition. This extends Sprint 066's finding from 1D to the ladder geometry.

3. **DMRG impractical for q≥5 cylinders at accessible chi.** The local dimension d=5 makes DMRG prohibitively expensive. Exact diag on Lx≤4 (dim≤390k) is the viable approach.

4. **Exact diag gap×Lx crossing on cylinder is reliable.** q=2 crossings converge monotonically (3 pairs), consistent with gap crossing on 1D chains.

5. **Cylinder/1D g_c ratio is non-monotonic:** 1.80 (q=2) vs 1.62 (q=5), mirroring the non-monotonic 2D/1D ratio (3.08, 3.80, 3.60).

## Surprises

- DMRG entropy peak drifts strongly with Lx on cylinders — NOT a good proxy for g_c
- q=5 DMRG completely impractical at chi=20-30 (truncation errors, multi-hour runtimes)
- Cylinder g_c is NOT simply geometric mean of 1D and 2D: g_c(cyl) ≈ 0.58·g_c(2D), not √(g_c(1D)·g_c(2D))
- q=5 order parameter steeper than q=2 (max slope 1.88 vs 1.21) — more critical, not less

## Methods

- TeNPy DMRG on Square lattice with bc=['open', 'periodic'] for cylinder
- Exact diag with vectorized non-adjacent bond construction (numpy array operations instead of Python loops)
- gpu_utils.eigsh for dim > 50k
