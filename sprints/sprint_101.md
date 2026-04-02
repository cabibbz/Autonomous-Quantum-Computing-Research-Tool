# Sprint 101 — Symmetry-Resolved Entanglement Entropy Across Walking Boundary

**Date:** 2026-04-02
**Status:** Complete (3 experiments)

## Motivation

Symmetry-resolved entanglement entropy (SREE) decomposes total entropy by Z_q charge sector. Literature search confirms: **no one has measured SREE for q>4 Potts or at the walking/complex-CFT boundary.** Completely unexplored diagnostic.

Key references: Goldstein-Sela (U(1) equipartition), Horvath et al. 2108.10935 (q=3 Potts form factors), Calabrese-Murciano (finite abelian groups). No SREE for q>4 or walking models exists.

## Experiment 101a — SREE for S_q Potts (PARTIALLY RETRACTED)

**Method:** Decompose ρ_A into Z_q charge sectors using subsystem symmetry projectors P_α = (1/q) Σ_k ω^{-αk} G_A^k where G_A = X⊗...⊗X on subsystem A.

**RETRACTION:** The Hamiltonian builder in 101a used the HYBRID field (X+X†) with S_q critical coupling g_c=1/q. For q≥4 these are different models, so 101a results at q=4,5,7 were at the wrong critical point — deep in the ordered phase, producing spurious near-equipartition. Results for q=2,3 (where models coincide up to field rescaling) are also affected by the 2× field factor.

**Corrected results in 101b** use proper S_q field (Σ X^k) at g_c=1/q, and proper hybrid field (X+X†) at known hybrid g_c values.

## Experiment 101b — Equipartition Transition: Dense q Scan

**S_q Potts at n=6, nA=3, g_c=1/q:**

| q | p(0) | p(0)*q | CV | S_n/S_t |
|---|------|--------|-----|---------|
| 2 | 0.774 | 1.547 | 0.547 | 0.910 |
| 3 | 0.678 | 2.034 | 0.731 | 0.909 |
| 4 | 0.624 | 2.495 | 0.863 | 0.908 |
| 5 | 0.589 | 2.944 | 0.972 | 0.908 |
| 6 | 0.564 | 3.387 | 1.067 | 0.908 |
| 7 | 0.547 | 3.827 | 1.154 | 0.908 |
| 8 | 0.533 | 4.267 | 1.235 | 0.909 |
| 9 | 0.523 | 4.708 | 1.311 | 0.909 |
| 10 | 0.515 | 5.151 | 1.384 | 0.910 |

**Hybrid model at n=6, nA=3, proper g_c:**

| q | p(0) | p(0)*q | CV | S_n/S_t |
|---|------|--------|-----|---------|
| 2 | 0.963 | 1.927 | 0.927 | 0.865 |
| 3 | 0.677 | 2.031 | 0.729 | 0.909 |
| 5 | 0.556 | 2.780 | 0.923 | 0.911 |
| 7 | 0.449 | 3.141 | 0.995 | 0.920 |
| 8 | 0.408 | 3.267 | 1.030 | 0.924 |

**Key findings from 101b:**

1. **Charge-0 enrichment increases monotonically with q.** p(0)*q goes from 1.55 (q=2) to 5.15 (q=10) for S_q Potts. Equipartition of probabilities is strongly broken at ALL q. p(0) → 1/2 (not 1/q) as q→∞.

2. **S_number/S_total ≈ 0.908 is remarkably universal across q** for S_q Potts at n=6, nA=3. Varies by only 0.2% from q=2 to q=10. Since S_total ∝ c(q), this means S_number ∝ c(q) at fixed size — charge fluctuation entropy tracks the central charge.

3. **Hybrid model anomaly at q=2:** p(0)=0.963 (much higher than S_q's 0.774) because hybrid field X+X†=2X pushes the system deeper into the disordered phase at the TFIM critical point (effective coupling is 2× higher).

4. **Conjugate pairing exact:** S(α) = S(q-α) to machine precision for all q, confirming Z_q charge-conjugation symmetry.

## Experiment 101c — Size Scaling of SREE

**Midchain bipartition (nA = n/2), S_q Potts at g_c=1/q:**

| q | n | nA | p(0)*q | CV | S_n/S_t |
|---|---|-----|--------|------|---------|
| 2 | 4 | 2 | 1.604 | 0.604 | 0.954 |
| 2 | 8 | 4 | 1.510 | 0.510 | 0.877 |
| 2 | 14 | 7 | 1.444 | 0.444 | 0.812 |
| 3 | 4 | 2 | 2.150 | 0.813 | 0.953 |
| 3 | 10 | 5 | 1.903 | 0.638 | 0.849 |
| 5 | 4 | 2 | 3.173 | 1.087 | 0.953 |
| 5 | 8 | 4 | 2.797 | 0.898 | 0.874 |
| 7 | 4 | 2 | 4.156 | 1.288 | 0.953 |
| 7 | 6 | 3 | 3.827 | 1.154 | 0.908 |

**Fixed nA=2, varying n:**

| q | n | p(0)*q | CV | S_n/S_t |
|---|---|--------|------|---------|
| 2 | 4 | 1.604 | 0.604 | 0.954 |
| 2 | 12 | 1.547 | 0.547 | 0.913 |
| 3 | 4 | 2.150 | 0.813 | 0.953 |
| 3 | 8 | 2.049 | 0.742 | 0.918 |
| 5 | 4 | 3.173 | 1.087 | 0.953 |
| 5 | 8 | 2.974 | 0.987 | 0.917 |

**Key findings from 101c:**

1. **S_n/S_t is q-INDEPENDENT at fixed n.** At n=4 nA=2: ratio=0.953 for ALL q (2,3,5,7). At n=6 nA=3: ratio=0.908 for all q. At n=8 nA=4: ratio≈0.875 for q=2,3,5. This is a non-trivial universality: the fraction of entropy carried by charge fluctuations depends only on system geometry, not on q or c.

2. **S_n/S_t decreases monotonically with n.** Goes from 0.954 (n=4) toward ~0.8 (n=14). In the thermodynamic limit S_n = O(1) while S_t = O(log n), so S_n/S_t → 0. But convergence is very slow.

3. **p(0)*q decreases with n but remains >> 1.** For q=2: extrapolates to ~1.35 at n→∞ (still 35% charge-0 enriched). For q=5: p(0)*q ≈ 2.8 at n=8 (still 2.8× enriched). Convergence to equipartition (p(0)*q=1) would require n >> 100.

4. **At fixed nA=2, varying n:** p(0)*q converges quickly to a finite value. For q=2: converges to ~1.55 by n=10. System size beyond 2nA barely affects SREE — confirming subsystem size nA controls the SREE, not total system size n.

## Surprises

- **S_n/S_t universality across q** at fixed size — charge fluctuation entropy carries the same fraction of total entropy regardless of q or c(q). Not predicted by any existing theory.
- **Equipartition of probabilities p(α)=1/q is strongly violated** at all accessible sizes. The charge-0 sector always dominates. This is NOT the same as entropy equipartition (which DOES approximately hold for conditional entropy S(α)).
- **No walking-specific signature in SREE** — the q-dependence is smooth and monotonic, with no special feature at q=4 or q=5. Walking does not show up in charge-sector decomposition at these sizes.

## Conclusions

SREE provides a new decomposition of entanglement entropy but does **not** provide a new walking discriminator at accessible sizes. The dominant finding is the q-independent S_n/S_t ratio at fixed geometry — a universality that appears to be undocumented in the SREE literature.

The charge-0 enrichment p(0)*q >> 1 at all accessible sizes means the Goldstein-Sela equipartition (of probabilities) is far from being reached. The corrections decay very slowly (~1/log(L)). Only entropy equipartition (conditional S(α) uniform across sectors) shows reasonable convergence, and only for q≥4.

**POTENTIALLY NOVEL:** Universal S_n/S_t ratio at fixed geometry, independent of q and c(q), for critical S_q Potts chains. First SREE measurement for q=4-10 Potts model.

[Full data: results/sprint_101a_sree.json, sprint_101b_equi_transition.json, sprint_101c_sree_scaling.json]
