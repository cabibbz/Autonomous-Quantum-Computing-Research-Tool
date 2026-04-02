# Sprint 085 — Rényi Entropies Across Walking Boundary

## Hypothesis
Walking breakdown (c_eff deviation from Re(c)) is visible in von Neumann entropy (α=1) but should vanish at α→∞ (probes only λ_max). Different Rényi indices α weight the entanglement spectrum differently:
- S_α = (1/(1-α)) ln(Σ λ_i^α)
- α→0: counts effective dimension
- α=1: von Neumann S = -Σ λ ln(λ) (standard entropy)
- α→∞: S_∞ = -ln(λ_max) (only largest eigenvalue)

CFT prediction for periodic chain (TWO entanglement cuts):
S_α = (c/6)(1 + 1/α) · ln[(N/π)sin(πℓ/N)] + const_α

Extract c_α from size pairs to cancel non-universal additive constant:
c_α = 6 × ΔS_α / ((1+1/α) × ln(N₂/N₁))

## Experiments

### 085a — Rényi entropies S_α for q=2-8 at g_c=1/q
**Status:** Complete. Computed S_α for α = 0.5, 1, 2, 3, 5, 10, ∞ across q=2,3,5,6,7,8 at multiple system sizes (periodic S_q Potts). Initial c_α extraction used wrong prefactor (c/12 instead of c/6 for periodic chain) — corrected in 085c.

Sizes: q=2 n=10-14; q=3 n=8,10; q=5 n=6,8; q=6 n=6,8; q=7 n=6,7; q=8 n=6,7.
Timing: q=6 n=8 (1.68M dim) = 204s, q=8 n=7 (2.1M dim) = 304s.

### 085b — c_α from size pairs (initial, wrong prefactor)
**Status:** Complete. Added smaller sizes (n=4-6) for additional pairs. Used (1+1/α)/12 prefactor → all c_α/Re(c) ≈ 2. Correct periodic formula has (1+1/α)/6 (two entanglement cuts). Corrected in 085c.

### 085c — Corrected analysis: c_α/Re(c) across walking boundary
**Status:** Complete.

**Corrected c_α/Re(c) from widest size pairs:**

| q | pair | α=0.5 | α=1 | α=2 | α=3 | α=5 | α=10 | α=∞ |
|---|------|-------|-----|-----|-----|-----|------|-----|
| 2 | (10,14) | 1.000 | 0.996 | 1.055 | 1.114 | 1.133 | 1.112 | 1.101 |
| 3 | (6,10) | 1.001 | 0.991 | 1.090 | 1.138 | 1.122 | 1.093 | 1.082 |
| 5 | (4,8) | 1.020 | **0.999** | 1.126 | 1.132 | 1.083 | 1.051 | 1.041 |
| 6 | (4,8) | 1.017 | 0.995 | 1.116 | 1.106 | 1.051 | 1.020 | 1.010 |
| 7 | (4,7) | 0.955 | 0.938 | 1.050 | **1.027** | 0.971 | 0.942 | 0.932 |
| 8 | (4,7) | 0.940 | 0.922 | 1.024 | **0.991** | 0.933 | 0.905 | 0.896 |

**Key results:**

1. **α=1 (von Neumann) is closest to Re(c) for walking q≤6.** q=5: 0.06% deviation. q=6: 0.47%.

2. **α=2 always OVERSHOOTS Re(c)** — by 5-13% across all q. Most robust: for q=8, c₂ still within 2.4% of Re(c) while c₁ is off by 7.8%.

3. **α=3 is the "magic" index for walking-broken regime.** For q=8: c₃/Re(c) = 0.991 (0.9% error). For q=7: 1.027 (2.7%). This α uniquely balances the spectral redistribution.

4. **Original hypothesis WRONG: c_∞ does NOT recover Re(c) at high q.** c_∞ deviates MORE than c₁ for q≤6, and comparably for q≥7. The spectral redistribution affects ALL eigenvalues including λ_max.

5. **Rényi spread (c₂ - c_∞)/Re(c) is a new monotonic walking discriminator:**
   - q=2: -0.046, q=3: +0.008, q=5: +0.086, q=6: +0.106, q=7: +0.118, q=8: +0.129
   - Changes sign near q=3, monotonically increases. Walking boundary at zero crossing.

6. **c_α profile shape shifts at walking boundary:**
   - q≤5: peak at α=3-5 (bump in high-α tail)
   - q≥6: peak shifts to α=2 (all higher α decline monotonically)
   - Peak position shifts systematically with walking breakdown

## Surprises
- α=3 (not α=∞) recovers Re(c) in walking-broken regime — nonintuitive
- von Neumann (α=1) is literally the MOST accurate Rényi entropy for q≤6
- α=2 (purity-related) always overshoots by ~5-13% — most biased but most stable
- Rényi spread is monotonic → better walking discriminator than c_eff/Re(c) itself
- Factor-of-2 error nearly led to incorrect physical conclusions (always check periodic vs open BC formula!)

## Methodological note
Periodic chain half-bipartition has TWO entanglement cuts. CFT formula:
- Periodic: S_α = (c/6)(1+1/α) ln(L_eff)  [two cuts]
- Open: S_α = (c/12)(1+1/α) ln(L_eff)  [one cut]
Size-pair differencing cancels additive constants. Widest baseline gives best c_α.

**POTENTIALLY NOVEL:** First systematic Rényi c_α(q,α) mapping across walking boundary for S_q Potts chain. Discovery that α=3 uniquely recovers Re(c) in walking-broken regime. Rényi spread (c₂-c_∞) as new monotonic walking discriminator.
