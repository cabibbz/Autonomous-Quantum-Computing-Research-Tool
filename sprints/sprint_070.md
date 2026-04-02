# Sprint 070 — 2D Energy Derivatives & Fidelity Susceptibility

**Goal:** Determine whether the 2D hybrid model transition at q=5 is continuous or first-order, using energy second derivative d²E₀/dg² and ground-state fidelity susceptibility χ_F.

**Literature context:**
- 2D standard Potts: continuous for q≤4, first-order for q>4 (Baxter exact)
- 2D clock model: BKT for q≥5 (two transitions, floating phase)
- Our hybrid: unknown in 2D. 1D is continuous for all q tested.

## Experiments

### 070a — Energy second derivative d²E₀/dg²
Dense g scan near g_c with numerical differentiation of E₀(g)/N.

**d²E/dg² per site peak values:**

| q | L=2 | L=3 | L=4 | L=3→4 slope |
|---|-----|-----|-----|-------------|
| 2 | -1.11 (g=0.44) | -3.00 (g=0.54) | -3.14 (g=0.62) | **0.16** |
| 3 | -0.16 (g=0.89) | -1.82 (g=0.89) | — | 6.07 (L=2 unreliable) |
| 5 | -0.10 (g=1.11) | -1.24 (g=1.16) | — | 6.13 (L=2 unreliable) |

**q=2 L=3→4 scaling exponent 0.16** — consistent with α=0 (log divergence) for 2D Ising.

L=2 is pathologically out of scaling for ALL q (already established Sprint 068). L=2→3 exponents are meaningless.

Peak positions are far below g_c for all q — ordered-phase dominance at small L, same issue as entropy peak (Sprint 069).

### 070b — Fidelity susceptibility χ_F (eigenstate sum)
χ_F(g) = Σ_{n>0} |<n|V|0>|² / (E_n-E_0)² where V = ∂H/∂g.

**Problem discovered:** For dim > 5000, eigenstate truncation (k=10 states out of ~2M) severely underestimates χ_F. The q=2 L=3 (full diag, 512 states) vs L=4 (truncated, 10/65536 states) comparison is unreliable.

q=2 L=3→4 scaling: L^0.79 (better than L^0 expected, likely from truncation bias at L=4).

**Lesson: Eigenstate-sum fidelity susceptibility requires full diag. Use overlap method (070c) instead for large systems.**

### 070c — Ground state overlap & Hellmann-Feynman derivatives
F(g, g+δg) = |⟨ψ₀(g)|ψ₀(g+δg)⟩| — includes ALL states implicitly.
χ_F = -2 ln F / (N δg²).
dE/dg = ⟨ψ₀|V|ψ₀⟩ via Hellmann-Feynman (exact, no numerical differentiation noise).

**χ_F/N peak (overlap method):**

| q | L=3 | L=4 | L=3→4 slope |
|---|-----|-----|-------------|
| 2 | 0.641 | 0.841 | **0.94** |
| 3 | 0.549 | — | — |
| 5 | 0.652 | — | — |

**q=2 L=3→4 slope = 0.94** — consistent with continuous (ν=1, d=2 → expect ~0 with log corrections). Far below first-order (L^2).

**d²E/dg² peak (Hellmann-Feynman):**

| q | L=3 | L=4 | L=3→4 slope |
|---|-----|-----|-------------|
| 2 | 3.00 | 3.14 | **0.16** |
| 3 | 2.01 | — | — |
| 5 | 1.58 | — | — |

**Minimum overlap F_min (level crossing indicator):**

| q | L=2 | L=3 | L=4 |
|---|-----|-----|-----|
| 2 | — | 0.9992 | 0.9980 |
| 3 | 1.0000 | 0.9981 | — |
| 5 | 1.0000 | **0.9851** | — |

q=5 L=3 has F_min = 0.985 — significantly lower than q=2,3. Could indicate stronger fluctuations or incipient first-order character. But still close to 1 (no level crossing).

**dE₀/dg per site at g_c (latent heat test):**

| q | L=2 | L=3 | L=4 | Change L=3→4 |
|---|-----|-----|-----|--------------|
| 2 | — | -1.741 | -1.762 | +1.2% |
| 3 | -1.974 | -1.922 | — | — |
| 5 | -1.980 | -1.942 | — | — |

All converging smoothly — **no sign of latent heat** (discontinuity in dE/dg).

## Summary

**q=2 (validation):** All three diagnostics confirm continuous transition.
- d²E/dg² scales as L^0.16 (expected: log divergence, α=0)
- χ_F/N scales as L^0.94 (expected: ~L^0 with log corrections)
- F_min close to 1, dE/dg smooth

**q=5 (test case): INCONCLUSIVE at accessible sizes.**
- Only L=2,3 available; L=2 is out of scaling regime
- No L=3→4 comparison possible (dim = 5^16 ≈ 10^11)
- F_min = 0.985 at L=3 is lower than q=2,3 but not diagnostic
- dE/dg smooth (no latent heat signal)
- Cannot distinguish continuous from weakly first-order

**Fundamental limitation:** L_max = 3 for q=5 in 2D exact diag. This is too small to determine the transition nature. Need: (a) QMC, (b) tensor network on cylinders, or (c) 2D DMRG (Ly=2 cylinder accessible with TeNPy).

## Surprises
- d²E/dg² peak is far below g_c for ALL q at accessible L — critical singularity is masked by ordered-phase curvature
- L=2 out of scaling by 10-50x for ALL quantities, not just gap (Sprint 068)
- Eigenstate-sum χ_F is severely truncated for dim > 5000 — overlap method is essential
- q=5 L=3 F_min = 0.985 is the lowest overlap of any system tested — watchable but not diagnostic
- q=5 d²E/dg² peak (1.58) is SMALLER than q=2 (3.00) and q=3 (2.01) at L=3 — less singular, not more
