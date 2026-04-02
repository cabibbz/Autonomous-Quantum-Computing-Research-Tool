# Sprint 070 — 2D Energy Derivatives & Fidelity Susceptibility

**Goal:** Determine whether the 2D hybrid model transition at q=5 is continuous or first-order, using energy second derivative d²E₀/dg² and ground-state fidelity susceptibility χ_F. These scale differently: continuous gives χ_F ~ L^{2/ν}, first-order gives χ_F ~ L^{2d}.

**Literature context:**
- 2D standard Potts: continuous for q≤4, first-order for q>4 (Baxter exact)
- 2D clock model: BKT for q≥5 (two transitions, floating phase)
- Our hybrid: unknown in 2D. 1D is continuous for all q tested.

## Experiments

### 070a — Energy second derivative d²E₀/dg²
- q=2 L=2,3,4 (validation: known continuous, α=0 → log divergence)
- q=3 L=2,3 and q=5 L=2,3
- Dense g scan near g_c, numerical differentiation
- **Results:** [pending]

### 070b — Fidelity susceptibility χ_F
- χ_F(g) = Σ_{n≠0} |<n|∂H/∂g|0>|² / (E_n - E_0)²
- Computed via second-order perturbation theory (no wavefunction overlap needed)
- Same systems as 070a
- **Results:** [pending]

### 070c — Peak scaling analysis
- Compare peak height vs L for q=2 (known) and q=5 (unknown)
- Continuous: peak ~ L^{2/ν - d} for d²E/dg², L^{2/ν} for χ_F
- First-order: peak ~ L^d for both
- **Results:** [pending]
