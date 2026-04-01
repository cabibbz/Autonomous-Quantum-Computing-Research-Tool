# Sprint 060 — OPE Coefficients: Structure Constants of the Novel q>4 Potts CFT

## Idea
We've confirmed genuine CFT with conformal towers (Sprint 059), measured c, x_1, operator content, and degeneracy patterns for q=2-10. The next step to fully characterize this CFT is to extract **OPE (Operator Product Expansion) coefficients** — the structure constants C_{ijk} that determine all correlation functions.

**Method:** On a periodic chain of N sites, the matrix element of a local operator O between CFT eigenstates gives:
- C_{I,sigma,sigma} proportional to N^{x_sigma} * <sigma|O_local(0)|0>
- C_{sigma,sigma,epsilon} proportional to N^{x_epsilon} * <epsilon|O_local(0)|sigma>

We use the Potts clock operator Z = diag(1, omega, omega^2, ..., omega^{q-1}) where omega = e^{2*pi*i/q} as the local spin field operator. Validation on q=2 (Ising, exact C_{sigma,sigma,epsilon} = 1/2) and q=3 (known from W_3 minimal model).

**Literature context:** Tang et al. (PRL 2024) measured OPE coefficients for non-Hermitian q=5. We measure them for the Hermitian model at q=2-10.

## Experiments

### Exp 060a — OPE validation on q=2 Ising
*Status: COMPLETE*

**Method validated.** The ratio C_{sigma,sigma,epsilon} = |<epsilon|Z_0|sigma>| / |<0|Z_0|sigma>| converges to exact 0.5:

| n | |<0|Z|sigma>| | |<eps|Z|sigma>| | C_ratio | error |
|---|---------------|-----------------|---------|-------|
| 6 | 0.8075 | 0.3967 | 0.4913 | 1.7% |
| 8 | 0.7791 | 0.3858 | 0.4951 | 1.0% |
| 10 | 0.7577 | 0.3765 | 0.4969 | 0.6% |
| 12 | 0.7407 | 0.3687 | 0.4979 | 0.4% |

Scaling check: |<0|Z|sigma>| * N^{x_sigma} = 1.0101-1.0104 (constant, confirming x_sigma=1/8).
Convergence is monotonic from below, error ~ 1/N^{1.5}.

### Exp 060b — OPE for q=3 Potts
*Status: pending*

### Exp 060c — OPE for q=4,5 (new territory)
*Status: pending*

### Exp 060d — OPE for q=7,10 (novel CFT regime)
*Status: pending*
