# Sprint 118 -- q=4 chi_F Extended + Model Identity Audit

**Status:** Complete

## Experiment 118a -- q=4 chi_F spectral decomposition (n=4 to n=11)

Model: **S_q Potts** (field = Σ_{k=1}^{q-1} X^k), g_c = 1/q = 0.25.

| n | dim | gap_m | |me|^2 | chi_F | dominant_frac | time |
|---|-----|-------|---------|-------|---------------|------|
| 4 | 256 | 0.9515 | 22.96 | 6.34 | 1.000 | 0.0s |
| 5 | 1,024 | 0.7449 | 26.43 | 9.53 | 1.000 | 0.0s |
| 6 | 4,096 | 0.6100 | 29.50 | 13.21 | 1.000 | 0.0s |
| 7 | 16,384 | 0.5155 | 32.33 | 17.38 | 1.000 | 0.2s |
| 8 | 65,536 | 0.4457 | 35.02 | 22.03 | 1.000 | 0.6s |
| 9 | 262,144 | 0.3923 | 37.59 | 27.14 | 1.000 | 0.9s |
| 10 | 1,048,576 | 0.3500 | 40.07 | 32.71 | 1.000 | 4.9s |
| 11 | 4,194,304 | 0.3158 | 42.48 | 38.73 | 1.000 | 19.4s |

Pairwise alpha convergence:
- (4,5): 1.825
- (5,6): 1.794
- (6,7): 1.780
- (7,8): 1.774
- (8,9): 1.772
- (9,10): 1.771
- (10,11): 1.771

**Alpha stabilizes at 1.771 +/- 0.001 (last 3 pairs).** Not converging to 2.0.

## Analysis (CORRECTED — Model Identity Audit)

### Original interpretation (WRONG):
~~"alpha=1.77 != 2.0 confirms hybrid model != standard Potts universality class."~~

### Corrected interpretation:
The code uses `build_sq_potts` with `for k in range(1, q)` — this IS the standard S_q Potts model, NOT the hybrid. The discrepancy between measured alpha=1.77 and exact alpha=2.0 (from nu=2/3) is most likely due to **multiplicative logarithmic corrections at the marginal q=4 Ashkin-Teller point**.

q=4 is exactly at the boundary between second-order (q<=4) and first-order (q>4) for the S_q Potts model. At this marginal point:
- The q=4 Potts transition has known logarithmic corrections
- The formula alpha = 2/nu - 1 with nu=2/3 gives alpha=2.0, but this assumes pure power-law FSS
- Multiplicative log corrections modify the effective exponent at finite size
- The apparent convergence at alpha=1.77 may be a plateau — convergence to 2.0 could require much larger sizes (n >> 100)
- This is analogous to known behavior at the q=4 BKT/Ashkin-Teller point where exponents converge logarithmically slowly

### What holds:
1. alpha(q=4) = 1.77 is a clean measurement with 8 converged sizes
2. Single-multiplet dominance at q=4 (frac=1.000) — same as q>=5
3. z_m = 1.079 at n=10-11

### What does NOT hold:
1. ~~This does NOT confirm "hybrid != S_q Potts"~~ — the code IS S_q Potts
2. The 11.5% discrepancy with exact alpha=2.0 is likely log corrections, NOT a different universality class

## Model Identity Audit (April 2026)

Code audit of ALL experiments from Sprint 076 onward confirms they use the **standard S_q Potts model**:
- Exact diag: `for k in range(1, q)` (Σ_{k=1}^{q-1} X^k)
- TeNPy DMRG: `SqField = ones(q,q) - eye(q)` (same operator)

Sprints 033-075 used the **Potts-clock hybrid** (X + X†). The switch happened at Sprint 076 and was never reversed. All six claimed novel findings (098, 102-118) are on the standard S_q Potts model — the same model studied by Gorbenko-Rychkov-Zan, Ma & He, Tang et al., and Jacobsen & Wiese.

See KNOWLEDGE.md Model Identity section for full details.
