# Sprint 118 -- q=4 chi_F Extended: Anomalous Exponent Confirmed

**Status:** Complete

## Motivation
Sprint 102c measured q=4 chi_F with only n=6,8, giving alpha=1.69. Exact q=4 Potts predicts alpha=2/nu-1=2.0 (nu=2/3). Extended to 8 sizes (n=4-11) to determine if alpha converges to 2.0 or stabilizes at a different value.

## Experiment 118a -- q=4 chi_F spectral decomposition (n=4 to n=11)

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

Progressive global fits converge: n=4-11 gives 1.786, n=8-11 gives 1.771.

## Analysis

1. **alpha(q=4) = 1.77 is NOT consistent with exact q=4 Potts (alpha=2.0).**
   - Discrepancy is 11.5%, stable across 8 sizes
   - Not a finite-size effect: pairwise alpha is flat for n>=8

2. **Our hybrid model at q=4 is NOT in the standard Potts universality class.**
   - This is consistent with Sprint 065: hybrid != S_q Potts for q>=4
   - The hybrid has Z_4 symmetry, standard Potts has S_4 symmetry
   - Different symmetry -> different critical exponents

3. **Logarithmic formula underpredicts: 1.62 vs measured 1.77 (9% off).**
   - The formula alpha(q) = 1.87*ln(q) - 0.97 was fit to q=5-30 data
   - q=4 is below the walking threshold and may follow different scaling

4. **q=4 is a crossover point.** The effective alpha places q=4 between:
   - Known exact q=3: alpha=1.40 (nu=5/6)
   - Walking q=5: alpha=2.09
   - q=4 alpha=1.77 falls naturally in this progression

5. **z_m = 1.079 at n=10-11.** Converging slowly from 1.097.

6. **Single-multiplet dominance universal at q=4** (frac=1.000 at all sizes).

## Conclusions
- q=4 chi_F gives a clean power-law with alpha=1.77 and 8 sizes
- The hybrid model q=4 exponent is distinct from standard q=4 Potts (alpha=2.0)
- This confirms q=4 is NOT in the standard Potts universality class for our hybrid
- The logarithmic formula for alpha(q) applies for q>=5 (walking regime) only
- alpha(q) can be extended down to q=4 with measured value 1.77
