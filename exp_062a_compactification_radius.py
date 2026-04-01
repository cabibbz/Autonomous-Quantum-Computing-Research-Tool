"""Sprint 062a: Extract effective compactification radius R(q) from x₁(q).

For a compact boson at radius R on a circle:
  - Vertex operators V_{n,m} have dimension x = n²/(4R²) + m²R²/4
  - Lowest: x₁ = min(1/(4R²), R²/4) depending on R vs 1
  - For Z_q orbifold: twist fields have x = k²/(4q²R²) or similar

We test several hypotheses:
  H1: x₁ = 1/(4R²)  → R = 1/(2√x₁)   (momentum mode dominates)
  H2: x₁ = R²/4     → R = 2√x₁         (winding mode dominates)

Then check:
  - Does R(q) explain c(q) via c = 1 + f(R)?
  - Does R² ~ ln(q)?
  - What is c/x₁ as a function of R?

Also: compute x₁ from n=4 descendant gap for all q=2-30 for consistency.
"""
import numpy as np
import json

# Collected x₁ data from Sprints 058, 061
# Best available values (largest n for each q)
x1_data = {
    2:  {"x1": 0.1246, "n_best": 12, "source": "Sprint 058, (10,12) pair"},
    3:  {"x1": 0.1337, "n_best": 10, "source": "Sprint 058, (8,10) pair"},
    4:  {"x1": 0.1172, "n_best": 10, "source": "Sprint 058, (8,10) pair"},
    5:  {"x1": 0.1015, "n_best": 8,  "source": "Sprint 058, (6,8) pair"},
    7:  {"x1": 0.0860, "n_best": 6,  "source": "Sprint 058, (4,6) pair"},
    10: {"x1": 0.0831, "n_best": 6,  "source": "Sprint 058, (4,6) pair"},
    15: {"x1": 0.071,  "n_best": 4,  "source": "Sprint 061, n=4 desc gap"},
    20: {"x1": 0.055,  "n_best": 4,  "source": "Sprint 061, n=4 desc gap"},
    25: {"x1": 0.045,  "n_best": 4,  "source": "Sprint 061, n=4 desc gap"},
    30: {"x1": 0.038,  "n_best": 4,  "source": "Sprint 061, n=4 desc gap"},
}

# Collected c data
c_data = {2: 0.500, 3: 0.800, 4: 1.000, 5: 1.10, 7: 1.30, 10: 1.40}

# c/x1 ratios
cx1_data = {2: 4.013, 3: 5.984, 4: 8.536, 5: 10.840, 7: 15.109, 10: 16.856}

print("=" * 70)
print("COMPACTIFICATION RADIUS ANALYSIS")
print("=" * 70)

# Hypothesis 1: x₁ = 1/(4R²) → R = 1/(2√x₁)
# This is the momentum mode interpretation
print("\n--- Hypothesis H1: x₁ = 1/(4R²) → R = 1/(2√x₁) ---")
print(f"{'q':>3} {'x₁':>8} {'R':>8} {'R²':>8} {'ln(q)':>7} {'R²/ln(q)':>9}")
R_h1 = {}
for q in sorted(x1_data.keys()):
    x1 = x1_data[q]["x1"]
    R = 1.0 / (2.0 * np.sqrt(x1))
    R_h1[q] = R
    lnq = np.log(q)
    ratio = R**2 / lnq if lnq > 0 else float('inf')
    print(f"{q:>3} {x1:>8.4f} {R:>8.4f} {R**2:>8.4f} {lnq:>7.4f} {ratio:>9.4f}")

# Hypothesis 2: x₁ = R²/4 → R = 2√x₁
# This is the winding mode interpretation
print("\n--- Hypothesis H2: x₁ = R²/4 → R = 2√x₁ ---")
print(f"{'q':>3} {'x₁':>8} {'R':>8} {'R²':>8} {'ln(q)':>7} {'R²/ln(q)':>9}")
R_h2 = {}
for q in sorted(x1_data.keys()):
    x1 = x1_data[q]["x1"]
    R = 2.0 * np.sqrt(x1)
    R_h2[q] = R
    lnq = np.log(q)
    ratio = R**2 / lnq if lnq > 0 else float('inf')
    print(f"{q:>3} {x1:>8.4f} {R:>8.4f} {R**2:>8.4f} {lnq:>7.4f} {ratio:>9.4f}")

# Test: does R² scale as ln(q)?
print("\n--- Testing R² ~ a·ln(q) + b ---")
from numpy.polynomial import polynomial as P

qs = np.array(sorted(x1_data.keys()), dtype=float)
lnqs = np.log(qs)

# H1
R2_h1 = np.array([R_h1[q]**2 for q in sorted(x1_data.keys())])
# Linear fit R² = a·ln(q) + b
coeffs_h1 = np.polyfit(lnqs, R2_h1, 1)
print(f"\nH1: R² = {coeffs_h1[0]:.4f}·ln(q) + {coeffs_h1[1]:.4f}")
pred_h1 = np.polyval(coeffs_h1, lnqs)
rms_h1 = np.sqrt(np.mean((R2_h1 - pred_h1)**2))
print(f"    RMS error: {rms_h1:.4f}")

# H2
R2_h2 = np.array([R_h2[q]**2 for q in sorted(x1_data.keys())])
coeffs_h2 = np.polyfit(lnqs, R2_h2, 1)
print(f"\nH2: R² = {coeffs_h2[0]:.4f}·ln(q) + {coeffs_h2[1]:.4f}")
pred_h2 = np.polyval(coeffs_h2, lnqs)
rms_h2 = np.sqrt(np.mean((R2_h2 - pred_h2)**2))
print(f"    RMS error: {rms_h2:.4f}")

# Alternative: x₁ ~ 1/ln(q)?
print("\n--- Testing x₁ ~ a/ln(q) + b ---")
x1_vals = np.array([x1_data[q]["x1"] for q in sorted(x1_data.keys())])
inv_lnq = 1.0 / lnqs
coeffs_x1 = np.polyfit(inv_lnq, x1_vals, 1)
print(f"x₁ = {coeffs_x1[0]:.4f}/ln(q) + {coeffs_x1[1]:.4f}")
pred_x1 = np.polyval(coeffs_x1, inv_lnq)
rms_x1 = np.sqrt(np.mean((x1_vals - pred_x1)**2))
print(f"RMS error: {rms_x1:.5f}")

# Alternative: x₁ ~ a·q^b (power law)?
print("\n--- Testing x₁ ~ a·q^b ---")
log_x1 = np.log(x1_vals)
log_q = np.log(qs)
coeffs_pow = np.polyfit(log_q, log_x1, 1)
b_pow = coeffs_pow[0]
a_pow = np.exp(coeffs_pow[1])
print(f"x₁ ≈ {a_pow:.4f}·q^{b_pow:.4f}")
pred_pow = a_pow * qs**b_pow
rms_pow = np.sqrt(np.mean((x1_vals - pred_pow)**2))
print(f"RMS error: {rms_pow:.5f}")

# Now the key test: does R(q) explain c(q)?
print("\n" + "=" * 70)
print("KEY TEST: Does R(q) explain c(q)?")
print("=" * 70)

# For a compact boson: c = 1 always
# For n free bosons: c = n
# For Z_q orbifold of compact boson: c = 1
# None of these give c ~ ln(q).
#
# Possibility: effective number of bosons n_eff(q) = c(q)
# Then x₁ per boson = x₁(q) (lowest excitation across all sectors)

print("\nIf the CFT has n_eff = c(q) independent boson sectors:")
print(f"{'q':>3} {'c':>6} {'x₁':>7} {'c·x₁':>7} {'c/x₁':>7} {'4R²(H1)':>8}")
for q in sorted(c_data.keys()):
    c = c_data[q]
    x1 = x1_data[q]["x1"]
    R = R_h1[q]
    print(f"{q:>3} {c:>6.3f} {x1:>7.4f} {c*x1:>7.4f} {c/x1:>7.3f} {4*R**2:>8.4f}")

# c·x₁ product
print("\n--- c·x₁ product vs q ---")
cx1_product = {}
for q in sorted(c_data.keys()):
    c = c_data[q]
    x1 = x1_data[q]["x1"]
    cx1_product[q] = c * x1
    print(f"  q={q}: c·x₁ = {c*x1:.5f}")

# If c·x₁ = const, then c = const/x₁ → x₁ determines c fully
# Check if c·x₁ is approximately 1/2q or some other simple function

# Z_q orbifold prediction: x_twist = 1/(4q²R²) for smallest twist field
# x₁ = x_twist when twist fields are lightest
print("\n--- Z_q orbifold test: x₁ = (q-1)/(4q²R²)? ---")
# If x₁ = 1/(4q²R²), then R² = 1/(4q²x₁)
print(f"{'q':>3} {'x₁':>7} {'q²x₁':>8} {'1/(4q²x₁)':>10}")
for q in sorted(x1_data.keys()):
    x1 = x1_data[q]["x1"]
    q2x1 = q**2 * x1
    R_orb = 1.0 / (4 * q**2 * x1)
    print(f"{q:>3} {x1:>7.4f} {q2x1:>8.3f} {R_orb:>10.5f}")

# Harmonic ratio test: if x_k/x_1 = k², then x_k = k² x_1
# Free boson: x_k = k²/(4R²) → x_k/x_1 = k²  ✓
# But we measured sub-quadratic ratios → corrections to free boson

# The self-dual radius is R=1: x_momentum = x_winding = 1/4
# For Ising (q=2): x₁ = 1/8 = 0.125. If x₁ = 1/(4R²), R² = 2 → R=√2
# For q=3 (3-state Potts): x₁ = 2/15 ≈ 0.133. R² = 1/(4·0.133) = 1.880
# The self-dual point R=1 gives x=0.25, which is too high.

print("\n--- Self-dual radius test ---")
print("Self-dual R=1 gives x₁ = 0.25 for compact boson")
print("All measured x₁ < 0.25, so R > 1 (decompactified direction)")
print(f"\n{'q':>3} {'R(H1)':>8} {'R(H1)/sqrt(2)':>14} {'comment':>20}")
for q in sorted(x1_data.keys()):
    R = R_h1[q]
    ratio = R / np.sqrt(2)
    comment = ""
    if q == 2:
        comment = f"Ising R=√2={np.sqrt(2):.3f}"
    elif q == 3:
        comment = f"3-Potts CFT"
    print(f"{q:>3} {R:>8.4f} {ratio:>14.4f} {comment:>20}")

# Key: does R grow without bound?
# If x₁ → 0, then R → ∞ under H1 (decompactification)
# This is consistent with c → ∞ if the boson decompactifies

# For a non-compact boson (R→∞): c=1, x₁→0
# But we see c>1 AND x₁→0. So it's NOT a single decompactifying boson.

print("\n" + "=" * 70)
print("CONCLUSION: Single compact boson FAILS")
print("=" * 70)
print("""
The data shows:
1. x₁ → 0 as q → ∞  (consistent with R → ∞, decompactification)
2. c → ∞ as q → ∞  (INCONSISTENT with single boson, always c=1)
3. Harmonic ratios → k²  (consistent with free boson)

This combination REQUIRES multiple effective degrees of freedom.
Possible interpretations:
  A. c(q) independent Luttinger liquids (one per Z_q generator?)
  B. A single boson with c=1 PLUS additional massive modes going critical
  C. An entirely new universality class

Test B via twisted boundary conditions in next experiment.
""")

# Fit: x₁ = a/(ln(q-1))^b
print("--- Refined fit: x₁ = a·(ln(q-1))^(-b) ---")
qs_fit = np.array([q for q in sorted(x1_data.keys()) if q >= 3])
x1_fit = np.array([x1_data[q]["x1"] for q in qs_fit])
ln_qm1 = np.log(qs_fit - 1)
log_x1_fit = np.log(x1_fit)

coeffs_ref = np.polyfit(np.log(ln_qm1), log_x1_fit, 1)
b_ref = -coeffs_ref[0]
a_ref = np.exp(coeffs_ref[1])
print(f"x₁ ≈ {a_ref:.4f}·(ln(q-1))^(-{b_ref:.4f})")
pred_ref = a_ref * ln_qm1**(-b_ref)
rms_ref = np.sqrt(np.mean((x1_fit - pred_ref)**2))
print(f"RMS error: {rms_ref:.5f}")

for q in qs_fit:
    x1 = x1_data[q]["x1"]
    pred = a_ref * np.log(q-1)**(-b_ref)
    err = abs(pred - x1) / x1 * 100
    print(f"  q={q}: x₁={x1:.4f}, predicted={pred:.4f}, error={err:.1f}%")

# Save results
results = {
    "x1_data": {str(q): v for q, v in x1_data.items()},
    "R_momentum_mode": {str(q): float(R_h1[q]) for q in sorted(x1_data.keys())},
    "R_winding_mode": {str(q): float(R_h2[q]) for q in sorted(x1_data.keys())},
    "fit_R2_vs_lnq_H1": {"a": float(coeffs_h1[0]), "b": float(coeffs_h1[1]), "rms": float(rms_h1)},
    "fit_x1_vs_inv_lnq": {"a": float(coeffs_x1[0]), "b": float(coeffs_x1[1]), "rms": float(rms_x1)},
    "fit_x1_power_law": {"a": float(a_pow), "b": float(b_pow), "rms": float(rms_pow)},
    "fit_x1_ln_power": {"a": float(a_ref), "b": float(b_ref), "rms": float(rms_ref)},
    "conclusion": "Single compact boson FAILS: x₁→0 and c>1 incompatible with c=1. Multiple DOF required.",
}
with open("results/sprint_062a_compactification.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_062a_compactification.json")
