#!/usr/bin/env python3
"""Sprint 052d: Test g_c = a*sqrt(q-1) + b hypothesis.

The power_offset fit gave exponent 0.516 ≈ 1/2, suggesting g_c ∝ sqrt(q-1).
Test this constrained form and compare to free fit.
Also test g_c = (1/q_c)*ln(q) form (motivated by self-dual q=2,3 structure).
"""
import numpy as np, json
from scipy.optimize import curve_fit

# Corrected data (from exp_052c)
q_data = np.array([2, 3, 4, 5, 7, 10], dtype=float)
gc_data = np.array([0.250, 0.333, 0.3917, 0.4409, 0.5350, 0.6835])
gc_err = np.array([0.001, 0.001, 0.012, 0.013, 0.026, 0.035])

def sqrt_q1(q, a, b):
    return a * np.sqrt(q - 1) + b

def sqrt_q(q, a, b):
    return a * np.sqrt(q) + b

def log_q(q, a, b):
    return a * np.log(q) + b

def power_free(q, a, b, c):
    return a * (q - 1)**b + c

def power_half(q, a, c):
    """Constrained: exponent = 1/2"""
    return a * np.sqrt(q - 1) + c

# Also test: g_c = 1/(2*sqrt(q)) type forms? No, g_c increases.
# Try: g_c = a * q^b (no offset)
def pure_power(q, a, b):
    return a * q**b

# Self-duality motivated: g_c = J/q for q=2,3. What about general?
# q=2: J/q = 0.5 ≠ 0.25. Actually g_c = J/(2q) for q=2 (0.25) and J/q for q=3 (0.333).
# Hmm, not a clean pattern.

models = [
    ('sqrt(q-1)', sqrt_q1, [0.15, 0.1], 2),
    ('sqrt(q)', sqrt_q, [0.15, 0.0], 2),
    ('ln(q)', log_q, [0.2, 0.1], 2),
    ('(q-1)^b free', power_free, [0.2, 0.5, 0.05], 3),
    ('q^b', pure_power, [0.15, 0.65], 2),
]

print("=== Testing functional forms for g_c(q) ===\n")
print(f"Data: q = {q_data.tolist()}")
print(f"g_c  = {gc_data.tolist()}\n")

results = {}
for name, func, p0, npar in models:
    try:
        popt, pcov = curve_fit(func, q_data, gc_data, p0=p0, sigma=gc_err,
                              absolute_sigma=True, maxfev=10000)
        perr = np.sqrt(np.diag(pcov))
        gc_fit = func(q_data, *popt)
        residuals = gc_data - gc_fit
        chi2 = np.sum((residuals / gc_err)**2)
        dof = len(q_data) - npar
        chi2_red = chi2 / dof if dof > 0 else float('inf')
        rms = np.sqrt(np.mean(residuals**2))

        # Predictions
        gc_20 = func(20, *popt)
        gc_50 = func(50, *popt)

        print(f"{name}:")
        print(f"  params = {[f'{p:.5f}' for p in popt]}")
        print(f"  χ²/dof = {chi2_red:.4f}, RMS = {rms:.5f}")
        print(f"  residuals = {[f'{r:.5f}' for r in residuals]}")
        print(f"  g_c(20) = {gc_20:.4f}, g_c(50) = {gc_50:.4f}")

        # Check at known exact points
        for q_test in [2, 3]:
            pred = func(q_test, *popt)
            exact = {2: 0.250, 3: 0.333}[q_test]
            print(f"  q={q_test}: pred={pred:.5f}, exact={exact:.3f}, err={abs(pred-exact)/exact*100:.2f}%")
        print()

        results[name] = {
            'params': popt.tolist(), 'chi2_red': float(chi2_red),
            'rms': float(rms), 'gc_20': float(gc_20), 'gc_50': float(gc_50),
            'residuals': residuals.tolist()
        }
    except Exception as e:
        print(f"{name}: FAILED - {e}\n")

# === Check if the data is consistent with self-duality pattern ===
print("=== Self-duality pattern check ===")
print("If self-duality gives g_c = J/q for q=2,3:")
print(f"  q=2: J/q = 0.5 but actual g_c = 0.25 = J/(2q)")
print(f"  q=3: J/q = 0.333 = actual g_c (J=1)")
print(f"  So q=2 has extra factor of 2 from Z_2 structure (X = X†)")
print()

# === Summary: which simple formula is best? ===
print("=== RANKING ===")
ranked = sorted(results.items(), key=lambda x: x[1]['chi2_red'])
for i, (name, v) in enumerate(ranked):
    print(f"  {i+1}. {name}: χ²/dof={v['chi2_red']:.4f}, g_c(20)={v['gc_20']:.4f}")

# Best simple (2-param) formula
best_2param = [(n, v) for n, v in ranked if len(v['params']) == 2][0]
print(f"\nBest 2-param formula: {best_2param[0]}")
print(f"  g_c(20) = {best_2param[1]['gc_20']:.4f}")

with open('results/sprint_052d_sqrt_test.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved to results/sprint_052d_sqrt_test.json")
