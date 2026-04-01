#!/usr/bin/env python3
"""Sprint 052a: Fit g_c(q) scaling law from corrected data points.

Data: q=2,3,4,5,7 with g_c = 0.250, 0.333, 0.392, 0.441, 0.524
Candidate forms:
  1. Power law: g_c = a * q^b
  2. Logarithmic: g_c = a * ln(q) + b
  3. Square root: g_c = a * sqrt(q) + b
  4. Rational: g_c = a * (q-1)/q + b  (asymptotes to a+b)
  5. Inverse: g_c = a - b/q  (asymptotes to a)
  6. Power with offset: g_c = a * (q-1)^b + c
"""
import numpy as np, json
from scipy.optimize import curve_fit

# Data
q_data = np.array([2, 3, 4, 5, 7], dtype=float)
gc_data = np.array([0.250, 0.333, 0.392, 0.441, 0.524])
# Uncertainties: exact for q=2,3 (~0), ~3% for q=4,5,7
gc_err = np.array([0.001, 0.001, 0.012, 0.013, 0.016])

# Candidate forms
def power_law(q, a, b):
    return a * q**b

def logarithmic(q, a, b):
    return a * np.log(q) + b

def sqrt_form(q, a, b):
    return a * np.sqrt(q) + b

def rational(q, a, b):
    return a * (q - 1) / q + b

def inverse(q, a, b):
    return a - b / q

def power_offset(q, a, b, c):
    return a * (q - 1)**b + c

models = {
    'power_law': (power_law, [0.15, 0.5], 2),
    'logarithmic': (logarithmic, [0.2, 0.1], 2),
    'sqrt': (sqrt_form, [0.15, 0.05], 2),
    'rational': (rational, [0.6, -0.05], 2),
    'inverse': (inverse, [0.7, 0.9], 2),
    'power_offset': (power_offset, [0.15, 0.5, 0.1], 3),
}

results = {'data': {'q': q_data.tolist(), 'gc': gc_data.tolist(), 'gc_err': gc_err.tolist()}}
fits = {}

print("=== Fitting g_c(q) to candidate forms ===\n")

for name, (func, p0, nparams) in models.items():
    try:
        popt, pcov = curve_fit(func, q_data, gc_data, p0=p0, sigma=gc_err, absolute_sigma=True, maxfev=10000)
        perr = np.sqrt(np.diag(pcov))
        gc_fit = func(q_data, *popt)
        residuals = gc_data - gc_fit
        chi2 = np.sum((residuals / gc_err)**2)
        dof = len(q_data) - nparams
        chi2_red = chi2 / dof if dof > 0 else float('inf')
        rms = np.sqrt(np.mean(residuals**2))

        # Predictions
        gc_10 = func(10, *popt)
        gc_20 = func(20, *popt)
        gc_100 = func(100, *popt)

        fit_result = {
            'params': popt.tolist(),
            'param_err': perr.tolist(),
            'chi2': float(chi2),
            'chi2_red': float(chi2_red),
            'rms': float(rms),
            'residuals': residuals.tolist(),
            'gc_10': float(gc_10),
            'gc_20': float(gc_20),
            'gc_100': float(gc_100),
        }
        fits[name] = fit_result

        print(f"{name}:")
        print(f"  params = {popt}")
        print(f"  χ²/dof = {chi2_red:.4f}  (RMS = {rms:.5f})")
        print(f"  residuals = {residuals}")
        print(f"  g_c(10) = {gc_10:.4f}, g_c(20) = {gc_20:.4f}, g_c(100) = {gc_100:.4f}")
        print()
    except Exception as e:
        print(f"{name}: FAILED - {e}\n")
        fits[name] = {'error': str(e)}

# Rank by chi2_red
print("\n=== Ranking (by χ²/dof) ===")
ranked = sorted([(k, v) for k, v in fits.items() if 'chi2_red' in v],
                key=lambda x: x[1]['chi2_red'])
for i, (name, v) in enumerate(ranked):
    print(f"  {i+1}. {name}: χ²/dof = {v['chi2_red']:.4f}, g_c(10) = {v['gc_10']:.4f}")

best_name = ranked[0][0]
best = ranked[0][1]
print(f"\nBest fit: {best_name}")
print(f"  g_c(10) prediction = {best['gc_10']:.4f}")
print(f"  g_c(20) prediction = {best['gc_20']:.4f}")

# Spread of predictions at q=10 across all good fits
gc10_all = [v['gc_10'] for _, v in ranked if v['chi2_red'] < 10]
print(f"\n  g_c(10) range across good fits: [{min(gc10_all):.4f}, {max(gc10_all):.4f}]")
print(f"  g_c(10) mean ± std: {np.mean(gc10_all):.4f} ± {np.std(gc10_all):.4f}")

# Check: does g_c → ∞ or saturate?
print(f"\n=== Asymptotic behavior ===")
for name, v in ranked:
    print(f"  {name}: g_c(100) = {v['gc_100']:.4f}")

results['fits'] = fits
results['ranking'] = [name for name, _ in ranked]
results['best_fit'] = best_name
results['gc10_prediction'] = {'mean': float(np.mean(gc10_all)), 'std': float(np.std(gc10_all)),
                               'range': [float(min(gc10_all)), float(max(gc10_all))]}

with open('results/sprint_052a_gc_scaling.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nSaved to results/sprint_052a_gc_scaling.json")
