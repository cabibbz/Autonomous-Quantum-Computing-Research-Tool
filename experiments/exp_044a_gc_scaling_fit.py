#!/usr/bin/env python3
"""Sprint 044a: Fit g_c vs q to candidate scaling laws.

Data points (Potts model, n=8,12 MI-CV crossings):
  q=2: g_c=1.0 (TFIM exact)
  q=3: g_c=1.0 (Potts self-dual exact)
  q=4: g_c≈0.893 (Sprint 040)
  q=5: g_c≈0.41 (Sprint 042)
  q=10: g_c≈0.246 (Sprint 043)

Candidate forms:
  1. Power law: g_c = A / q^alpha
  2. Inverse: g_c = A / (q - q0)
  3. Exponential: g_c = A * exp(-B * q)
  4. Logarithmic: g_c = A / log(q)^beta
  5. Rational: g_c = A / (1 + B*(q-2)^C)
"""
import numpy as np
from scipy.optimize import curve_fit, minimize
import json

# Data
q_data = np.array([2, 3, 4, 5, 10], dtype=float)
gc_data = np.array([1.0, 1.0, 0.893, 0.41, 0.246])
gc_err = np.array([0.01, 0.01, 0.03, 0.02, 0.02])  # rough uncertainties from finite-size shift

# Candidate functional forms
def power_law(q, A, alpha):
    return A / q**alpha

def inverse_shifted(q, A, q0):
    return A / (q - q0)

def exponential(q, A, B):
    return A * np.exp(-B * q)

def log_law(q, A, beta):
    return A / np.log(q)**beta

def rational(q, A, B, C):
    return A / (1 + B * (q - 2)**C)

# Also test the known 1D quantum Potts self-dual point formula:
# For q=2,3 (Ising/3-state Potts), g_c = 1/√(q-1) at self-dual point
# Let's check: q=2 -> 1/1=1 ✓, q=3 -> 1/√2=0.707 ✗ (actual is 1.0)
# Actually the self-dual point for 1D quantum q-state Potts is g_c = 1 for q=2,3
# and the formula differs. Let's try g_c = (q-1)^(-gamma)
def potts_power(q, gamma):
    return (q - 1)**(-gamma)

# Fit each model
results = {}

# 1. Power law: g_c = A / q^alpha
try:
    popt, pcov = curve_fit(power_law, q_data, gc_data, p0=[2.0, 0.5], sigma=gc_err)
    pred = power_law(q_data, *popt)
    resid = np.sum(((gc_data - pred) / gc_err)**2)
    results['power_law'] = {'params': {'A': popt[0], 'alpha': popt[1]},
                            'chi2': resid, 'formula': f'g_c = {popt[0]:.3f} / q^{popt[1]:.3f}',
                            'predictions': {int(q): float(p) for q, p in zip(q_data, pred)}}
    print(f"Power law: g_c = {popt[0]:.4f} / q^{popt[1]:.4f}, χ²={resid:.3f}")
except Exception as e:
    print(f"Power law fit failed: {e}")

# 2. Inverse shifted: g_c = A / (q - q0)
try:
    popt, pcov = curve_fit(inverse_shifted, q_data, gc_data, p0=[3.0, -1.0], sigma=gc_err)
    pred = inverse_shifted(q_data, *popt)
    resid = np.sum(((gc_data - pred) / gc_err)**2)
    results['inverse_shifted'] = {'params': {'A': popt[0], 'q0': popt[1]},
                                   'chi2': resid, 'formula': f'g_c = {popt[0]:.3f} / (q - {popt[1]:.3f})',
                                   'predictions': {int(q): float(p) for q, p in zip(q_data, pred)}}
    print(f"Inverse shifted: g_c = {popt[0]:.4f} / (q - {popt[1]:.4f}), χ²={resid:.3f}")
except Exception as e:
    print(f"Inverse shifted fit failed: {e}")

# 3. Exponential: g_c = A * exp(-B*q)
try:
    popt, pcov = curve_fit(exponential, q_data, gc_data, p0=[2.0, 0.2], sigma=gc_err)
    pred = exponential(q_data, *popt)
    resid = np.sum(((gc_data - pred) / gc_err)**2)
    results['exponential'] = {'params': {'A': popt[0], 'B': popt[1]},
                               'chi2': resid, 'formula': f'g_c = {popt[0]:.3f} * exp(-{popt[1]:.3f}*q)',
                               'predictions': {int(q): float(p) for q, p in zip(q_data, pred)}}
    print(f"Exponential: g_c = {popt[0]:.4f} * exp(-{popt[1]:.4f}*q), χ²={resid:.3f}")
except Exception as e:
    print(f"Exponential fit failed: {e}")

# 4. Log law: g_c = A / log(q)^beta
try:
    popt, pcov = curve_fit(log_law, q_data, gc_data, p0=[1.0, 2.0], sigma=gc_err)
    pred = log_law(q_data, *popt)
    resid = np.sum(((gc_data - pred) / gc_err)**2)
    results['log_law'] = {'params': {'A': popt[0], 'beta': popt[1]},
                           'chi2': resid, 'formula': f'g_c = {popt[0]:.3f} / log(q)^{popt[1]:.3f}',
                           'predictions': {int(q): float(p) for q, p in zip(q_data, pred)}}
    print(f"Log law: g_c = {popt[0]:.4f} / log(q)^{popt[1]:.4f}, χ²={resid:.3f}")
except Exception as e:
    print(f"Log law fit failed: {e}")

# 5. (q-1)^(-gamma) power law
try:
    popt, pcov = curve_fit(potts_power, q_data, gc_data, p0=[0.5], sigma=gc_err)
    pred = potts_power(q_data, *popt)
    resid = np.sum(((gc_data - pred) / gc_err)**2)
    results['potts_power'] = {'params': {'gamma': popt[0]},
                               'chi2': resid, 'formula': f'g_c = (q-1)^(-{popt[0]:.3f})',
                               'predictions': {int(q): float(p) for q, p in zip(q_data, pred)}}
    print(f"Potts power: g_c = (q-1)^(-{popt[0]:.4f}), χ²={resid:.3f}")
except Exception as e:
    print(f"Potts power fit failed: {e}")

# 6. Rational: g_c = A / (1 + B*(q-2)^C)  (3 params)
try:
    popt, pcov = curve_fit(rational, q_data, gc_data, p0=[1.0, 0.5, 1.5], sigma=gc_err,
                           maxfev=5000)
    pred = rational(q_data, *popt)
    resid = np.sum(((gc_data - pred) / gc_err)**2)
    # AIC-like: penalize extra params
    n_pts = len(q_data)
    aic = resid + 2 * 3  # 3 params
    results['rational'] = {'params': {'A': popt[0], 'B': popt[1], 'C': popt[2]},
                            'chi2': resid, 'aic': aic,
                            'formula': f'g_c = {popt[0]:.3f} / (1 + {popt[1]:.3f}*(q-2)^{popt[2]:.3f})',
                            'predictions': {int(q): float(p) for q, p in zip(q_data, pred)}}
    print(f"Rational: g_c = {popt[0]:.4f} / (1 + {popt[1]:.4f}*(q-2)^{popt[2]:.4f}), χ²={resid:.3f}")
except Exception as e:
    print(f"Rational fit failed: {e}")

# Rank by chi-squared
print("\n=== RANKING (by χ²) ===")
ranked = sorted(results.items(), key=lambda x: x[1]['chi2'])
for i, (name, r) in enumerate(ranked):
    n_params = len(r['params'])
    # BIC-like penalty
    bic = r['chi2'] + n_params * np.log(len(q_data))
    print(f"  {i+1}. {name}: χ²={r['chi2']:.3f}, n_params={n_params}, BIC={bic:.3f}")
    print(f"     {r['formula']}")
    print(f"     Predictions: {r['predictions']}")

# Predictions for untested q values
print("\n=== PREDICTIONS for untested q ===")
q_pred = [6, 7, 8, 15, 20, 50, 100]
best_name, best = ranked[0]
print(f"Best model: {best_name}")
print(f"Formula: {best['formula']}")
for q in q_pred:
    # Evaluate each model
    preds = {}
    for name, r in results.items():
        p = r['params']
        if name == 'power_law':
            val = power_law(q, p['A'], p['alpha'])
        elif name == 'inverse_shifted':
            val = inverse_shifted(q, p['A'], p['q0'])
        elif name == 'exponential':
            val = exponential(q, p['A'], p['B'])
        elif name == 'log_law':
            val = log_law(q, p['A'], p['beta'])
        elif name == 'potts_power':
            val = potts_power(q, p['gamma'])
        elif name == 'rational':
            val = rational(q, p['A'], p['B'], p['C'])
        preds[name] = float(val)
    vals = list(preds.values())
    spread = max(vals) - min(vals)
    print(f"  q={q:3d}: best={preds[best_name]:.4f}, range=[{min(vals):.4f}, {max(vals):.4f}], spread={spread:.4f}")

# Save results
with open('results/sprint_044a_gc_scaling.json', 'w') as f:
    json.dump({
        'sprint': '044a',
        'description': 'g_c vs q scaling law fits',
        'data': {'q': q_data.tolist(), 'gc': gc_data.tolist(), 'gc_err': gc_err.tolist()},
        'fits': {k: {kk: vv for kk, vv in v.items()} for k, v in results.items()},
        'ranking': [(name, r['chi2']) for name, r in ranked],
        'predictions': {q: {name: float(
            power_law(q, results[name]['params']['A'], results[name]['params']['alpha']) if name == 'power_law' else
            potts_power(q, results[name]['params']['gamma']) if name == 'potts_power' else 0
        ) for name in ['power_law', 'potts_power'] if name in results} for q in q_pred}
    }, f, indent=2, default=str)

print("\nResults saved to results/sprint_044a_gc_scaling.json")
