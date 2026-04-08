"""Sprint 117b: Refit alpha(q) with q=30 data point.

11 data points (q=5-30). q=30 has only n=3,4 (no n=5).
Using pairwise (3,4) alpha=5.384 as the estimate.
Note: this is a lower bound -- global alpha with more sizes would be higher.
"""
import numpy as np
import json, time, os
from scipy.optimize import curve_fit

results = {
    'experiment': '117b_alpha_refit_q30',
    'sprint': 117,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results', 'sprint_117b_alpha_refit_q30.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

# Use pairwise (3,4) for q=30 since no global available
# For consistency, also show what happens using pairwise (3,4) for ALL q
data_global = {
    5:  2.091, 6:  2.377, 7:  2.650, 8:  2.897, 9:  3.159,
    10: 3.349, 12: 3.631, 15: 3.921, 20: 4.590, 25: 5.170,
}
# q=30 only has pair (3,4) = 5.384
data_global[30] = 5.384  # lower bound

data_pair34 = {
    5:  1.956, 6:  2.241, 7:  2.506, 8:  2.750, 9:  3.003,
    10: 3.160, 12: 3.417, 15: 3.815, 20: 4.416, 25: 4.932, 30: 5.384,
}

qs = np.array(sorted(data_global.keys()), dtype=float)
alphas_global = np.array([data_global[int(q)] for q in qs])

print("=" * 70)
print(f"Sprint 117b: alpha(q) refit with q=30 (11 data points)")
print("  NOTE: q=30 uses pair(3,4) alpha -- lower bound on true global")
print("=" * 70)

def logarithmic(q, a, b):
    return a * np.log(q) + b

def log_plus_loglog(q, a, b, c):
    return a * np.log(q) + b * np.log(np.log(q)) + c

def power_law(q, a, b):
    return a * q**b

def sqrt_law(q, a, b):
    return a * np.sqrt(q) + b

def quadratic(q, a, b, c):
    return a * q**2 + b * q + c

def linear(q, a, b):
    return a * q + b

models = {
    'linear':        (linear, 2, ['a', 'b']),
    'quadratic':     (quadratic, 3, ['a', 'b', 'c']),
    'power_law':     (power_law, 2, ['a', 'b']),
    'logarithmic':   (logarithmic, 2, ['a', 'b']),
    'sqrt':          (sqrt_law, 2, ['a', 'b']),
    'log+loglog':    (log_plus_loglog, 3, ['a', 'b', 'c']),
}

n_data = len(qs)
fit_results = {}

print(f"\nModel fits (11 points, q=30 uses pair(3,4)):")
print("-" * 70)

for name, (func, n_params, param_names) in models.items():
    try:
        popt, pcov = curve_fit(func, qs, alphas_global, maxfev=10000)
        predicted = func(qs, *popt)
        residuals = alphas_global - predicted
        rms = np.sqrt(np.mean(residuals**2))
        ss_res = np.sum(residuals**2)
        aic = n_data * np.log(ss_res / n_data) + 2 * n_params

        params_str = ", ".join(f"{pn}={pv:.4f}" for pn, pv in zip(param_names, popt))
        print(f"  {name:14s}: RMS={rms:.4f}, AIC={aic:7.2f}  ({params_str})")

        q50 = func(50.0, *popt)
        q100 = func(100.0, *popt)
        print(f"    -> q=50: {q50:.2f}, q=100: {q100:.2f}")

        fit_results[name] = {
            'params': {pn: float(pv) for pn, pv in zip(param_names, popt)},
            'rms': float(rms),
            'aic': float(aic),
            'residuals': {int(q): float(r) for q, r in zip(qs, residuals)},
        }
    except Exception as e:
        print(f"  {name:14s}: FAILED -- {e}")

# AIC ranking
print(f"\n  AIC ranking:")
sorted_models = sorted(fit_results.items(), key=lambda x: x[1]['aic'])
best_aic = sorted_models[0][1]['aic']
for name, fr in sorted_models:
    delta = fr['aic'] - best_aic
    print(f"    {name:14s}: AIC={fr['aic']:7.2f}, dAIC={delta:+.1f}")

results['fit_results'] = fit_results
results['best_model'] = sorted_models[0][0]
results['qs'] = [int(q) for q in qs]
results['alphas'] = [float(a) for a in alphas_global]
results['note'] = 'q=30 uses pairwise(3,4) alpha -- lower bound on global'

# Leave-one-out cross-validation for top 2 models
print(f"\n{'='*70}")
print("Leave-one-out cross-validation (top models):")
print("-" * 70)
for name in ['logarithmic', 'log+loglog']:
    func, n_params, _ = models[name]
    loo_errors = []
    for i in range(n_data):
        mask = np.ones(n_data, dtype=bool)
        mask[i] = False
        try:
            popt_loo, _ = curve_fit(func, qs[mask], alphas_global[mask], maxfev=10000)
            pred_i = func(qs[i], *popt_loo)
            loo_errors.append(float(alphas_global[i] - pred_i))
        except:
            loo_errors.append(float('nan'))
    loo_rms = np.sqrt(np.nanmean(np.array(loo_errors)**2))
    print(f"  {name:14s}: LOO-RMS={loo_rms:.4f}")
    for i, e in enumerate(loo_errors):
        print(f"    q={int(qs[i]):2d}: error={e:+.4f}")
    results[f'loo_{name}'] = {'rms': float(loo_rms), 'errors': loo_errors}

save()
print("\nResults saved.")
