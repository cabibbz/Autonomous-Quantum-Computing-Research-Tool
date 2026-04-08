"""Sprint 116b: Refit alpha(q) with q=25 data point.

10 data points (q=5-25). q=25 measured in exp_116a: alpha_global=5.170.
Logarithmic predicted 4.94 — off by +0.23 (underpredicts).
"""
import numpy as np
import json, time, os
from scipy.optimize import curve_fit

results = {
    'experiment': '116b_alpha_refit_q25',
    'sprint': 116,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results', 'sprint_116b_alpha_refit_q25.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

# Collected alpha data from sprints 103, 110-116
data = {
    5:  {'alpha_global': 2.091, 'alpha_last_pair': 2.100, 'z_m': 1.028, 'beta_me': 0.125},
    6:  {'alpha_global': 2.377, 'alpha_last_pair': 2.402, 'z_m': 1.072, 'beta_me': 0.260},
    7:  {'alpha_global': 2.650, 'alpha_last_pair': 2.670, 'z_m': 1.238, 'beta_me': 0.418},
    8:  {'alpha_global': 2.897, 'alpha_last_pair': 2.897, 'z_m': 1.320, 'beta_me': 0.564},
    9:  {'alpha_global': 3.159, 'alpha_last_pair': 3.159, 'z_m': 1.438, 'beta_me': 0.709},
    10: {'alpha_global': 3.349, 'alpha_last_pair': 3.417, 'z_m': 1.566, 'beta_me': 1.285},
    12: {'alpha_global': 3.631, 'alpha_last_pair': 3.734, 'z_m': 1.662, 'beta_me': 1.307},
    15: {'alpha_global': 3.921, 'alpha_last_pair': 4.070, 'z_m': 1.784, 'beta_me': 1.353},
    20: {'alpha_global': 4.590, 'alpha_last_pair': 4.835, 'z_m': 2.047, 'beta_me': 1.496},
    25: {'alpha_global': 5.170, 'alpha_last_pair': 5.504, 'z_m': 2.286, 'beta_me': 1.598},
}

qs = np.array(sorted(data.keys()), dtype=float)
alphas_global = np.array([data[int(q)]['alpha_global'] for q in qs])
alphas_last = np.array([data[int(q)]['alpha_last_pair'] for q in qs])
z_ms = np.array([data[int(q)]['z_m'] for q in qs])
beta_mes = np.array([data[int(q)]['beta_me'] for q in qs])

print("=" * 70)
print("Sprint 116b: alpha(q) refit with q=25 (10 data points)")
print("=" * 70)

print("\nData (global alpha / last-pair alpha):")
for q_val in qs:
    d = data[int(q_val)]
    print(f"  q={int(q_val):2d}: alpha_global={d['alpha_global']:.3f}, "
          f"alpha_last={d['alpha_last_pair']:.3f}, z_m={d['z_m']:.3f}, beta_me={d['beta_me']:.3f}")

def linear(q, a, b):
    return a * q + b

def quadratic(q, a, b, c):
    return a * q**2 + b * q + c

def power_law(q, a, b):
    return a * q**b

def logarithmic(q, a, b):
    return a * np.log(q) + b

def sqrt_law(q, a, b):
    return a * np.sqrt(q) + b

def log_plus_loglog(q, a, b, c):
    """a*ln(q) + b*ln(ln(q)) + c — tests subleading log correction."""
    return a * np.log(q) + b * np.log(np.log(q)) + c

# Fit both global and last-pair
for alphas, label in [(alphas_global, "global"), (alphas_last, "last-pair")]:
    print(f"\n{'='*70}")
    print(f"Model fits ({label} alpha, {len(qs)} points):")
    print("-" * 70)

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
    for name, (func, n_params, param_names) in models.items():
        try:
            popt, pcov = curve_fit(func, qs, alphas, maxfev=10000)
            predicted = func(qs, *popt)
            residuals = alphas - predicted
            rms = np.sqrt(np.mean(residuals**2))
            ss_res = np.sum(residuals**2)
            aic = n_data * np.log(ss_res / n_data) + 2 * n_params

            params_str = ", ".join(f"{pn}={pv:.4f}" for pn, pv in zip(param_names, popt))
            print(f"  {name:14s}: RMS={rms:.4f}, AIC={aic:7.2f}  ({params_str})")

            # Predictions at q=30, 50
            q30 = func(30.0, *popt)
            q50 = func(50.0, *popt)
            print(f"    -> q=30: {q30:.2f}, q=50: {q50:.2f}")

            if name == 'quadratic' and popt[0] < 0:
                q_peak = -popt[1] / (2 * popt[0])
                alpha_peak = func(q_peak, *popt)
                print(f"    *** PEAK at q={q_peak:.1f}, alpha_max={alpha_peak:.2f} — UNPHYSICAL")

            fit_results[name] = {
                'params': {pn: float(pv) for pn, pv in zip(param_names, popt)},
                'rms': float(rms),
                'aic': float(aic),
                'residuals': {int(q): float(r) for q, r in zip(qs, residuals)},
                'q30_pred': float(q30),
                'q50_pred': float(q50),
            }
        except Exception as e:
            print(f"  {name:14s}: FAILED — {e}")

    # AIC ranking
    print(f"\n  AIC ranking ({label}):")
    sorted_models = sorted(fit_results.items(), key=lambda x: x[1]['aic'])
    best_aic = sorted_models[0][1]['aic']
    for name, fr in sorted_models:
        delta = fr['aic'] - best_aic
        print(f"    {name:14s}: AIC={fr['aic']:7.2f}, dAIC={delta:+.1f}")

    results[f'fit_results_{label}'] = fit_results
    results[f'best_model_{label}'] = sorted_models[0][0]
    results[f'aic_ranking_{label}'] = [(name, fr['aic']) for name, fr in sorted_models]

# Component analysis
print(f"\n{'='*70}")
print("Component z_m(q) and beta_me(q) fits:")
print("-" * 70)

for comp_name, comp_vals in [('z_m', z_ms), ('beta_me', beta_mes)]:
    print(f"\n  {comp_name}:")
    for func, fname, n_p in [(logarithmic, 'log', 2), (power_law, 'power', 2), (linear, 'linear', 2)]:
        try:
            popt, _ = curve_fit(func, qs, comp_vals, p0=[0.5, 0.5], maxfev=10000)
            pred = func(qs, *popt)
            rms = np.sqrt(np.mean((comp_vals - pred)**2))
            ss = np.sum((comp_vals - pred)**2)
            aic = len(qs) * np.log(ss / len(qs)) + 2 * n_p
            print(f"    {fname:8s}: {popt[0]:.4f}, {popt[1]:.4f}  (RMS={rms:.4f}, AIC={aic:.2f})")
        except Exception as e:
            print(f"    {fname:8s}: FAILED — {e}")

results['qs'] = [int(q) for q in qs]
results['alphas_global'] = [float(a) for a in alphas_global]
results['alphas_last_pair'] = [float(a) for a in alphas_last]
results['z_ms'] = [float(z) for z in z_ms]
results['beta_mes'] = [float(b) for b in beta_mes]

save()
print("\nResults saved.")
