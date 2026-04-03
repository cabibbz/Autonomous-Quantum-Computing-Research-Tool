"""Sprint 114b: Refit alpha(q) with q=15 data point.

8 data points (q=5-15). Tests linear, quadratic, power-law, logarithmic, sqrt.
q=15 measured in exp_114a: alpha_global=3.921, alpha_last_pair=4.070.
All prior models overpredicted — alpha(q) sublinearity is stronger than expected.
"""
import numpy as np
import json, time, os
from scipy.optimize import curve_fit

results = {
    'experiment': '114b_alpha_refit_q15',
    'sprint': 114,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_114b_alpha_refit_q15.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

# Collected alpha data: q, alpha_global, alpha_last_pair
# From sprints 103, 110, 111, 112a, 114a
data = {
    5:  {'alpha_global': 2.091, 'alpha_last_pair': 2.100, 'z_m': 1.028, 'beta_me': 0.125},
    6:  {'alpha_global': 2.377, 'alpha_last_pair': 2.402, 'z_m': 1.072, 'beta_me': 0.260},
    7:  {'alpha_global': 2.650, 'alpha_last_pair': 2.670, 'z_m': 1.238, 'beta_me': 0.418},
    8:  {'alpha_global': 2.897, 'alpha_last_pair': 2.897, 'z_m': 1.320, 'beta_me': 0.564},
    9:  {'alpha_global': 3.159, 'alpha_last_pair': 3.159, 'z_m': 1.438, 'beta_me': 0.709},
    10: {'alpha_global': 3.349, 'alpha_last_pair': 3.417, 'z_m': 1.566, 'beta_me': 1.285},
    12: {'alpha_global': 3.631, 'alpha_last_pair': 3.734, 'z_m': 1.662, 'beta_me': 1.307},
    15: {'alpha_global': 3.921, 'alpha_last_pair': 4.070, 'z_m': 1.784, 'beta_me': 1.353},
}

qs = np.array(sorted(data.keys()), dtype=float)
alphas_global = np.array([data[int(q)]['alpha_global'] for q in qs])
alphas_last = np.array([data[int(q)]['alpha_last_pair'] for q in qs])
z_ms = np.array([data[int(q)]['z_m'] for q in qs])
beta_mes = np.array([data[int(q)]['beta_me'] for q in qs])

print("=" * 70)
print("Sprint 114b: alpha(q) refit with q=15")
print("=" * 70)

print("\nData (global alpha / last-pair alpha):")
for q_val in qs:
    d = data[int(q_val)]
    print(f"  q={int(q_val):2d}: alpha_global={d['alpha_global']:.3f}, "
          f"alpha_last={d['alpha_last_pair']:.3f}, z_m={d['z_m']:.3f}, beta_me={d['beta_me']:.3f}")

# Fit both global and last-pair
for alphas, label in [(alphas_global, "global"), (alphas_last, "last-pair")]:
    print(f"\n{'='*70}")
    print(f"Model fits ({label} alpha, {len(qs)} points):")
    print("-" * 70)

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

    models = {
        'linear':      (linear, 2, ['a', 'b']),
        'quadratic':   (quadratic, 3, ['a', 'b', 'c']),
        'power_law':   (power_law, 2, ['a', 'b']),
        'logarithmic': (logarithmic, 2, ['a', 'b']),
        'sqrt':        (sqrt_law, 2, ['a', 'b']),
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
            print(f"  {name:12s}: RMS={rms:.4f}, AIC={aic:7.2f}  ({params_str})")

            # Predictions at q=20, 25
            q20 = func(20.0, *popt)
            q25 = func(25.0, *popt)
            print(f"    -> q=20: {q20:.2f}, q=25: {q25:.2f}")

            fit_results[name] = {
                'params': {pn: float(pv) for pn, pv in zip(param_names, popt)},
                'rms': float(rms),
                'aic': float(aic),
                'residuals': {int(q): float(r) for q, r in zip(qs, residuals)},
                'q20_pred': float(q20),
                'q25_pred': float(q25),
            }
        except Exception as e:
            print(f"  {name:12s}: FAILED — {e}")

    # AIC ranking
    print(f"\n  AIC ranking ({label}):")
    sorted_models = sorted(fit_results.items(), key=lambda x: x[1]['aic'])
    best_aic = sorted_models[0][1]['aic']
    for name, fr in sorted_models:
        delta = fr['aic'] - best_aic
        print(f"    {name:12s}: AIC={fr['aic']:7.2f}, ΔAIC={delta:+.1f}")

    results[f'fit_results_{label}'] = fit_results
    results[f'best_model_{label}'] = sorted_models[0][0]
    results[f'aic_ranking_{label}'] = [(name, fr['aic']) for name, fr in sorted_models]

# Component analysis (using global alpha)
print(f"\n{'='*70}")
print("Component z_m(q) and beta_me(q) fits:")
print("-" * 70)

# z_m
pz_lin = np.polyfit(qs, z_ms, 1)
pz_pred = np.polyval(pz_lin, qs)
pz_rms = np.sqrt(np.mean((z_ms - pz_pred)**2))
print(f"  z_m linear:    {pz_lin[0]:.4f}*q + {pz_lin[1]:.4f}  (RMS={pz_rms:.4f})")

try:
    pz_pow, _ = curve_fit(power_law, qs, z_ms)
    pz_pow_pred = power_law(qs, *pz_pow)
    pz_pow_rms = np.sqrt(np.mean((z_ms - pz_pow_pred)**2))
    print(f"  z_m power-law: {pz_pow[0]:.4f}*q^{pz_pow[1]:.4f}  (RMS={pz_pow_rms:.4f})")
except:
    pass

try:
    pz_log, _ = curve_fit(logarithmic, qs, z_ms)
    pz_log_pred = logarithmic(qs, *pz_log)
    pz_log_rms = np.sqrt(np.mean((z_ms - pz_log_pred)**2))
    print(f"  z_m log:       {pz_log[0]:.4f}*ln(q) + {pz_log[1]:.4f}  (RMS={pz_log_rms:.4f})")
except:
    pass

# beta_me
pb_lin = np.polyfit(qs, beta_mes, 1)
pb_pred = np.polyval(pb_lin, qs)
pb_rms = np.sqrt(np.mean((beta_mes - pb_pred)**2))
print(f"  beta_me linear:    {pb_lin[0]:.4f}*q + {pb_lin[1]:.4f}  (RMS={pb_rms:.4f})")

try:
    pb_pow, _ = curve_fit(power_law, qs, beta_mes, p0=[0.1, 1.0])
    pb_pow_pred = power_law(qs, *pb_pow)
    pb_pow_rms = np.sqrt(np.mean((beta_mes - pb_pow_pred)**2))
    print(f"  beta_me power-law: {pb_pow[0]:.4f}*q^{pb_pow[1]:.4f}  (RMS={pb_pow_rms:.4f})")
except:
    pass

try:
    pb_log, _ = curve_fit(logarithmic, qs, beta_mes)
    pb_log_pred = logarithmic(qs, *pb_log)
    pb_log_rms = np.sqrt(np.mean((beta_mes - pb_log_pred)**2))
    print(f"  beta_me log:       {pb_log[0]:.4f}*ln(q) + {pb_log[1]:.4f}  (RMS={pb_log_rms:.4f})")
except:
    pass

# Reconstructed alpha from components
alpha_recon = pb_lin[0]*qs + pb_lin[1] + 2*(pz_lin[0]*qs + pz_lin[1]) - 1
recon_rms = np.sqrt(np.mean((alphas_global - alpha_recon)**2))
print(f"\n  Reconstructed alpha (linear components): "
      f"{2*pz_lin[0]+pb_lin[0]:.4f}*q + {2*pz_lin[1]+pb_lin[1]-1:.4f}")
print(f"  Reconstruction RMS vs measured: {recon_rms:.4f}")

results['qs'] = [int(q) for q in qs]
results['alphas_global'] = [float(a) for a in alphas_global]
results['alphas_last_pair'] = [float(a) for a in alphas_last]
results['z_ms'] = [float(z) for z in z_ms]
results['beta_mes'] = [float(b) for b in beta_mes]
results['component_fits'] = {
    'z_m_linear': {'slope': float(pz_lin[0]), 'intercept': float(pz_lin[1]), 'rms': float(pz_rms)},
    'beta_me_linear': {'slope': float(pb_lin[0]), 'intercept': float(pb_lin[1]), 'rms': float(pb_rms)},
}

save()
print("\nResults saved.")
