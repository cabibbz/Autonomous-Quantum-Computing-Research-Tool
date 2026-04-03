"""Sprint 115b: Refit alpha(q) with q=20 data point.

9 data points (q=5-20). q=20 measured in exp_115a: alpha_global=4.590.
Quadratic predicted 3.82 — DECISIVELY RULED OUT (off by 0.77).
Logarithmic predicted 4.49 — closest (off by 0.10).
"""
import numpy as np
import json, time, os
from scipy.optimize import curve_fit

results = {
    'experiment': '115b_alpha_refit_q20',
    'sprint': 115,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_115b_alpha_refit_q20.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

# Collected alpha data from sprints 103, 110, 111, 112a, 114a, 115a
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
}

qs = np.array(sorted(data.keys()), dtype=float)
alphas_global = np.array([data[int(q)]['alpha_global'] for q in qs])
alphas_last = np.array([data[int(q)]['alpha_last_pair'] for q in qs])
z_ms = np.array([data[int(q)]['z_m'] for q in qs])
beta_mes = np.array([data[int(q)]['beta_me'] for q in qs])

print("=" * 70)
print("Sprint 115b: alpha(q) refit with q=20 (9 data points)")
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

# Fit both global and last-pair
for alphas, label in [(alphas_global, "global"), (alphas_last, "last-pair")]:
    print(f"\n{'='*70}")
    print(f"Model fits ({label} alpha, {len(qs)} points):")
    print("-" * 70)

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

            # Predictions at q=25, 30
            q25 = func(25.0, *popt)
            q30 = func(30.0, *popt)
            print(f"    -> q=25: {q25:.2f}, q=30: {q30:.2f}")

            # Check for unphysical behavior
            if name == 'quadratic' and popt[0] < 0:
                q_peak = -popt[1] / (2 * popt[0])
                alpha_peak = func(q_peak, *popt)
                print(f"    *** PEAK at q={q_peak:.1f}, alpha_max={alpha_peak:.2f} — {'UNPHYSICAL' if q_peak > 0 else 'OK'}")

            fit_results[name] = {
                'params': {pn: float(pv) for pn, pv in zip(param_names, popt)},
                'rms': float(rms),
                'aic': float(aic),
                'residuals': {int(q): float(r) for q, r in zip(qs, residuals)},
                'q25_pred': float(q25),
                'q30_pred': float(q30),
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

# Component analysis
print(f"\n{'='*70}")
print("Component z_m(q) and beta_me(q) fits:")
print("-" * 70)

for comp_name, comp_vals in [('z_m', z_ms), ('beta_me', beta_mes)]:
    print(f"\n  {comp_name}:")

    # Linear
    p_lin = np.polyfit(qs, comp_vals, 1)
    p_pred = np.polyval(p_lin, qs)
    rms_lin = np.sqrt(np.mean((comp_vals - p_pred)**2))
    ss_lin = np.sum((comp_vals - p_pred)**2)
    aic_lin = len(qs) * np.log(ss_lin / len(qs)) + 2 * 2
    print(f"    linear:    {p_lin[0]:.4f}*q + {p_lin[1]:.4f}  (RMS={rms_lin:.4f}, AIC={aic_lin:.2f})")

    # Logarithmic
    try:
        popt_log, _ = curve_fit(logarithmic, qs, comp_vals)
        pred_log = logarithmic(qs, *popt_log)
        rms_log = np.sqrt(np.mean((comp_vals - pred_log)**2))
        ss_log = np.sum((comp_vals - pred_log)**2)
        aic_log = len(qs) * np.log(ss_log / len(qs)) + 2 * 2
        print(f"    log:       {popt_log[0]:.4f}*ln(q) + {popt_log[1]:.4f}  (RMS={rms_log:.4f}, AIC={aic_log:.2f})")
    except:
        pass

    # Power-law
    try:
        popt_pow, _ = curve_fit(power_law, qs, comp_vals, p0=[0.5, 0.5])
        pred_pow = power_law(qs, *popt_pow)
        rms_pow = np.sqrt(np.mean((comp_vals - pred_pow)**2))
        ss_pow = np.sum((comp_vals - pred_pow)**2)
        aic_pow = len(qs) * np.log(ss_pow / len(qs)) + 2 * 2
        print(f"    power-law: {popt_pow[0]:.4f}*q^{popt_pow[1]:.4f}  (RMS={rms_pow:.4f}, AIC={aic_pow:.2f})")
    except:
        pass

    results[f'{comp_name}_fits'] = {
        'linear': {'slope': float(p_lin[0]), 'intercept': float(p_lin[1]), 'rms': float(rms_lin)},
    }

# Reconstructed alpha from log components
print(f"\n{'='*70}")
print("Reconstructed alpha from component fits:")
print("-" * 70)

try:
    popt_z_log, _ = curve_fit(logarithmic, qs, z_ms)
    popt_b_log, _ = curve_fit(logarithmic, qs, beta_mes)
    alpha_recon_log = popt_b_log[0]*np.log(qs) + popt_b_log[1] + 2*(popt_z_log[0]*np.log(qs) + popt_z_log[1]) - 1
    recon_rms = np.sqrt(np.mean((alphas_global - alpha_recon_log)**2))
    a_coeff = popt_b_log[0] + 2*popt_z_log[0]
    b_coeff = popt_b_log[1] + 2*popt_z_log[1] - 1
    print(f"  alpha_recon = {a_coeff:.4f}*ln(q) + {b_coeff:.4f}")
    print(f"  vs direct log fit of alpha")
    print(f"  Reconstruction RMS: {recon_rms:.4f}")
    print(f"  q=25 prediction (recon): {a_coeff*np.log(25) + b_coeff:.3f}")
    print(f"  q=30 prediction (recon): {a_coeff*np.log(30) + b_coeff:.3f}")
except Exception as e:
    print(f"  FAILED: {e}")

results['qs'] = [int(q) for q in qs]
results['alphas_global'] = [float(a) for a in alphas_global]
results['alphas_last_pair'] = [float(a) for a in alphas_last]
results['z_ms'] = [float(z) for z in z_ms]
results['beta_mes'] = [float(b) for b in beta_mes]

save()
print("\nResults saved.")
