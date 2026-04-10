"""Sprint 128c: Asymptotic extrapolation of chi_F exponents.

Fit alpha_eff(N) = alpha_inf + c/N^p to extract converged alpha.
Also refit alpha(q) with q=6 added.

Uses all exact chi_F data from Sprints 126-128.
"""
import numpy as np
import json, time, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from fss_utils import fit_power_law, pairwise_exponents
from db_utils import record

results = {
    'experiment': '128c_asymptotic_extrapolation',
    'sprint': 128,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'extrapolations': {},
    'alpha_q_refit': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results',
                           'sprint_128c_asymptotic_extrapolation.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)


# All exact chi_F data (Sprints 126-128)
all_data = {
    ('sq', 3): {4: 4.086420, 6: 7.717202, 8: 11.868345, 10: 16.456744, 12: 21.427516, 14: 26.749478},
    ('sq', 4): {4: 12.839321, 6: 27.138030, 8: 45.527791, 9: 56.571892, 10: 67.784201, 11: 80.413},
    ('sq', 5): {4: 30.943214, 6: 72.624695, 8: 132.424150, 9: 171.217, 10: 215.741},
    ('sq', 6): {4: 63.565710, 5: 107.344523, 6: 164.724878, 7: 237.086100, 8: 325.804672, 9: 432.259953},
    ('sq', 7): {4: 117.243073, 6: 333.985930, 7: 501.054584, 8: 703.901},
    ('hybrid', 3): {4: 4.094578, 6: 7.732527, 8: 11.891721, 10: 16.488781, 12: 21.468582},
    ('hybrid', 4): {4: 4.312760, 6: 8.354260, 8: 13.064844, 9: 15.617, 10: 18.341038, 11: 21.219},
    ('hybrid', 5): {4: 4.003316, 6: 7.390136, 8: 10.988427, 9: 12.687, 10: 14.417},
    ('hybrid', 6): {4: 3.465667, 5: 4.702946, 6: 5.889064, 7: 7.011443, 8: 8.066003},
    ('hybrid', 7): {4: 2.775864, 6: 4.330642, 7: 4.943292, 8: 5.519},
}

# Load Sprint 127 data from results file to get exact values
try:
    with open(os.path.join(os.path.dirname(__file__), '..', 'results',
              'sprint_127a_extended_exact_chif.json')) as f:
        s127 = json.load(f)
    # Update from actual Sprint 127 data
    for label, fit_data in s127.get('combined_fits', {}).items():
        if 'sizes' in fit_data and 'chi_F' in fit_data:
            parts = label.split('_')
            if len(parts) == 2 and parts[1].startswith('q'):
                model = parts[0]
                q = int(parts[1][1:])
                key = (model, q)
                if key in all_data:
                    for s, c in zip(fit_data['sizes'], fit_data['chi_F']):
                        all_data[key][s] = c
    print("Loaded Sprint 127 data from results file.")
except Exception as e:
    print(f"Note: Could not load Sprint 127 results: {e}")
    print("Using hardcoded values.")

# Load Sprint 128a/b data
for fname, model, q_val in [
    ('sprint_128a_sq_q6_chif.json', 'sq', 6),
    ('sprint_128b_hybrid_q6_gc_chif.json', 'hybrid', 6),
]:
    try:
        with open(os.path.join(os.path.dirname(__file__), '..', 'results', fname)) as f:
            data = json.load(f)
        key = (model, q_val)
        if key not in all_data:
            all_data[key] = {}
        for d in data.get('data', []) + data.get('chi_F_data', []):
            all_data[key][d['n']] = d['chi_F_exact']
        print(f"Loaded {fname}")
    except Exception as e:
        print(f"Note: Could not load {fname}: {e}")


print("=" * 80)
print("Sprint 128c: Asymptotic extrapolation of chi_F exponents")
print("=" * 80)

# Phase 1: Pairwise alpha_eff and extrapolation
print("\nPhase 1: Pairwise alpha_eff(N) analysis")
print("-" * 80)

from scipy.optimize import curve_fit

def alpha_correction_model(N, alpha_inf, c, p):
    """alpha_eff(N) = alpha_inf + c / N^p"""
    return alpha_inf + c / np.power(N, p)

for (model, q), data in sorted(all_data.items()):
    sizes = np.array(sorted(data.keys()), dtype=float)
    chi_vals = np.array([data[int(n)] for n in sizes])

    if len(sizes) < 4:
        continue

    label = f"{model}_q{q}"
    print(f"\n  {model.upper()} q={q}: sizes={[int(s) for s in sizes]}")

    # Compute pairwise alpha_eff at midpoint N
    pairs = pairwise_exponents(sizes, chi_vals)
    N_mid = []
    alpha_eff = []
    for p in pairs:
        mid = (p['n1'] + p['n2']) / 2.0
        N_mid.append(mid)
        alpha_eff.append(p['alpha'])
        print(f"    ({p['n1']},{p['n2']}): N_mid={mid:.1f}, alpha_eff={p['alpha']:.4f}")

    N_mid = np.array(N_mid)
    alpha_eff = np.array(alpha_eff)

    # Try extrapolation: alpha_eff = alpha_inf + c/N^p
    if len(N_mid) >= 3:
        try:
            # Good initial guess: alpha_inf = last pairwise, c from trend
            a0 = alpha_eff[-1]
            c0 = (alpha_eff[0] - alpha_eff[-1]) * N_mid[-1]
            popt, pcov = curve_fit(alpha_correction_model, N_mid, alpha_eff,
                                    p0=[a0, c0, 1.0],
                                    bounds=([0, -100, 0.1], [10, 100, 5]),
                                    maxfev=5000)
            perr = np.sqrt(np.diag(pcov))
            alpha_inf, c_corr, p_corr = popt
            alpha_inf_err = perr[0]

            predicted = alpha_correction_model(N_mid, *popt)
            ss_res = np.sum((alpha_eff - predicted)**2)
            ss_tot = np.sum((alpha_eff - alpha_eff.mean())**2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

            print(f"    Extrapolation: alpha_inf = {alpha_inf:.4f} +/- {alpha_inf_err:.4f}")
            print(f"      correction: {c_corr:.4f} / N^{p_corr:.2f}, R2={r2:.6f}")

            results['extrapolations'][label] = {
                'alpha_inf': float(alpha_inf),
                'alpha_inf_err': float(alpha_inf_err),
                'c': float(c_corr),
                'p': float(p_corr),
                'r_squared': float(r2),
                'N_mid': [float(n) for n in N_mid],
                'alpha_eff': [float(a) for a in alpha_eff],
            }
        except Exception as e:
            print(f"    Extrapolation failed: {e}")

            # Fallback: linear extrapolation in 1/N
            try:
                inv_N = 1.0 / N_mid
                coeffs = np.polyfit(inv_N, alpha_eff, 1)
                alpha_inf_lin = coeffs[1]
                print(f"    Linear 1/N extrapolation: alpha_inf ~ {alpha_inf_lin:.4f}")
                results['extrapolations'][label] = {
                    'alpha_inf_linear': float(alpha_inf_lin),
                    'method': 'linear_1/N',
                }
            except:
                pass

    # Also do standard power-law fit
    fit = fit_power_law(sizes, chi_vals)
    print(f"    Power-law fit: alpha = {fit['alpha']:.4f} +/- {fit['alpha_err']:.4f}")

save()

# Phase 2: alpha(q) refit with q=6
print(f"\n{'='*80}")
print("Phase 2: alpha(q) refit with q=6 added")
print(f"{'='*80}")

for model_name in ['sq', 'hybrid']:
    qs = []
    alphas = []
    errs = []
    alpha_infs = []  # from extrapolation

    for (m, q), data in sorted(all_data.items()):
        if m != model_name:
            continue
        sizes = np.array(sorted(data.keys()), dtype=float)
        chi_vals = np.array([data[int(n)] for n in sizes])
        if len(sizes) >= 3:
            fit = fit_power_law(sizes, chi_vals)
            qs.append(q)
            alphas.append(fit['alpha'])
            errs.append(fit['alpha_err'])

            # Get extrapolated alpha if available
            label = f"{model_name}_q{q}"
            if label in results['extrapolations'] and 'alpha_inf' in results['extrapolations'][label]:
                alpha_infs.append(results['extrapolations'][label]['alpha_inf'])
            else:
                alpha_infs.append(None)

    qs = np.array(qs, dtype=float)
    alphas = np.array(alphas)
    errs = np.array(errs)

    print(f"\n  {model_name.upper()} alpha(q) with q=6:")
    for i, q in enumerate(qs):
        extrap_str = f", alpha_inf={alpha_infs[i]:.4f}" if alpha_infs[i] is not None else ""
        print(f"    q={int(q)}: alpha={alphas[i]:.4f}+/-{errs[i]:.4f}{extrap_str}")

    if len(qs) >= 3:
        log_q = np.log(qs)
        weights = 1.0 / errs

        # Weighted linear regression: alpha = a*ln(q) + b
        W = np.diag(weights**2)
        X = np.column_stack([log_q, np.ones_like(log_q)])
        beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ alphas)
        residuals = alphas - X @ beta
        chi2 = np.sum((residuals * weights)**2)
        dof = len(qs) - 2

        print(f"\n  Log fit: alpha(q) = {beta[0]:.4f} * ln(q) + ({beta[1]:.4f})")
        print(f"  chi2/dof = {chi2/dof:.4f}")
        for i, q in enumerate(qs):
            pred = X[i] @ beta
            print(f"    q={int(q)}: meas={alphas[i]:.4f}, pred={pred:.4f}, resid={residuals[i]:+.4f}")

        results['alpha_q_refit'][f'{model_name}_log'] = {
            'a': float(beta[0]), 'b': float(beta[1]),
            'formula': f'alpha(q) = {beta[0]:.4f} * ln(q) + ({beta[1]:.4f})',
            'chi2_dof': float(chi2/dof),
            'qs': [int(q) for q in qs],
            'alphas': [float(a) for a in alphas],
            'errors': [float(e) for e in errs],
        }

        # Quadratic in ln(q): alpha = a*ln(q)^2 + b*ln(q) + c
        X3 = np.column_stack([log_q**2, log_q, np.ones_like(log_q)])
        if len(qs) > 3:
            beta3 = np.linalg.solve(X3.T @ W @ X3, X3.T @ W @ alphas)
            pred3 = X3 @ beta3
            resid3 = alphas - pred3
            chi2_3 = np.sum((resid3 * weights)**2)
            dof3 = len(qs) - 3
            print(f"\n  Quadratic in ln(q): alpha(q) = {beta3[0]:.4f}*ln(q)^2 + {beta3[1]:.4f}*ln(q) + ({beta3[2]:.4f})")
            print(f"  chi2/dof = {chi2_3/dof3:.4f}")

            results['alpha_q_refit'][f'{model_name}_quadlog'] = {
                'a': float(beta3[0]), 'b': float(beta3[1]), 'c': float(beta3[2]),
                'chi2_dof': float(chi2_3/dof3) if dof3 > 0 else None,
            }

save()

# Summary
print(f"\n{'='*80}")
print("SUMMARY")
print(f"{'='*80}")
print(f"\n{'Model':>8} {'q':>3} {'alpha':>14} {'alpha_inf':>14} {'drift':>10}")
for (model, q), data in sorted(all_data.items()):
    sizes = np.array(sorted(data.keys()), dtype=float)
    chi_vals = np.array([data[int(n)] for n in sizes])
    if len(sizes) >= 3:
        fit = fit_power_law(sizes, chi_vals)
        label = f"{model}_q{q}"
        alpha_str = f"{fit['alpha']:.4f}+/-{fit['alpha_err']:.4f}"
        extrap = results['extrapolations'].get(label, {})
        if 'alpha_inf' in extrap:
            inf_str = f"{extrap['alpha_inf']:.4f}+/-{extrap['alpha_inf_err']:.4f}"
        else:
            inf_str = "---"

        # Drift direction from pairwise
        pairs = pairwise_exponents(sizes, chi_vals)
        if len(pairs) >= 2:
            drift = pairs[-1]['alpha'] - pairs[0]['alpha']
            drift_str = f"{drift:+.4f}"
        else:
            drift_str = "---"

        print(f"{model:>8} {q:>3} {alpha_str:>14} {inf_str:>14} {drift_str:>10}")

save()
print("\nResults saved.")
