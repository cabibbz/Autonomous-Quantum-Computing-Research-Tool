"""Sprint 122a: Log correction fits for chi_F at q=4.

S_q q=4 has a marginal operator → multiplicative log corrections:
  chi_F = A * N^alpha_inf * (ln N)^{-p}

Compare three models:
  1. Pure power law: chi_F = A * N^alpha
  2. Log-corrected: chi_F = A * N^2 * (ln N)^{-p}  (fixed alpha=2)
  3. Free log-corrected: chi_F = A * N^alpha * (ln N)^{-p}

Apply to both S_q and hybrid data at q=4.
"""
import numpy as np
import json
from lmfit import Model, Parameters

# === DATA ===
# S_q q=4 from Sprint 121b
sq_n = np.array([4, 5, 6, 7, 8, 9, 10, 11], dtype=float)
sq_chi = np.array([6.340, 9.527, 13.211, 17.383, 22.030, 27.142, 32.710, 38.725])

# Hybrid q=4 from Sprint 120b
hy_n = np.array([4, 5, 6, 7, 8, 9, 10, 11], dtype=float)
hy_chi = np.array([2.071, 2.984, 3.978, 5.047, 6.184, 7.385, 8.647, 9.965])

def pure_power(N, A, alpha):
    return A * N**alpha

def log_corrected_fixed(N, A, p):
    """chi_F = A * N^2 * (ln N)^{-p}"""
    return A * N**2 * np.log(N)**(-p)

def log_corrected_free(N, A, alpha, p):
    """chi_F = A * N^alpha * (ln N)^{-p}"""
    return A * N**alpha * np.log(N)**(-p)

def fit_and_report(label, sizes, chi_F):
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    print(f"  N = {sizes.tolist()}")
    print(f"  chi_F = {[f'{x:.4f}' for x in chi_F]}")

    log_chi = np.log(chi_F)
    log_N = np.log(sizes)
    results = {}

    # 1. Pure power law
    m = Model(pure_power)
    p = m.make_params(A=0.5, alpha=1.7)
    p['A'].min = 0.01
    r = m.fit(chi_F, p, N=sizes)
    ss_res = np.sum(r.residual**2)
    ss_tot = np.sum((chi_F - chi_F.mean())**2)
    r2 = 1 - ss_res/ss_tot
    # Also compute in log space
    pred = r.best_fit
    log_ss_res = np.sum((np.log(pred) - log_chi)**2)
    log_ss_tot = np.sum((log_chi - log_chi.mean())**2)
    log_r2 = 1 - log_ss_res/log_ss_tot

    results['pure_power'] = {
        'alpha': float(r.params['alpha'].value),
        'alpha_err': float(r.params['alpha'].stderr) if r.params['alpha'].stderr else 0,
        'A': float(r.params['A'].value),
        'r2': float(r2),
        'log_r2': float(log_r2),
        'aic': float(r.aic),
        'residuals': [float(x) for x in r.residual],
    }
    print(f"\n  [1] Pure power: chi_F = {r.params['A'].value:.4f} * N^{r.params['alpha'].value:.4f}")
    print(f"      R² = {r2:.8f}, log-R² = {log_r2:.8f}, AIC = {r.aic:.2f}")

    # 2. Log-corrected, alpha=2 fixed
    m = Model(log_corrected_fixed)
    p = m.make_params(A=0.5, p=1.0)
    p['A'].min = 0.001
    p['p'].min = -5
    p['p'].max = 10
    r = m.fit(chi_F, p, N=sizes)
    ss_res = np.sum(r.residual**2)
    r2 = 1 - ss_res/ss_tot
    pred = r.best_fit
    log_ss_res = np.sum((np.log(pred) - log_chi)**2)
    log_r2 = 1 - log_ss_res/log_ss_tot

    results['log_fixed_alpha2'] = {
        'p': float(r.params['p'].value),
        'p_err': float(r.params['p'].stderr) if r.params['p'].stderr else 0,
        'A': float(r.params['A'].value),
        'r2': float(r2),
        'log_r2': float(log_r2),
        'aic': float(r.aic),
    }
    print(f"\n  [2] Log-corrected (alpha=2 fixed): chi_F = {r.params['A'].value:.4f} * N^2 * (ln N)^{{{-r.params['p'].value:.4f}}}")
    print(f"      p = {r.params['p'].value:.4f} ± {r.params['p'].stderr if r.params['p'].stderr else 0:.4f}")
    print(f"      R² = {r2:.8f}, log-R² = {log_r2:.8f}, AIC = {r.aic:.2f}")

    # 3. Log-corrected, alpha free
    m = Model(log_corrected_free)
    p = m.make_params(A=0.5, alpha=1.8, p=0.5)
    p['A'].min = 0.001
    p['alpha'].min = 0.1
    p['alpha'].max = 4.0
    p['p'].min = -10
    p['p'].max = 20
    r = m.fit(chi_F, p, N=sizes)
    ss_res = np.sum(r.residual**2)
    r2 = 1 - ss_res/ss_tot
    pred = r.best_fit
    log_ss_res = np.sum((np.log(pred) - log_chi)**2)
    log_r2 = 1 - log_ss_res/log_ss_tot

    results['log_free'] = {
        'alpha': float(r.params['alpha'].value),
        'alpha_err': float(r.params['alpha'].stderr) if r.params['alpha'].stderr else 0,
        'p': float(r.params['p'].value),
        'p_err': float(r.params['p'].stderr) if r.params['p'].stderr else 0,
        'A': float(r.params['A'].value),
        'r2': float(r2),
        'log_r2': float(log_r2),
        'aic': float(r.aic),
    }
    print(f"\n  [3] Log-corrected (alpha free): chi_F = {r.params['A'].value:.4f} * N^{r.params['alpha'].value:.4f} * (ln N)^{{{-r.params['p'].value:.4f}}}")
    print(f"      alpha = {r.params['alpha'].value:.4f} ± {r.params['alpha'].stderr if r.params['alpha'].stderr else 0:.4f}")
    print(f"      p = {r.params['p'].value:.4f} ± {r.params['p'].stderr if r.params['p'].stderr else 0:.4f}")
    print(f"      R² = {r2:.8f}, log-R² = {log_r2:.8f}, AIC = {r.aic:.2f}")

    # 4. Power-law with 1/N^2 correction: chi_F = A*N^alpha*(1 + B/N^2)
    def power_corr(N, A, alpha, B):
        return A * N**alpha * (1 + B/N**2)
    m = Model(power_corr)
    p = m.make_params(A=0.5, alpha=1.8, B=-0.5)
    p['A'].min = 0.001
    r = m.fit(chi_F, p, N=sizes)
    ss_res = np.sum(r.residual**2)
    r2 = 1 - ss_res/ss_tot
    pred = r.best_fit
    log_ss_res = np.sum((np.log(pred) - log_chi)**2)
    log_r2 = 1 - log_ss_res/log_ss_tot

    results['power_corr'] = {
        'alpha': float(r.params['alpha'].value),
        'alpha_err': float(r.params['alpha'].stderr) if r.params['alpha'].stderr else 0,
        'B': float(r.params['B'].value),
        'A': float(r.params['A'].value),
        'r2': float(r2),
        'log_r2': float(log_r2),
        'aic': float(r.aic),
    }
    print(f"\n  [4] Power + 1/N² corr: chi_F = {r.params['A'].value:.4f} * N^{r.params['alpha'].value:.4f} * (1 + {r.params['B'].value:.4f}/N²)")
    print(f"      alpha = {r.params['alpha'].value:.4f} ± {r.params['alpha'].stderr if r.params['alpha'].stderr else 0:.4f}")
    print(f"      R² = {r2:.8f}, log-R² = {log_r2:.8f}, AIC = {r.aic:.2f}")

    # Summary
    print(f"\n  --- Model comparison (AIC) ---")
    all_models = sorted(results.items(), key=lambda x: x[1]['aic'])
    best_aic = all_models[0][1]['aic']
    for name, r in all_models:
        daic = r['aic'] - best_aic
        alpha_str = f"alpha={r.get('alpha', 2.0):.4f}" if 'alpha' in r else "alpha=2(fixed)"
        print(f"    {name:25s}: AIC={r['aic']:8.2f} (dAIC={daic:6.2f}), R²={r['r2']:.8f}, {alpha_str}")

    return results

# Run fits
print("SPRINT 122a: LOG CORRECTION FITS AT q=4")
print("="*60)

sq_results = fit_and_report("S_q POTTS q=4 (g_c=0.25)", sq_n, sq_chi)
hy_results = fit_and_report("HYBRID q=4 (g_c=0.393)", hy_n, hy_chi)

# Compare
print(f"\n{'='*60}")
print("COMPARISON: S_q vs Hybrid at q=4")
print("="*60)

print(f"\n  S_q best alpha (pure power): {sq_results['pure_power']['alpha']:.4f}")
print(f"  Hybrid best alpha (pure power): {hy_results['pure_power']['alpha']:.4f}")

if 'p' in sq_results.get('log_fixed_alpha2', {}):
    print(f"\n  S_q log correction p (alpha=2 fixed): {sq_results['log_fixed_alpha2']['p']:.4f}")
if 'p' in hy_results.get('log_free', {}):
    print(f"  Hybrid log correction p (alpha free): {hy_results['log_free']['p']:.4f}, alpha={hy_results['log_free']['alpha']:.4f}")

print(f"\n  S_q: log model (alpha=2) preferred? AIC comparison:")
sq_log_aic = sq_results['log_fixed_alpha2']['aic']
sq_power_aic = sq_results['pure_power']['aic']
print(f"    Pure power AIC = {sq_power_aic:.2f}")
print(f"    Log (alpha=2) AIC = {sq_log_aic:.2f}")
print(f"    dAIC = {sq_power_aic - sq_log_aic:.2f} (positive = log wins)")

print(f"\n  Hybrid: log model vs pure power?")
hy_log_aic = hy_results['log_free']['aic']
hy_power_aic = hy_results['pure_power']['aic']
print(f"    Pure power AIC = {hy_power_aic:.2f}")
print(f"    Log (free) AIC = {hy_log_aic:.2f}")
print(f"    dAIC = {hy_power_aic - hy_log_aic:.2f} (positive = log wins)")

# Key question: does S_q alpha extrapolate to 2.0?
print(f"\n  KEY QUESTION: Does S_q alpha -> 2.0?")
if 'alpha' in sq_results['log_free']:
    print(f"    Free fit: alpha = {sq_results['log_free']['alpha']:.4f} ± {sq_results['log_free']['alpha_err']:.4f}")
    print(f"    Consistent with 2.0? {abs(sq_results['log_free']['alpha'] - 2.0) < 3*sq_results['log_free']['alpha_err'] if sq_results['log_free']['alpha_err'] > 0 else 'N/A'}")
if 'alpha' in sq_results['power_corr']:
    print(f"    Power+1/N² correction: alpha = {sq_results['power_corr']['alpha']:.4f} ± {sq_results['power_corr']['alpha_err']:.4f}")

# Save
all_results = {
    'experiment': '122a_log_corrections_q4',
    'sprint': 122,
    'sq_q4': sq_results,
    'hybrid_q4': hy_results,
}
with open('../results/sprint_122a_log_corrections.json', 'w') as f:
    json.dump(all_results, f, indent=2, default=str)
print(f"\nResults saved to results/sprint_122a_log_corrections.json")
