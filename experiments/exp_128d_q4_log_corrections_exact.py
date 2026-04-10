"""Sprint 128d: Log correction fits for S_q q=4 chi_F using EXACT data.

Sprint 122a did this with SPECTRAL chi_F (systematic negative bias, Sprint 126).
This redo uses exact chi_F from the database.

Four models compared:
  1. Pure power law: chi_F = A * N^alpha
  2. Log-corrected (alpha=2 fixed): chi_F = A * N^2 * (ln N)^{-p}  [Salas-Sokal]
  3. Free log-corrected: chi_F = A * N^alpha * (ln N)^{-p}
  4. Power + 1/N^2: chi_F = A * N^alpha * (1 + B/N^2)

Key question: Salas-Sokal predicts p=3/2 for marginal operator corrections.
Does exact data (2x spectral) change the model ranking?
"""
import numpy as np
import json, time, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from lmfit import Model
from db_utils import query, record
from fss_utils import pairwise_exponents

results = {
    'experiment': '128d_q4_log_corrections_exact',
    'sprint': 128,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'models': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results',
                           'sprint_128d_q4_log_corrections_exact.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)


# Load exact chi_F from DB
rows = query(quantity='chi_F_exact', model='sq', q=4)
data = {}
for r in rows:
    n, val, sprint = r[4], r[6], r[1]
    if n is not None and val is not None:
        if n not in data or sprint > data[n][1]:
            data[n] = (val, sprint)

sizes = np.array(sorted(data.keys()), dtype=float)
chi_F = np.array([data[int(n)][0] for n in sizes])

print("=" * 80)
print("Sprint 128d: S_q q=4 log corrections with EXACT chi_F")
print("=" * 80)
print(f"\nData from DB ({len(sizes)} sizes):")
for i, n in enumerate(sizes):
    print(f"  n={int(n):2d}: chi_F = {chi_F[i]:.6f}")

# Sprint 122a spectral values for comparison
spectral_chi = {4: 6.340, 5: 9.527, 6: 13.211, 7: 17.383,
                8: 22.030, 9: 27.142, 10: 32.710, 11: 38.725}
print(f"\nComparison with Sprint 122a spectral chi_F:")
for n in sorted(set(spectral_chi.keys()) & set(int(s) for s in sizes)):
    exact = data[n][0]
    spec = spectral_chi[n]
    ratio = exact / spec
    print(f"  n={n}: exact={exact:.4f}, spectral={spec:.4f}, ratio={ratio:.3f}")

results['data'] = {'sizes': [int(s) for s in sizes],
                   'chi_F_exact': [float(v) for v in chi_F]}

# Model definitions
def pure_power(N, A, alpha):
    return A * N**alpha

def log_corrected_fixed(N, A, p):
    """alpha=2 fixed, Salas-Sokal form"""
    return A * N**2 * np.log(N)**(-p)

def log_corrected_free(N, A, alpha, p):
    return A * N**alpha * np.log(N)**(-p)

def power_plus_correction(N, A, alpha, B):
    return A * N**alpha * (1.0 + B / N**2)


def fit_model(name, func, sizes, chi_F, param_dict, bounds_dict=None):
    """Fit a model, return results dict."""
    m = Model(func)
    params = m.make_params(**param_dict)
    if bounds_dict:
        for pname, (lo, hi) in bounds_dict.items():
            if lo is not None: params[pname].min = lo
            if hi is not None: params[pname].max = hi

    r = m.fit(chi_F, params, N=sizes)

    ss_res = np.sum(r.residual**2)
    ss_tot = np.sum((chi_F - chi_F.mean())**2)
    r2 = 1 - ss_res / ss_tot

    log_pred = np.log(r.best_fit)
    log_obs = np.log(chi_F)
    log_ss_res = np.sum((log_pred - log_obs)**2)
    log_ss_tot = np.sum((log_obs - log_obs.mean())**2)
    log_r2 = 1 - log_ss_res / log_ss_tot

    result = {
        'params': {},
        'r_squared': float(r2),
        'log_r_squared': float(log_r2),
        'aic': float(r.aic),
        'bic': float(r.bic),
        'n_params': r.nvarys,
        'residuals': [float(x) for x in r.residual],
    }
    for pname in r.params:
        p = r.params[pname]
        result['params'][pname] = {
            'value': float(p.value),
            'stderr': float(p.stderr) if p.stderr else None,
        }

    return result


print(f"\n{'='*80}")
print("Model comparison (4 models)")
print(f"{'='*80}")

# 1. Pure power law
print("\n[1] Pure power law: chi_F = A * N^alpha")
r1 = fit_model('pure_power', pure_power, sizes, chi_F,
               {'A': 0.5, 'alpha': 1.8}, {'A': (0.001, None)})
a = r1['params']['alpha']
print(f"    alpha = {a['value']:.4f} +/- {a['stderr']:.4f}")
print(f"    R2 = {r1['r_squared']:.8f}, log-R2 = {r1['log_r_squared']:.8f}")
print(f"    AIC = {r1['aic']:.2f}")
results['models']['pure_power'] = r1

# 2. Log-corrected alpha=2 fixed (Salas-Sokal)
print("\n[2] Log-corrected (alpha=2 fixed): chi_F = A * N^2 * (ln N)^{-p}")
print("    Salas-Sokal prediction: p = 3/2 = 1.5")
r2 = fit_model('log_fixed', log_corrected_fixed, sizes, chi_F,
               {'A': 0.5, 'p': 1.5}, {'A': (0.001, None), 'p': (-5, 10)})
p_val = r2['params']['p']
print(f"    p = {p_val['value']:.4f} +/- {p_val['stderr']:.4f}  (predicted: 1.500)")
print(f"    R2 = {r2['r_squared']:.8f}, log-R2 = {r2['log_r_squared']:.8f}")
print(f"    AIC = {r2['aic']:.2f}")
results['models']['log_corrected_alpha2'] = r2

# 3. Free log-corrected
print("\n[3] Free log-corrected: chi_F = A * N^alpha * (ln N)^{-p}")
r3 = fit_model('log_free', log_corrected_free, sizes, chi_F,
               {'A': 0.5, 'alpha': 1.8, 'p': 0.5},
               {'A': (0.001, None), 'p': (-5, 10)})
a3 = r3['params']['alpha']
p3 = r3['params']['p']
print(f"    alpha = {a3['value']:.4f} +/- {a3['stderr']:.4f}")
print(f"    p = {p3['value']:.4f} +/- {p3['stderr']:.4f}")
print(f"    R2 = {r3['r_squared']:.8f}, log-R2 = {r3['log_r_squared']:.8f}")
print(f"    AIC = {r3['aic']:.2f}")
results['models']['log_corrected_free'] = r3

# 4. Power + 1/N^2 correction
print("\n[4] Power + 1/N^2: chi_F = A * N^alpha * (1 + B/N^2)")
r4 = fit_model('power_corr', power_plus_correction, sizes, chi_F,
               {'A': 0.5, 'alpha': 1.8, 'B': 1.0},
               {'A': (0.001, None)})
a4 = r4['params']['alpha']
B4 = r4['params']['B']
print(f"    alpha = {a4['value']:.4f} +/- {a4['stderr']:.4f}")
print(f"    B = {B4['value']:.4f} +/- {B4['stderr']:.4f}")
print(f"    R2 = {r4['r_squared']:.8f}, log-R2 = {r4['log_r_squared']:.8f}")
print(f"    AIC = {r4['aic']:.2f}")
results['models']['power_plus_correction'] = r4

# AIC comparison
print(f"\n{'='*80}")
print("AIC COMPARISON (lower = better)")
print(f"{'='*80}")
models = [
    ('Pure power', r1),
    ('Log alpha=2', r2),
    ('Log free', r3),
    ('Power+1/N^2', r4),
]
aics = [(name, m['aic']) for name, m in models]
aics.sort(key=lambda x: x[1])
best_aic = aics[0][1]
print(f"\n  {'Model':<15} {'AIC':>10} {'dAIC':>10} {'R2':>12} {'log-R2':>12}")
for name, aic in aics:
    m = [m for n, m in models if n == name][0]
    print(f"  {name:<15} {aic:>10.2f} {aic-best_aic:>10.2f} {m['r_squared']:>12.8f} {m['log_r_squared']:>12.8f}")

results['aic_ranking'] = [{'model': name, 'aic': aic, 'daic': aic - best_aic}
                          for name, aic in aics]

# Comparison with Sprint 122a (spectral data)
print(f"\n{'='*80}")
print("COMPARISON: Sprint 122a (spectral) vs Sprint 128d (exact)")
print(f"{'='*80}")
s122a_results = {
    'pure_power': {'alpha': 1.777, 'aic': -50.13},
    'log_alpha2': {'p': 0.448, 'aic': -30.28},
    'log_free': {'alpha': 1.693, 'p': -0.172, 'aic': -66.48},
    'power_corr': {'alpha': 1.757, 'aic': -72.20},
}
print(f"\n  {'Model':<15} {'S122a AIC':>10} {'S128d AIC':>10} {'Winner':>10}")
for name_122, name_128, label in [
    ('pure_power', r1, 'Pure power'),
    ('log_alpha2', r2, 'Log alpha=2'),
    ('log_free', r3, 'Log free'),
    ('power_corr', r4, 'Power+1/N^2'),
]:
    old_aic = s122a_results[name_122]['aic']
    new_aic = name_128['aic']
    print(f"  {label:<15} {old_aic:>10.2f} {new_aic:>10.2f}")

print(f"\n  S122a used SPECTRAL chi_F (biased low by factor ~2)")
print(f"  S128d uses EXACT chi_F (finite-difference, gold standard)")
print(f"  Same ranking = model comparison robust to data source")
print(f"  Different ranking = spectral bias affected conclusions")

# Pairwise exponents
print(f"\n{'='*80}")
print("Pairwise exponents (exact chi_F)")
print(f"{'='*80}")
pairs = pairwise_exponents(sizes, chi_F)
for p in pairs:
    print(f"  ({p['n1']},{p['n2']}): alpha = {p['alpha']:.4f}")

# Record to DB
record(sprint=128, model='sq', q=4, n=int(max(sizes)),
       quantity='alpha_log_analysis', value=float(r1['params']['alpha']['value']),
       method='exact_chif_4model_comparison',
       notes=f"AIC ranking: {[a[0] for a in aics]}")

save()
print(f"\nResults saved.")
