"""Sprint 121a: z_m(q) continuous fit for hybrid model — find q_cross where z_m=1.

Data from Sprints 119-120:
  hybrid z_m: q=[3,4,5,7,10] -> z_m=[1.030, 1.004, 0.909, 0.770, 0.690]
  S_q z_m:    q=[3,4,5,7]    -> z_m=[1.030, 1.050, 1.310, 1.400]

Fit multiple functional forms: linear, log, power-law, rational.
Find q_cross (hybrid z_m = 1) by interpolation and best-fit inversion.
Also fit S_q z_m(q) for comparison.
"""
import numpy as np
import json, time, os
from lmfit import Model, Parameters
from db_utils import record

results = {
    'experiment': '121a_zm_q_crossing',
    'sprint': 121,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results', 'sprint_121a_zm_q_crossing.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

# ---- DATA FROM SPRINTS 119-120 ----
# Hybrid z_m (global fits from spectral decomposition)
hybrid_q   = np.array([3, 4, 5, 7, 10], dtype=float)
hybrid_zm  = np.array([1.030, 1.004, 0.909, 0.770, 0.690])
hybrid_zm_err = np.array([0.010, 0.005, 0.008, 0.015, 0.020])  # approximate from fit errors

# S_q z_m
sq_q  = np.array([3, 4, 5, 7], dtype=float)
sq_zm = np.array([1.030, 1.050, 1.310, 1.400])
sq_zm_err = np.array([0.010, 0.015, 0.015, 0.030])

# Hybrid alpha (for reference)
hybrid_alpha = np.array([1.40, 1.55, 1.41, 0.95, 0.55])
sq_alpha     = np.array([1.40, 1.77, 2.09, 2.65])

print("=" * 70)
print("EXPERIMENT 121a: z_m(q) CROSSING ANALYSIS")
print("=" * 70)

# ---- FIT MODELS ----
def linear_model(x, a, b):
    return a + b * x

def log_model(x, a, b):
    return a + b * np.log(x)

def power_model(x, a, b, c):
    return a + b * np.power(x, c)

def rational_model(x, a, b, c):
    return a + b / (x + c)

# ---- FIT HYBRID z_m(q) ----
print("\n--- HYBRID z_m(q) fits ---\n")

fit_results = {}

# 1. Linear: z_m = a + b*q
m_lin = Model(linear_model)
p_lin = m_lin.make_params(a=1.3, b=-0.05)
r_lin = m_lin.fit(hybrid_zm, x=hybrid_q, params=p_lin, weights=1/hybrid_zm_err)
ss_res = np.sum(r_lin.residual**2 * hybrid_zm_err**2)
ss_tot = np.sum((hybrid_zm - hybrid_zm.mean())**2)
r2_lin = 1 - ss_res/ss_tot if ss_tot > 0 else 0
q_cross_lin = (1.0 - r_lin.params['a'].value) / r_lin.params['b'].value if r_lin.params['b'].value != 0 else None
print(f"Linear: z_m = {r_lin.params['a'].value:.4f} + {r_lin.params['b'].value:.4f}*q")
print(f"  R² = {r2_lin:.5f}, AIC = {r_lin.aic:.2f}")
print(f"  q_cross = {q_cross_lin:.3f}" if q_cross_lin else "  q_cross = undefined")
fit_results['linear'] = {'a': r_lin.params['a'].value, 'b': r_lin.params['b'].value,
                          'r2': r2_lin, 'aic': r_lin.aic, 'q_cross': q_cross_lin}

# 2. Log: z_m = a + b*ln(q)
m_log = Model(log_model)
p_log = m_log.make_params(a=1.5, b=-0.2)
r_log = m_log.fit(hybrid_zm, x=hybrid_q, params=p_log, weights=1/hybrid_zm_err)
ss_res = np.sum(r_log.residual**2 * hybrid_zm_err**2)
r2_log = 1 - ss_res/ss_tot if ss_tot > 0 else 0
q_cross_log = np.exp((1.0 - r_log.params['a'].value) / r_log.params['b'].value) if r_log.params['b'].value != 0 else None
print(f"\nLog: z_m = {r_log.params['a'].value:.4f} + {r_log.params['b'].value:.4f}*ln(q)")
print(f"  R² = {r2_log:.5f}, AIC = {r_log.aic:.2f}")
print(f"  q_cross = {q_cross_log:.3f}" if q_cross_log else "  q_cross = undefined")
fit_results['log'] = {'a': r_log.params['a'].value, 'b': r_log.params['b'].value,
                       'r2': r2_log, 'aic': r_log.aic, 'q_cross': q_cross_log}

# 3. Power: z_m = a + b * q^c
m_pow = Model(power_model)
p_pow = m_pow.make_params(a=0.5, b=1.5, c=-1.0)
p_pow['c'].min = -5
p_pow['c'].max = 0
r_pow = m_pow.fit(hybrid_zm, x=hybrid_q, params=p_pow, weights=1/hybrid_zm_err)
ss_res = np.sum(r_pow.residual**2 * hybrid_zm_err**2)
r2_pow = 1 - ss_res/ss_tot if ss_tot > 0 else 0
# Solve a + b*q^c = 1 numerically
from scipy.optimize import brentq
try:
    a_p, b_p, c_p = r_pow.params['a'].value, r_pow.params['b'].value, r_pow.params['c'].value
    q_cross_pow = brentq(lambda x: a_p + b_p * x**c_p - 1.0, 2.5, 15.0)
except:
    q_cross_pow = None
print(f"\nPower: z_m = {r_pow.params['a'].value:.4f} + {r_pow.params['b'].value:.4f}*q^{r_pow.params['c'].value:.3f}")
print(f"  R² = {r2_pow:.5f}, AIC = {r_pow.aic:.2f}")
print(f"  q_cross = {q_cross_pow:.3f}" if q_cross_pow else "  q_cross = undefined")
fit_results['power'] = {'a': a_p, 'b': b_p, 'c': c_p,
                         'r2': r2_pow, 'aic': r_pow.aic, 'q_cross': q_cross_pow}

# 4. Rational: z_m = a + b/(q+c)
m_rat = Model(rational_model)
p_rat = m_rat.make_params(a=0.5, b=2.0, c=1.0)
p_rat['c'].min = -2
p_rat['c'].max = 20
r_rat = m_rat.fit(hybrid_zm, x=hybrid_q, params=p_rat, weights=1/hybrid_zm_err)
ss_res = np.sum(r_rat.residual**2 * hybrid_zm_err**2)
r2_rat = 1 - ss_res/ss_tot if ss_tot > 0 else 0
a_r, b_r, c_r = r_rat.params['a'].value, r_rat.params['b'].value, r_rat.params['c'].value
q_cross_rat = b_r / (1.0 - a_r) - c_r if abs(1.0 - a_r) > 1e-10 else None
print(f"\nRational: z_m = {a_r:.4f} + {b_r:.4f}/(q+{c_r:.3f})")
print(f"  R² = {r2_rat:.5f}, AIC = {r_rat.aic:.2f}")
print(f"  q_cross = {q_cross_rat:.3f}" if q_cross_rat else "  q_cross = undefined")
fit_results['rational'] = {'a': a_r, 'b': b_r, 'c': c_r,
                            'r2': r2_rat, 'aic': r_rat.aic, 'q_cross': q_cross_rat}

# ---- BEST FIT SELECTION ----
print("\n--- FIT COMPARISON ---\n")
print(f"{'Model':>10} {'R²':>8} {'AIC':>8} {'q_cross':>8}")
print("-" * 40)
best_aic = min(f['aic'] for f in fit_results.values())
for name, f in fit_results.items():
    daic = f['aic'] - best_aic
    marker = " <-- BEST" if daic < 0.01 else f" (dAIC={daic:.1f})"
    qc = f"{f['q_cross']:.3f}" if f['q_cross'] is not None else "—"
    print(f"{name:>10} {f['r2']:>8.5f} {f['aic']:>8.2f} {qc:>8}{marker}")

# Collect q_cross estimates
q_crosses = [f['q_cross'] for f in fit_results.values() if f['q_cross'] is not None]
q_cross_mean = np.mean(q_crosses)
q_cross_std = np.std(q_crosses)
print(f"\nq_cross consensus: {q_cross_mean:.3f} ± {q_cross_std:.3f}")
print(f"  range: [{min(q_crosses):.3f}, {max(q_crosses):.3f}]")

results['data']['hybrid_fits'] = fit_results
results['data']['q_cross_consensus'] = {'mean': q_cross_mean, 'std': q_cross_std,
                                         'range': [min(q_crosses), max(q_crosses)]}

# ---- SIMPLE INTERPOLATION ----
print("\n--- LINEAR INTERPOLATION between q=3 and q=5 ---\n")
# z_m crosses 1 between q=3 (1.030) and q=4 (1.004) and q=5 (0.909)
# Linear interpolation between q=3 and q=4:
q_cross_34 = 3 + (1.0 - 1.030) / (1.004 - 1.030)  # where z_m=1 between q=3,4
q_cross_45 = 4 + (1.0 - 1.004) / (0.909 - 1.004)  # where z_m=1 between q=4,5
print(f"Linear interp q=3..4: q_cross = {q_cross_34:.3f}")
print(f"Linear interp q=4..5: q_cross = {q_cross_45:.3f}")
print(f"Average: q_cross ~ {(q_cross_34 + q_cross_45)/2:.3f}")
results['data']['interpolation'] = {'q_cross_34': q_cross_34, 'q_cross_45': q_cross_45}

# ---- S_Q z_m(q) FIT ----
print("\n--- S_Q POTTS z_m(q) fits ---\n")
ss_tot_sq = np.sum((sq_zm - sq_zm.mean())**2)

# Linear
r_sq_lin = m_lin.fit(sq_zm, x=sq_q, params=m_lin.make_params(a=0.7, b=0.1), weights=1/sq_zm_err)
ss_res_sq = np.sum(r_sq_lin.residual**2 * sq_zm_err**2)
r2_sq_lin = 1 - ss_res_sq/ss_tot_sq if ss_tot_sq > 0 else 0
print(f"Linear: z_m = {r_sq_lin.params['a'].value:.4f} + {r_sq_lin.params['b'].value:.4f}*q")
print(f"  R² = {r2_sq_lin:.5f}")

# Log
r_sq_log = m_log.fit(sq_zm, x=sq_q, params=m_log.make_params(a=0.5, b=0.4), weights=1/sq_zm_err)
ss_res_sq = np.sum(r_sq_log.residual**2 * sq_zm_err**2)
r2_sq_log = 1 - ss_res_sq/ss_tot_sq if ss_tot_sq > 0 else 0
print(f"Log: z_m = {r_sq_log.params['a'].value:.4f} + {r_sq_log.params['b'].value:.4f}*ln(q)")
print(f"  R² = {r2_sq_log:.5f}")

results['data']['sq_fits'] = {
    'linear': {'a': r_sq_lin.params['a'].value, 'b': r_sq_lin.params['b'].value, 'r2': r2_sq_lin},
    'log': {'a': r_sq_log.params['a'].value, 'b': r_sq_log.params['b'].value, 'r2': r2_sq_log},
}

# ---- CONVERGENCE BEHAVIOR ----
print("\n--- z_m BEHAVIOR COMPARISON ---\n")
print(f"{'q':>3} {'hybrid z_m':>10} {'S_q z_m':>10} {'divergence':>10}")
for i, q in enumerate([3, 4, 5, 7]):
    div = abs(hybrid_zm[i] - sq_zm[i]) / sq_zm[i] * 100
    print(f"{q:>3} {hybrid_zm[i]:>10.3f} {sq_zm[i]:>10.3f} {div:>9.1f}%")
print(f"\nHybrid: z_m DECREASING -- crosses 1 at q~{q_cross_mean:.1f}")
print(f"S_q:    z_m INCREASING -- stays above 1 (walking persists)")

# ---- RECORD KEY QUANTITIES ----
print("\n--- RECORDING TO DB ---\n")
record(sprint=121, model='hybrid', q=0, n=0, quantity='q_cross_zm',
       value=q_cross_mean, error=q_cross_std, method='zm_q_fit',
       notes=f'q where z_m=1, consensus of 4 fits')
print(f"Recorded: q_cross = {q_cross_mean:.3f} ± {q_cross_std:.3f}")

# Record hybrid z_m values
for i, q in enumerate(hybrid_q):
    record(sprint=121, model='hybrid', q=int(q), n=0, quantity='z_m',
           value=float(hybrid_zm[i]), error=float(hybrid_zm_err[i]),
           method='chif_spectral', notes='from sprints 119-120')

save()
print("\nResults saved to results/sprint_121a_zm_q_crossing.json")
print("\nDone.")
