"""Sprint 121a: Fit z_m(q) for hybrid model to find q_cross where z_m=1.

Uses data from Sprints 119-120. Five data points: q=3,4,5,7,10.
Tests: linear, logarithmic, power-law, and rational forms.
Finds q_cross = q where z_m(q) = 1 for each model.
"""
import numpy as np
import json
from lmfit import Model

# Hybrid z_m data from Sprints 119-120 (global fits)
q_vals = np.array([3, 4, 5, 7, 10], dtype=float)
zm_hybrid = np.array([1.03, 1.004, 0.91, 0.77, 0.69])
# S_q z_m data (from Sprints 116-118, 119)
zm_sq = np.array([1.03, 1.05, 1.31, 1.4, np.nan])  # q=10 not measured for S_q

# Hybrid alpha data
alpha_hybrid = np.array([1.40, 1.55, 1.41, 0.95, 0.55])
alpha_sq = np.array([1.40, 1.77, 2.09, 2.65, 3.2])

# --- Fit models for z_m(q) ---

def linear(q, a, b):
    return a + b * q

def logarithmic(q, a, b):
    return a + b * np.log(q)

def power_law(q, a, b, c):
    return a + b * q**c

def rational(q, a, b):
    return a / (1 + b * q)

def inv_sqrt(q, a, b):
    return a + b / np.sqrt(q)

results = {}

# Fit each model to hybrid z_m data
print("=" * 60)
print("HYBRID z_m(q) FITTING")
print("=" * 60)
print(f"\nData: q = {q_vals.tolist()}")
print(f"       z_m = {zm_hybrid.tolist()}")

fits = {}

# Linear
m = Model(linear)
p = m.make_params(a=1.3, b=-0.05)
r = m.fit(zm_hybrid, p, q=q_vals)
ss_res = np.sum(r.residual**2)
ss_tot = np.sum((zm_hybrid - zm_hybrid.mean())**2)
r2 = 1 - ss_res / ss_tot
aic = r.aic
q_cross_lin = (1 - r.params['a'].value) / r.params['b'].value if r.params['b'].value != 0 else np.nan
fits['linear'] = {'r2': r2, 'aic': aic, 'q_cross': q_cross_lin,
                  'params': {k: v.value for k, v in r.params.items()}}
print(f"\nLinear: z_m = {r.params['a'].value:.4f} + {r.params['b'].value:.4f}*q")
print(f"  R² = {r2:.6f}, AIC = {aic:.2f}, q_cross = {q_cross_lin:.3f}")

# Logarithmic
m = Model(logarithmic)
p = m.make_params(a=1.8, b=-0.4)
r = m.fit(zm_hybrid, p, q=q_vals)
ss_res = np.sum(r.residual**2)
ss_tot = np.sum((zm_hybrid - zm_hybrid.mean())**2)
r2 = 1 - ss_res / ss_tot
aic = r.aic
q_cross_log = np.exp((1 - r.params['a'].value) / r.params['b'].value) if r.params['b'].value != 0 else np.nan
fits['logarithmic'] = {'r2': r2, 'aic': aic, 'q_cross': q_cross_log,
                       'params': {k: v.value for k, v in r.params.items()}}
print(f"\nLogarithmic: z_m = {r.params['a'].value:.4f} + {r.params['b'].value:.4f}*ln(q)")
print(f"  R² = {r2:.6f}, AIC = {aic:.2f}, q_cross = {q_cross_log:.3f}")

# 1/sqrt(q)
m = Model(inv_sqrt)
p = m.make_params(a=0.3, b=0.4)
r = m.fit(zm_hybrid, p, q=q_vals)
ss_res = np.sum(r.residual**2)
ss_tot = np.sum((zm_hybrid - zm_hybrid.mean())**2)
r2 = 1 - ss_res / ss_tot
aic = r.aic
# Solve a + b/sqrt(q) = 1 => sqrt(q) = b/(1-a) => q = (b/(1-a))^2
a_v, b_v = r.params['a'].value, r.params['b'].value
q_cross_isq = (b_v / (1 - a_v))**2 if (1 - a_v) != 0 else np.nan
fits['inv_sqrt'] = {'r2': r2, 'aic': aic, 'q_cross': q_cross_isq,
                    'params': {k: v.value for k, v in r.params.items()}}
print(f"\nInv sqrt: z_m = {a_v:.4f} + {b_v:.4f}/sqrt(q)")
print(f"  R² = {r2:.6f}, AIC = {aic:.2f}, q_cross = {q_cross_isq:.3f}")

# Rational: a / (1 + b*q)
m = Model(rational)
p = m.make_params(a=1.5, b=0.1)
r = m.fit(zm_hybrid, p, q=q_vals)
ss_res = np.sum(r.residual**2)
ss_tot = np.sum((zm_hybrid - zm_hybrid.mean())**2)
r2 = 1 - ss_res / ss_tot
aic = r.aic
a_v, b_v = r.params['a'].value, r.params['b'].value
q_cross_rat = (a_v - 1) / b_v if b_v != 0 else np.nan
fits['rational'] = {'r2': r2, 'aic': aic, 'q_cross': q_cross_rat,
                    'params': {k: v.value for k, v in r.params.items()}}
print(f"\nRational: z_m = {a_v:.4f} / (1 + {b_v:.4f}*q)")
print(f"  R² = {r2:.6f}, AIC = {aic:.2f}, q_cross = {q_cross_rat:.3f}")

# Power-law: a + b * q^c
m = Model(power_law)
p = m.make_params(a=1.5, b=-0.3, c=0.5)
p['c'].min = -3
p['c'].max = 3
r = m.fit(zm_hybrid, p, q=q_vals)
ss_res = np.sum(r.residual**2)
ss_tot = np.sum((zm_hybrid - zm_hybrid.mean())**2)
r2 = 1 - ss_res / ss_tot
aic = r.aic
a_v, b_v, c_v = r.params['a'].value, r.params['b'].value, r.params['c'].value
# Numerical solve for q_cross
from scipy.optimize import brentq
try:
    q_cross_pl = brentq(lambda q: a_v + b_v * q**c_v - 1, 2, 20)
except:
    q_cross_pl = np.nan
fits['power_law'] = {'r2': r2, 'aic': aic, 'q_cross': q_cross_pl,
                     'params': {k: v.value for k, v in r.params.items()}}
print(f"\nPower-law: z_m = {a_v:.4f} + {b_v:.4f}*q^{c_v:.4f}")
print(f"  R² = {r2:.6f}, AIC = {aic:.2f}, q_cross = {q_cross_pl:.3f}")

# Summary
print("\n" + "=" * 60)
print("SUMMARY: q_cross estimates")
print("=" * 60)
q_crosses = []
for name, f in sorted(fits.items(), key=lambda x: -x[1]['r2']):
    print(f"  {name:15s}: q_cross = {f['q_cross']:.3f}, R² = {f['r2']:.6f}, AIC = {f['aic']:.2f}")
    if not np.isnan(f['q_cross']):
        q_crosses.append(f['q_cross'])

q_cross_mean = np.mean(q_crosses)
q_cross_std = np.std(q_crosses)
print(f"\n  Mean q_cross = {q_cross_mean:.3f} ± {q_cross_std:.3f}")

# Also fit hybrid alpha(q) to characterize the non-monotonic peak
print("\n" + "=" * 60)
print("HYBRID alpha(q) — characterizing the peak")
print("=" * 60)
print(f"q = {q_vals.tolist()}")
print(f"alpha = {alpha_hybrid.tolist()}")

# Quadratic in log(q)
def quad_logq(q, a, b, c):
    lq = np.log(q)
    return a + b * lq + c * lq**2

m = Model(quad_logq)
p = m.make_params(a=-1, b=3, c=-1)
r = m.fit(alpha_hybrid, p, q=q_vals)
ss_res = np.sum(r.residual**2)
ss_tot = np.sum((alpha_hybrid - alpha_hybrid.mean())**2)
r2 = 1 - ss_res / ss_tot
a_v, b_v, c_v = r.params['a'].value, r.params['b'].value, r.params['c'].value
q_peak = np.exp(-b_v / (2 * c_v)) if c_v != 0 else np.nan
alpha_peak = a_v + b_v * np.log(q_peak) + c_v * np.log(q_peak)**2
print(f"\nalpha = {a_v:.4f} + {b_v:.4f}*ln(q) + {c_v:.4f}*ln(q)²")
print(f"  R² = {r2:.6f}")
print(f"  Peak at q = {q_peak:.3f}, alpha_peak = {alpha_peak:.3f}")

# Compare z_m and alpha at the crossing
print("\n" + "=" * 60)
print("z_m vs alpha RELATIONSHIP")
print("=" * 60)
print("q    z_m    alpha   alpha_predicted(=beta_me+2*z_m-1)")
for i, q in enumerate(q_vals):
    # From decomposition: alpha = beta_me + 2*z_m - 1
    # beta_me is fitted independently. Let's see how z_m correlates with alpha directly
    print(f"{int(q):3d}  {zm_hybrid[i]:.3f}  {alpha_hybrid[i]:.3f}")

# z_m asymptote: what does z_m approach as q->inf?
print(f"\nz_m large-q trend: z_m(7)={zm_hybrid[3]:.3f}, z_m(10)={zm_hybrid[4]:.3f}")
print(f"  If z_m -> 0: the spectral gap closes slower than 1/N")
print(f"  If z_m -> const>0: gap closes as N^(-z_m_inf)")

# Save results
results = {
    'q_vals': q_vals.tolist(),
    'zm_hybrid': zm_hybrid.tolist(),
    'zm_sq': [float(x) if not np.isnan(x) else None for x in zm_sq],
    'alpha_hybrid': alpha_hybrid.tolist(),
    'alpha_sq': alpha_sq.tolist(),
    'fits': {},
    'q_cross_mean': float(q_cross_mean),
    'q_cross_std': float(q_cross_std),
    'alpha_peak_q': float(q_peak),
}
for name, f in fits.items():
    results['fits'][name] = {
        'r2': float(f['r2']),
        'aic': float(f['aic']),
        'q_cross': float(f['q_cross']) if not np.isnan(f['q_cross']) else None,
        'params': {k: float(v) for k, v in f['params'].items()},
    }

with open('../results/sprint_121a_zm_q_fit.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to results/sprint_121a_zm_q_fit.json")
