"""Sprint 111b: Refit alpha(q) with q=10 data point included.

Tests linear vs quadratic vs power-law discrimination with 6 walking q values (q=5-10).
Sprint 110 had 5 values — q=10 is the critical test point.
Linear predicts 3.44, quadratic predicts 3.38.
"""
import numpy as np
import json, time, os

results = {
    'experiment': '111b_alpha_refit_q10',
    'sprint': 111,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_111b_alpha_refit_q10.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def load_json(path):
    with open(path) as f:
        return json.load(f)

# ---- Collect all per-size data ----

# q=5: Sprint 107a (n=6,7,8,9) + 107a2 (n=10)
d107a = load_json('results/sprint_107a_extend_spectral.json')
d107a2 = load_json('results/sprint_107a2_gpu_extend.json')
q5_data = {}
for nk, nd in d107a['data']['q5']['sizes'].items():
    q5_data[int(nk)] = nd
q5_data[10] = d107a2['data']['q5_n10']

# q=6: Sprint 110a (n=6,7,8,9)
d110a = load_json('results/sprint_110a_q6_extend.json')
q6_data = {int(k): v for k, v in d110a['data'].items()}

# q=7: Sprint 106b + 107a2
d106b = load_json('results/sprint_106b_scaling_decomp.json')
q7_data = {}
for nk, nd in d106b['data']['q7']['sizes'].items():
    q7_data[int(nk)] = nd
q7_data[8] = d107a2['data']['q7_n8']

# q=8, q=9: Sprint 110b
d110b = load_json('results/sprint_110b_q8q9_new.json')
q8_data, q9_data = {}, {}
for key, v in d110b['data'].items():
    if v['q'] == 8:
        q8_data[v['n']] = v
    elif v['q'] == 9:
        q9_data[v['n']] = v

# q=10: Sprint 111a (NEW)
d111a = load_json('results/sprint_111a_q10_chif.json')
q10_data = {}
for key, v in d111a['data'].items():
    q10_data[v['n']] = v

# ---- Compute pairwise exponents ----
def compute_pairwise(q_data_dict):
    ns = sorted(q_data_dict.keys())
    chis = [q_data_dict[n]['chi_F'] for n in ns]
    gaps = [q_data_dict[n]['gap_multiplet'] for n in ns]
    mes = [q_data_dict[n]['me_sq'] for n in ns]

    ns_arr = np.array(ns, dtype=float)
    log_N = np.log(ns_arr)
    log_chi = np.log(chis)
    log_gap = np.log(gaps)
    log_me = np.log(mes)

    pw_alpha, pw_z, pw_beta = [], [], []
    for i in range(len(ns) - 1):
        dln = log_N[i+1] - log_N[i]
        pw_alpha.append((log_chi[i+1] - log_chi[i]) / dln)
        pw_z.append(-(log_gap[i+1] - log_gap[i]) / dln)
        pw_beta.append((log_me[i+1] - log_me[i]) / dln)

    p = np.polyfit(log_N, log_chi, 1)

    return {
        'n_sizes': ns,
        'pw_alpha': pw_alpha,
        'pw_z': pw_z,
        'pw_beta': pw_beta,
        'alpha_last': pw_alpha[-1],
        'alpha_global': p[0],
        'z_m_last': pw_z[-1],
        'beta_last': pw_beta[-1],
    }

all_q_data = {5: q5_data, 6: q6_data, 7: q7_data, 8: q8_data, 9: q9_data, 10: q10_data}
all_results = {}
for q, qd in all_q_data.items():
    all_results[q] = compute_pairwise(qd)

# ---- Display ----
print("=" * 70)
print("ALPHA(q) DATA TABLE (q=5-10)")
print("=" * 70)
print(f"{'q':>3} {'n_sizes':>24} {'alpha_last':>12} {'alpha_glob':>12} {'z_m_last':>10} {'beta_last':>10}")
for q in sorted(all_results.keys()):
    d = all_results[q]
    print(f"{q:>3} {str(d['n_sizes']):>24} {d['alpha_last']:>12.4f} {d['alpha_global']:>12.4f} "
          f"{d['z_m_last']:>10.4f} {d['beta_last']:>10.4f}")

# ---- Fits with q=5-10 (6 points) ----
print(f"\n{'=' * 70}")
print("ALPHA(q) FITS — q=5-10 (6 data points)")
print("=" * 70)

qs = np.array([5, 6, 7, 8, 9, 10], dtype=float)
alphas_last = np.array([all_results[q]['alpha_last'] for q in [5, 6, 7, 8, 9, 10]])
alphas_glob = np.array([all_results[q]['alpha_global'] for q in [5, 6, 7, 8, 9, 10]])

n_pts = len(qs)

# 1. Linear fit (last-pair)
p_lin = np.polyfit(qs, alphas_last, 1)
res_lin = alphas_last - np.polyval(p_lin, qs)
rms_lin = np.sqrt(np.mean(res_lin**2))
ss_lin = np.sum(res_lin**2)
print(f"\n  1. Linear: alpha = {p_lin[0]:.4f}*q + {p_lin[1]:.4f}")
print(f"     Residuals: {['%+.4f' % r for r in res_lin]}")
print(f"     RMS: {rms_lin:.5f}")

# 2. Quadratic fit (last-pair)
p_quad = np.polyfit(qs, alphas_last, 2)
res_quad = alphas_last - np.polyval(p_quad, qs)
rms_quad = np.sqrt(np.mean(res_quad**2))
ss_quad = np.sum(res_quad**2)
print(f"\n  2. Quadratic: alpha = {p_quad[0]:.5f}*q² + {p_quad[1]:.4f}*q + {p_quad[2]:.4f}")
print(f"     Residuals: {['%+.4f' % r for r in res_quad]}")
print(f"     RMS: {rms_quad:.5f}")

# 3. Power-law: alpha = a * q^b + c (fit log(alpha) vs log(q) won't work if intercept)
# Try alpha = a * (q - q0)^b — too many params
# Simpler: alpha = a * q^b
from scipy.optimize import curve_fit
def power_law(q, a, b):
    return a * q**b

try:
    popt, pcov = curve_fit(power_law, qs, alphas_last, p0=[0.3, 1.0])
    res_pow = alphas_last - power_law(qs, *popt)
    rms_pow = np.sqrt(np.mean(res_pow**2))
    ss_pow = np.sum(res_pow**2)
    print(f"\n  3. Power-law: alpha = {popt[0]:.4f} * q^{popt[1]:.4f}")
    print(f"     Residuals: {['%+.4f' % r for r in res_pow]}")
    print(f"     RMS: {rms_pow:.5f}")
except Exception as e:
    print(f"\n  3. Power-law fit failed: {e}")
    rms_pow = None
    ss_pow = None

# 4. sqrt fit: alpha = a * sqrt(q) + b
def sqrt_fit(q, a, b):
    return a * np.sqrt(q) + b
try:
    popt_sq, pcov_sq = curve_fit(sqrt_fit, qs, alphas_last, p0=[1.0, 0.0])
    res_sq = alphas_last - sqrt_fit(qs, *popt_sq)
    rms_sq = np.sqrt(np.mean(res_sq**2))
    ss_sq = np.sum(res_sq**2)
    print(f"\n  4. Sqrt: alpha = {popt_sq[0]:.4f} * sqrt(q) + {popt_sq[1]:.4f}")
    print(f"     Residuals: {['%+.4f' % r for r in res_sq]}")
    print(f"     RMS: {rms_sq:.5f}")
except Exception as e:
    print(f"\n  4. Sqrt fit failed: {e}")
    rms_sq = None
    ss_sq = None

# AIC comparison
print(f"\n{'=' * 70}")
print("AIC MODEL COMPARISON")
print("=" * 70)
models = []
if rms_lin and rms_lin > 0:
    aic_lin = n_pts * np.log(ss_lin / n_pts) + 2 * 2
    models.append(('Linear', aic_lin, rms_lin, 2))
if rms_quad and rms_quad > 0:
    aic_quad = n_pts * np.log(ss_quad / n_pts) + 2 * 3
    models.append(('Quadratic', aic_quad, rms_quad, 3))
if rms_pow and rms_pow > 0:
    aic_pow = n_pts * np.log(ss_pow / n_pts) + 2 * 2
    models.append(('Power-law', aic_pow, rms_pow, 2))
if rms_sq and rms_sq > 0:
    aic_sq = n_pts * np.log(ss_sq / n_pts) + 2 * 2
    models.append(('Sqrt', aic_sq, rms_sq, 2))

models.sort(key=lambda x: x[1])
best_aic = models[0][1]
print(f"  {'Model':<15} {'AIC':>8} {'ΔAIC':>8} {'RMS':>10} {'params':>6}")
for name, aic, rms, npar in models:
    print(f"  {name:<15} {aic:>8.2f} {aic-best_aic:>8.2f} {rms:>10.5f} {npar:>6}")

# ---- Component fits z_m(q) and beta_me(q) ----
print(f"\n{'=' * 70}")
print("COMPONENT FITS: z_m(q) AND beta_me(q) — q=5-10")
print("=" * 70)

z_ms = np.array([all_results[q]['z_m_last'] for q in [5, 6, 7, 8, 9, 10]])
betas = np.array([all_results[q]['beta_last'] for q in [5, 6, 7, 8, 9, 10]])

p_z = np.polyfit(qs, z_ms, 1)
p_b = np.polyfit(qs, betas, 1)

print(f"\n  z_m(q) = {p_z[0]:.4f}*q + {p_z[1]:.4f}")
print(f"    Sprint 110: 0.082*q + 0.741")
print(f"    Values: {['%.4f' % z for z in z_ms]}")
res_z = z_ms - np.polyval(p_z, qs)
print(f"    Residuals: {['%+.5f' % r for r in res_z]}, RMS={np.sqrt(np.mean(res_z**2)):.5f}")

print(f"\n  beta_me(q) = {p_b[0]:.4f}*q + {p_b[1]:.4f}")
print(f"    Sprint 110: 0.098*q + 0.333")
print(f"    Values: {['%.4f' % b for b in betas]}")
res_b = betas - np.polyval(p_b, qs)
print(f"    Residuals: {['%+.5f' % r for r in res_b]}, RMS={np.sqrt(np.mean(res_b**2)):.5f}")

slope_recon = p_b[0] + 2 * p_z[0]
intercept_recon = p_b[1] + 2 * p_z[1] - 1
print(f"\n  Reconstructed: alpha = {slope_recon:.4f}*q + {intercept_recon:.4f}")
print(f"  Direct linear: alpha = {p_lin[0]:.4f}*q + {p_lin[1]:.4f}")

# ---- Finite-size convergence ----
print(f"\n{'=' * 70}")
print("CONVERGENCE ANALYSIS")
print("=" * 70)
for q in [5, 6, 7, 8, 9, 10]:
    d = all_results[q]
    pw = d['pw_alpha']
    ns = d['n_sizes']
    pairs = [f"({ns[i]},{ns[i+1]})" for i in range(len(pw))]
    print(f"  q={q:>2}: {' '.join(f'{p}={a:.3f}' for p, a in zip(pairs, pw))}")
    if len(pw) >= 2:
        drift = pw[-1] - pw[-2]
        print(f"        Last drift: {drift:+.4f} ({'converging UP' if drift > 0 else 'converging DOWN'})")

# ---- Bias correction for q=8,9,10 ----
print(f"\n{'=' * 70}")
print("BIAS CORRECTION")
print("=" * 70)

# q values with 3+ pairs show upward convergence
# Estimate bias from small-size pairs
# q=5: (9,10)=2.100, global=2.091 — well-converged
# q=6: (8,9)=2.402, global=2.377 — converging
# q=7: (7,8)=2.670, global=2.651 — converging
# q=10: (6,7)=3.417, global=3.349 — still early

# Use q=5 convergence pattern as template
q5_pw = all_results[5]['pw_alpha']
q5_sizes = all_results[5]['n_sizes']
# At q=5, the (6,7) pair was 2.076 and converged to (9,10)=2.100
# Upward shift from (6,7) to last: +0.024
# At q=6, (6,7) pair was 2.358, last (8,9)=2.402 → shift +0.044

# For q=10, only have (5,6)=3.295 and (6,7)=3.417
# Trend is strongly upward (+0.122)
# Estimate: if convergence pattern like q=5-7, true alpha is ~0.02-0.05 above last pair

# Conservative: use last-pair as lower bound, last-pair + 2*drift as upper bound
for q in [8, 9, 10]:
    d = all_results[q]
    pw = d['pw_alpha']
    if len(pw) >= 2:
        drift = pw[-1] - pw[-1]  # no further pairs
        est_lo = pw[-1]
        est_hi = pw[-1] + 0.05  # typical upward bias from q=5-7
        print(f"  q={q}: last pair = {pw[-1]:.4f}, estimated range [{est_lo:.3f}, {est_hi:.3f}]")
    else:
        print(f"  q={q}: only {len(pw)} pair(s), last = {pw[-1]:.4f}")

# ---- Best estimate ----
# Use last-pair for q=5-7 (well-converged), global for q=8,9,10
print(f"\n{'=' * 70}")
print("BEST ESTIMATES")
print("=" * 70)
best_alpha = {}
for q in [5, 6, 7]:
    best_alpha[q] = all_results[q]['alpha_last']
    print(f"  q={q}: alpha = {best_alpha[q]:.4f} (last-pair, well-converged)")
for q in [8, 9, 10]:
    # For q with only 2 sizes, use last pair but note it's a lower bound
    best_alpha[q] = all_results[q]['alpha_last']
    n_pairs = len(all_results[q]['pw_alpha'])
    convergence = "converging UP" if n_pairs >= 2 and all_results[q]['pw_alpha'][-1] > all_results[q]['pw_alpha'][-2] else "2 sizes only"
    print(f"  q={q}: alpha = {best_alpha[q]:.4f} (last-pair, {convergence})")

qs_best = np.array(sorted(best_alpha.keys()), dtype=float)
alphas_best = np.array([best_alpha[q] for q in sorted(best_alpha.keys())])

p_best = np.polyfit(qs_best, alphas_best, 1)
res_best = alphas_best - np.polyval(p_best, qs_best)
rms_best = np.sqrt(np.mean(res_best**2))
print(f"\n  Linear fit (best estimates): alpha = {p_best[0]:.4f}*q + {p_best[1]:.4f}")
print(f"  RMS: {rms_best:.5f}")

p_best_q = np.polyfit(qs_best, alphas_best, 2)
res_best_q = alphas_best - np.polyval(p_best_q, qs_best)
rms_best_q = np.sqrt(np.mean(res_best_q**2))
print(f"  Quadratic fit: alpha = {p_best_q[0]:.5f}*q² + {p_best_q[1]:.4f}*q + {p_best_q[2]:.4f}")
print(f"  RMS: {rms_best_q:.5f}")

# Save everything
results['all_q_data'] = {}
for q, d in all_results.items():
    results['all_q_data'][str(q)] = {
        'n_sizes': d['n_sizes'],
        'alpha_last': float(d['alpha_last']),
        'alpha_global': float(d['alpha_global']),
        'z_m_last': float(d['z_m_last']),
        'beta_last': float(d['beta_last']),
        'pw_alpha': [float(a) for a in d['pw_alpha']],
        'pw_z': [float(z) for z in d['pw_z']],
        'pw_beta': [float(b) for b in d['pw_beta']],
    }

results['fits_6pt'] = {
    'linear': {'slope': float(p_lin[0]), 'intercept': float(p_lin[1]), 'rms': float(rms_lin)},
    'quadratic': {'a': float(p_quad[0]), 'b': float(p_quad[1]), 'c': float(p_quad[2]), 'rms': float(rms_quad)},
    'z_m': {'slope': float(p_z[0]), 'intercept': float(p_z[1])},
    'beta_me': {'slope': float(p_b[0]), 'intercept': float(p_b[1])},
    'best_linear': {'slope': float(p_best[0]), 'intercept': float(p_best[1]), 'rms': float(rms_best)},
    'best_quadratic': {'a': float(p_best_q[0]), 'b': float(p_best_q[1]), 'c': float(p_best_q[2]), 'rms': float(rms_best_q)},
}

save()
print("\nResults saved.")
