"""Sprint 110c: Refit alpha(q) with all available data including q=6 (4 sizes), q=8, q=9.

Tests: (1) linear vs quadratic fit, (2) effect of using last-pair vs global alpha,
(3) updated z_m(q) and beta_me(q) linear fits.
"""
import numpy as np
import json, time, os

results = {
    'experiment': '110c_alpha_refit',
    'sprint': 110,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_110c_alpha_refit.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def load_json(path):
    with open(path) as f:
        return json.load(f)

# Collect per-size chi_F, gap_m, me_sq for each q in walking regime

# q=5: Sprint 107a (n=6,7,8,9) + 107a2 (n=10)
d107a = load_json('results/sprint_107a_extend_spectral.json')
d107a2 = load_json('results/sprint_107a2_gpu_extend.json')
q5_data = {}
for nk, nd in d107a['data']['q5']['sizes'].items():
    q5_data[int(nk)] = nd
# Add n=10 from 107a2
nd10 = d107a2['data']['q5_n10']
q5_data[10] = nd10

# q=6: Sprint 110a (n=6,7,8,9)
d110a = load_json('results/sprint_110a_q6_extend.json')
q6_data = {int(k): v for k, v in d110a['data'].items()}

# q=7: Sprint 106b (n=6,8,9) + 106a (n=6) + 107a2 (n=8)
# Best approach: use 106b for n=6,8,9 and 107a2 for n=8 (most sizes)
d106b = load_json('results/sprint_106b_scaling_decomp.json')
q7_data = {}
for nk, nd in d106b['data']['q7']['sizes'].items():
    q7_data[int(nk)] = nd
# Add n=8 from 107a2 if not already present
if 8 not in q7_data:
    q7_data[8] = d107a2['data']['q7_n8']
else:
    # Update with 107a2 data (more recent, potentially more accurate)
    q7_data[8] = d107a2['data']['q7_n8']

# q=8, q=9: Sprint 110b
d110b = load_json('results/sprint_110b_q8q9_new.json')
q8_data, q9_data = {}, {}
for key, v in d110b['data'].items():
    if v['q'] == 8:
        q8_data[v['n']] = v
    elif v['q'] == 9:
        q9_data[v['n']] = v

# Compute pairwise exponents
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

    # Global fit
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

all_results = {}
for q, qd in [(5, q5_data), (6, q6_data), (7, q7_data), (8, q8_data), (9, q9_data)]:
    all_results[q] = compute_pairwise(qd)

print("=" * 70)
print("ALPHA(q) COMPREHENSIVE DATA TABLE")
print("=" * 70)
print(f"{'q':>3} {'n_sizes':>20} {'alpha_last':>12} {'alpha_glob':>12} {'z_m_last':>10} {'beta_last':>10} {'old_formula':>12}")
for q in sorted(all_results.keys()):
    d = all_results[q]
    formula = 0.315 * q + 0.469
    print(f"{q:>3} {str(d['n_sizes']):>20} {d['alpha_last']:>12.4f} {d['alpha_global']:>12.4f} "
          f"{d['z_m_last']:>10.4f} {d['beta_last']:>10.4f} {formula:>12.3f}")

# Fit alpha(q)
print(f"\n{'=' * 70}")
print("ALPHA(q) FITS — WALKING REGIME (q>=5)")
print("=" * 70)

qs = np.array([5, 6, 7, 8, 9], dtype=float)
alphas_last = np.array([all_results[q]['alpha_last'] for q in [5, 6, 7, 8, 9]])
alphas_glob = np.array([all_results[q]['alpha_global'] for q in [5, 6, 7, 8, 9]])

# Linear fit
p_lin = np.polyfit(qs, alphas_last, 1)
res_lin = alphas_last - np.polyval(p_lin, qs)
rms_lin = np.sqrt(np.mean(res_lin**2))
print(f"\n  Linear (last-pair): alpha = {p_lin[0]:.4f}*q + {p_lin[1]:.4f}")
print(f"    Residuals: {['%+.4f' % r for r in res_lin]}")
print(f"    RMS: {rms_lin:.5f}")

# Quadratic fit
p_quad = np.polyfit(qs, alphas_last, 2)
res_quad = alphas_last - np.polyval(p_quad, qs)
rms_quad = np.sqrt(np.mean(res_quad**2))
print(f"\n  Quadratic (last-pair): alpha = {p_quad[0]:.5f}*q² + {p_quad[1]:.4f}*q + {p_quad[2]:.4f}")
print(f"    Residuals: {['%+.4f' % r for r in res_quad]}")
print(f"    RMS: {rms_quad:.5f}")

# AIC comparison (3 vs 2 params, 5 data points)
n_pts = len(qs)
if rms_lin > 0 and rms_quad > 0:
    aic_lin = n_pts * np.log(rms_lin**2) + 2 * 2
    aic_quad = n_pts * np.log(rms_quad**2) + 2 * 3
    print(f"\n  AIC linear: {aic_lin:.2f}")
    print(f"  AIC quadratic: {aic_quad:.2f}")
    print(f"  {'Linear' if aic_lin < aic_quad else 'Quadratic'} preferred (ΔAIC={abs(aic_lin-aic_quad):.2f})")

# Linear fit on global alpha
p_lin_g = np.polyfit(qs, alphas_glob, 1)
print(f"\n  Linear (global): alpha = {p_lin_g[0]:.4f}*q + {p_lin_g[1]:.4f}")

# Compare old formula
print(f"\n  Old formula (Sprint 103): alpha = 0.315*q + 0.469")
print(f"  New formula (last-pair):  alpha = {p_lin[0]:.4f}*q + {p_lin[1]:.4f}")
print(f"  New formula (global):     alpha = {p_lin_g[0]:.4f}*q + {p_lin_g[1]:.4f}")

# z_m(q) and beta_me(q)
print(f"\n{'=' * 70}")
print("COMPONENT FITS: z_m(q) AND beta_me(q)")
print("=" * 70)

z_ms = np.array([all_results[q]['z_m_last'] for q in [5, 6, 7, 8, 9]])
betas = np.array([all_results[q]['beta_last'] for q in [5, 6, 7, 8, 9]])

p_z = np.polyfit(qs, z_ms, 1)
p_b = np.polyfit(qs, betas, 1)

print(f"\n  z_m(q) = {p_z[0]:.4f}*q + {p_z[1]:.4f}")
print(f"    Old (Sprint 107): 0.065*q + 0.845")
print(f"    Values: {['%.4f' % z for z in z_ms]}")
res_z = z_ms - np.polyval(p_z, qs)
print(f"    Residuals: {['%+.5f' % r for r in res_z]}, RMS={np.sqrt(np.mean(res_z**2)):.5f}")

print(f"\n  beta_me(q) = {p_b[0]:.4f}*q + {p_b[1]:.4f}")
print(f"    Old (Sprint 107): 0.188*q - 0.238")
print(f"    Values: {['%.4f' % b for b in betas]}")
res_b = betas - np.polyval(p_b, qs)
print(f"    Residuals: {['%+.5f' % r for r in res_b]}, RMS={np.sqrt(np.mean(res_b**2)):.5f}")

# Reconstructed
slope_recon = p_b[0] + 2 * p_z[0]
intercept_recon = p_b[1] + 2 * p_z[1] - 1
print(f"\n  Reconstructed: alpha = {slope_recon:.4f}*q + {intercept_recon:.4f}")
print(f"  Direct fit:    alpha = {p_lin[0]:.4f}*q + {p_lin[1]:.4f}")
print(f"  Difference: slope {abs(slope_recon - p_lin[0]):.4f}, intercept {abs(intercept_recon - p_lin[1]):.4f}")

# Finite-size convergence
print(f"\n{'=' * 70}")
print("FINITE-SIZE CONVERGENCE")
print("=" * 70)
for q in [5, 6, 7, 8, 9]:
    d = all_results[q]
    pw = d['pw_alpha']
    n_sizes = d['n_sizes']
    pairs = [f"({n_sizes[i]},{n_sizes[i+1]})" for i in range(len(pw))]
    print(f"  q={q}: {' '.join(f'{p}={a:.4f}' for p, a in zip(pairs, pw))}")
    if len(pw) >= 2:
        drift = pw[-1] - pw[-2]
        print(f"        Last drift: {drift:+.4f} ({'UP' if drift > 0 else 'DOWN'})")

# Extrapolation: if q=8,9 are converging upward like q=5,6,7
print(f"\n{'=' * 70}")
print("BIAS CORRECTION ESTIMATE")
print("=" * 70)
# q=5,6 show upward convergence. Estimate bias from n=6,7 pair vs converged value
# q=5: pair(6,7)=2.076, converged=2.100 → bias at (6,7) = -0.024
# q=6: pair(6,7)=2.358, last pair(8,9)=2.402 → bias at (6,7) ≈ -0.044
# q=7: pair(6,7)=2.630, pair(7,8)=2.670 → bias at (6,7) ≈ -0.040
q5_67_pair = all_results[5]['pw_alpha'][0]  # (6,7) pair
q5_last = all_results[5]['alpha_last']
q6_67_pair = all_results[6]['pw_alpha'][0]
q6_last = all_results[6]['alpha_last']
q7_67_pair = all_results[7]['pw_alpha'][0]
q7_last = all_results[7]['alpha_last']

bias_5 = q5_67_pair - q5_last
bias_6 = q6_67_pair - q6_last
bias_7 = q7_67_pair - q7_last

print(f"  q=5: (6,7) pair = {q5_67_pair:.4f}, converged = {q5_last:.4f}, bias = {bias_5:+.4f}")
print(f"  q=6: (6,7) pair = {q6_67_pair:.4f}, last pair = {q6_last:.4f}, bias = {bias_6:+.4f}")
print(f"  q=7: (6,7) pair = {q7_67_pair:.4f}, last pair = {q7_last:.4f}, bias = {bias_7:+.4f}")

avg_bias = np.mean([bias_5, bias_6, bias_7])
print(f"\n  Average (6,7)-pair bias: {avg_bias:+.4f}")

# Corrected q=8,9 estimates
q8_corrected = all_results[8]['alpha_last'] - avg_bias
q9_corrected = all_results[9]['alpha_last'] - avg_bias
print(f"\n  q=8: measured {all_results[8]['alpha_last']:.4f}, bias-corrected → {q8_corrected:.4f}")
print(f"  q=9: measured {all_results[9]['alpha_last']:.4f}, bias-corrected → {q9_corrected:.4f}")

# Refit with corrected values
alphas_corrected = alphas_last.copy()
alphas_corrected[3] = q8_corrected  # q=8
alphas_corrected[4] = q9_corrected  # q=9
p_corr = np.polyfit(qs, alphas_corrected, 1)
res_corr = alphas_corrected - np.polyval(p_corr, qs)
rms_corr = np.sqrt(np.mean(res_corr**2))
print(f"\n  Linear fit (bias-corrected): alpha = {p_corr[0]:.4f}*q + {p_corr[1]:.4f}")
print(f"    RMS: {rms_corr:.5f}")

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

results['fits'] = {
    'linear_last': {'slope': float(p_lin[0]), 'intercept': float(p_lin[1]), 'rms': float(rms_lin)},
    'quadratic_last': {'a': float(p_quad[0]), 'b': float(p_quad[1]), 'c': float(p_quad[2]), 'rms': float(rms_quad)},
    'linear_global': {'slope': float(p_lin_g[0]), 'intercept': float(p_lin_g[1])},
    'z_m': {'slope': float(p_z[0]), 'intercept': float(p_z[1])},
    'beta_me': {'slope': float(p_b[0]), 'intercept': float(p_b[1])},
    'linear_bias_corrected': {'slope': float(p_corr[0]), 'intercept': float(p_corr[1]), 'rms': float(rms_corr)},
    'bias_avg': float(avg_bias),
}
save()
print("\nResults saved.")
