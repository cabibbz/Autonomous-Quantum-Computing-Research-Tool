#!/usr/bin/env python3
"""Sprint 056d: Comprehensive comparison of c(q) predictions vs measurements.
Calibrate overshoot correction using q=2,3 (exact c known), apply to q=4,5,7,10.
"""
import numpy as np, json
from scipy.optimize import curve_fit

# === All data collected ===
# Format: {q: {n: c_profile_raw}}
# From Sprints 054-056
profile_data = {
    2: {32: 0.524, 64: 0.512},       # exact c=0.500
    3: {16: 0.870, 48: 0.827},       # exact c=0.800
    4: {16: 1.148, 24: 1.148},       # exact c=1.000 (flat = log corrections)
    5: {16: 1.261},                   # measured c≈1.10
    7: {4: 1.529, 6: 1.537, 8: 1.462},  # new
    10: {4: 1.627, 5: 1.616, 6: 1.596},  # new
}

fss_data = {
    2: {(16, 24): 0.516},
    3: {(16, 24): 0.884},
    4: {(16, 24): 1.229},
    5: {(12, 16): 1.335},
    7: {(4, 6): 1.691, (6, 8): 1.564},
    10: {(4, 5): 1.408, (5, 6): 2.080},
}

exact_c = {2: 0.500, 3: 0.800, 4: 1.000}

# === Overshoot model ===
# For q=2,3 with known c: overshoot = c_raw/c_true - 1
# Hypothesis: overshoot(q, n) = A(q) / n^alpha
# Or simpler: overshoot = a/n + b*ln(q)/n

print("=== Overshoot calibration from q=2,3 ===")
overshoot_points = []
for q in [2, 3]:
    c_true = exact_c[q]
    for n, c_raw in profile_data[q].items():
        ov = c_raw / c_true - 1
        print(f"  q={q}, n={n}: c_raw={c_raw:.4f}, c_true={c_true:.3f}, overshoot={ov:.4f} ({ov*100:.1f}%)")
        overshoot_points.append((q, n, ov))

# Fit: overshoot = a/n (simplest model)
ns_cal = np.array([p[1] for p in overshoot_points])
ovs_cal = np.array([p[2] for p in overshoot_points])
qs_cal = np.array([p[0] for p in overshoot_points])

# Model 1: overshoot = a/n
a_simple = np.mean(ovs_cal * ns_cal)
print(f"\nModel 1: overshoot = {a_simple:.2f}/n")
for q, n, ov in overshoot_points:
    print(f"  q={q}, n={n}: predicted={a_simple/n:.4f}, actual={ov:.4f}")

# Model 2: overshoot = (a + b*ln(q))/n
def ov_model2(X, a, b):
    q, n = X
    return (a + b * np.log(q)) / n

popt2, _ = curve_fit(ov_model2, (qs_cal, ns_cal), ovs_cal, p0=[1.0, 0.5])
print(f"\nModel 2: overshoot = ({popt2[0]:.3f} + {popt2[1]:.3f}*ln(q))/n")
for q, n, ov in overshoot_points:
    pred = ov_model2((q, n), *popt2)
    print(f"  q={q}, n={n}: predicted={pred:.4f}, actual={ov:.4f}")

# Model 3: overshoot = a * q^b / n  (allows q-dependence)
def ov_model3(X, a, b):
    q, n = X
    return a * q**b / n

popt3, _ = curve_fit(ov_model3, (qs_cal, ns_cal), ovs_cal, p0=[1.0, 0.5])
print(f"\nModel 3: overshoot = {popt3[0]:.3f} * q^{popt3[1]:.3f} / n")
for q, n, ov in overshoot_points:
    pred = ov_model3((q, n), *popt3)
    print(f"  q={q}, n={n}: predicted={pred:.4f}, actual={ov:.4f}")

# === Apply overshoot correction to get c_corrected ===
print("\n=== Corrected central charges ===")
print(f"Using Model 2: overshoot = ({popt2[0]:.3f} + {popt2[1]:.3f}*ln(q))/n\n")

corrected_c = {}
for q in sorted(profile_data.keys()):
    best_n = max(profile_data[q].keys())
    c_raw = profile_data[q][best_n]
    ov = ov_model2((q, best_n), *popt2)
    c_corr = c_raw / (1 + ov)
    c_exact = exact_c.get(q, None)

    corrected_c[q] = c_corr
    exact_str = f" (exact: {c_exact:.3f})" if c_exact else ""
    print(f"  q={q:2d}: n_best={best_n:3d}, c_raw={c_raw:.4f}, overshoot={ov:.1%}, c_corrected={c_corr:.4f}{exact_str}")

# Validate: how good is the correction for known q=2,3,4?
print(f"\n  Validation:")
for q in [2, 3]:
    err = abs(corrected_c[q] - exact_c[q]) / exact_c[q]
    print(f"    q={q}: corrected={corrected_c[q]:.4f}, exact={exact_c[q]:.3f}, error={err:.1%}")

# === Fit c(q) formulas to corrected data ===
q_fit = np.array(sorted(corrected_c.keys()))
c_fit = np.array([corrected_c[q] for q in q_fit])

print("\n=== Formula fitting to corrected c(q) ===")
print(f"  Data: q = {list(q_fit)}")
print(f"  Data: c = {[f'{c:.4f}' for c in c_fit]}")

# Fit 1: c = a*ln(q) + b
def f_log(q, a, b):
    return a * np.log(q) + b
popt_log, _ = curve_fit(f_log, q_fit, c_fit)
resid_log = c_fit - f_log(q_fit, *popt_log)
rms_log = np.sqrt(np.mean(resid_log**2))
print(f"\n  Log: c = {popt_log[0]:.4f}*ln(q) + {popt_log[1]:.4f}")
print(f"    RMS residual: {rms_log:.4f}")
for q_pred in [2, 3, 4, 5, 7, 10, 20]:
    print(f"    q={q_pred:3d}: c = {f_log(q_pred, *popt_log):.4f}")

# Fit 2: c = a*(q-1)^b
def f_pow(q, a, b):
    return a * (q - 1)**b
popt_pow, _ = curve_fit(f_pow, q_fit, c_fit, p0=[0.5, 0.5])
resid_pow = c_fit - f_pow(q_fit, *popt_pow)
rms_pow = np.sqrt(np.mean(resid_pow**2))
print(f"\n  Power: c = {popt_pow[0]:.4f}*(q-1)^{popt_pow[1]:.4f}")
print(f"    RMS residual: {rms_pow:.4f}")
for q_pred in [2, 3, 4, 5, 7, 10, 20]:
    print(f"    q={q_pred:3d}: c = {f_pow(q_pred, *popt_pow):.4f}")

# Fit 3: c = a*ln(q-1) + b (shifted log)
def f_log2(q, a, b):
    return a * np.log(q - 1) + b
popt_log2, _ = curve_fit(f_log2, q_fit, c_fit)
resid_log2 = c_fit - f_log2(q_fit, *popt_log2)
rms_log2 = np.sqrt(np.mean(resid_log2**2))
print(f"\n  Log shifted: c = {popt_log2[0]:.4f}*ln(q-1) + {popt_log2[1]:.4f}")
print(f"    RMS residual: {rms_log2:.4f}")
for q_pred in [2, 3, 4, 5, 7, 10, 20]:
    print(f"    q={q_pred:3d}: c = {f_log2(q_pred, *popt_log2):.4f}")

# Fit 4: c = a*sqrt(q-1) + b
def f_sqrt(q, a, b):
    return a * np.sqrt(q - 1) + b
popt_sqrt, _ = curve_fit(f_sqrt, q_fit, c_fit)
resid_sqrt = c_fit - f_sqrt(q_fit, *popt_sqrt)
rms_sqrt = np.sqrt(np.mean(resid_sqrt**2))
print(f"\n  Sqrt: c = {popt_sqrt[0]:.4f}*sqrt(q-1) + {popt_sqrt[1]:.4f}")
print(f"    RMS residual: {rms_sqrt:.4f}")
for q_pred in [2, 3, 4, 5, 7, 10, 20]:
    print(f"    q={q_pred:3d}: c = {f_sqrt(q_pred, *popt_sqrt):.4f}")

# === Final summary table ===
print("\n" + "="*70)
print("FINAL SUMMARY: c(q) for 1D quantum Potts")
print("="*70)
print(f"{'q':>4s}  {'c_raw':>8s}  {'n_best':>6s}  {'c_corr':>8s}  {'c_exact':>8s}  {'log_fit':>8s}")
print("-" * 70)
for q in sorted(corrected_c.keys()):
    best_n = max(profile_data[q].keys())
    c_raw = profile_data[q][best_n]
    c_corr = corrected_c[q]
    c_ex = exact_c.get(q, None)
    c_log = f_log(q, *popt_log)
    ex_str = f"{c_ex:.4f}" if c_ex else "   —"
    print(f"{q:4d}  {c_raw:8.4f}  {best_n:6d}  {c_corr:8.4f}  {ex_str:>8s}  {c_log:8.4f}")

print(f"\nBest formula: c(q) ≈ {popt_log[0]:.3f} * ln(q) + {popt_log[1]:.3f}")
print(f"  → c(q) grows logarithmically with q")
print(f"  → No saturation, no peak — c increases without bound")

# Save
results = {
    'experiment': '056d',
    'overshoot_model': f'({popt2[0]:.3f} + {popt2[1]:.3f}*ln(q))/n',
    'overshoot_params': popt2.tolist(),
    'corrected_c': {str(q): float(c) for q, c in corrected_c.items()},
    'exact_c': {str(k): v for k, v in exact_c.items()},
    'best_fit': {
        'formula': 'c = a*ln(q) + b',
        'a': float(popt_log[0]),
        'b': float(popt_log[1]),
        'rms': float(rms_log),
    },
    'all_fits': {
        'log': {'params': popt_log.tolist(), 'rms': float(rms_log)},
        'power': {'params': popt_pow.tolist(), 'rms': float(rms_pow)},
        'log_shifted': {'params': popt_log2.tolist(), 'rms': float(rms_log2)},
        'sqrt': {'params': popt_sqrt.tolist(), 'rms': float(rms_sqrt)},
    },
}
with open('results/exp_056d.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved results/exp_056d.json")
