#!/usr/bin/env python3
"""Sprint 094c3: Exponential fit to BW breakdown vs nA.

Power law fits poorly (R²=0.42-0.80). Try exponential: 1-R² = A * exp(B*nA).
The Frobenius norm is in q^nA dimensional space — exponential growth expected
from dimensional counting (BW-compatible operators grow linearly in nA,
total operator space grows as q^{2nA}).
"""
import numpy as np
from scipy.optimize import curve_fit
from db_utils import record
import json

# All data from experiments 094a,b,c
data = {
    2: [(3, 2.9294e-04), (4, 5.2547e-04), (5, 1.0004e-02), (6, 2.4014e-01), (7, 3.6235e-01)],
    3: [(3, 5.6425e-04), (4, 8.4161e-04), (5, 2.1924e-01), (6, 4.4226e-01)],
    4: [(3, 8.0548e-04), (4, 1.7538e-03), (5, 3.1168e-01)],
    5: [(3, 1.0243e-03), (4, 3.3679e-02)],
}

results = {'experiment': '094c3_expfit', 'fits': {}}

def exp_model(x, A, B):
    return A * np.exp(B * x)

print("=== Exponential fits: 1-R² = A * exp(B*nA) ===\n")
print(f"{'q':>3} {'B':>8} {'A':>12} {'fit_R²':>8} {'pts':>4}")
print("-" * 40)

for q in sorted(data.keys()):
    pts = data[q]
    nA_arr = np.array([p[0] for p in pts], dtype=float)
    omr_arr = np.array([p[1] for p in pts])

    # Log-linear fit: ln(1-R²) = ln(A) + B*nA
    log_omr = np.log(omr_arr)
    coeffs = np.polyfit(nA_arr, log_omr, 1)
    B = coeffs[0]
    A = np.exp(coeffs[1])

    pred = A * np.exp(B * nA_arr)
    ss_res = np.sum((omr_arr - pred)**2)
    ss_tot = np.sum((omr_arr - np.mean(omr_arr))**2)
    fit_R2 = 1 - ss_res / ss_tot if ss_tot > 1e-30 else 1.0

    print(f"{q:>3} {B:>8.3f} {A:>12.2e} {fit_R2:>8.4f} {len(pts):>4}")
    for nA, omr in pts:
        p = A * np.exp(B * nA)
        print(f"    nA={nA}: measured {omr:.4e}, predicted {p:.4e}, ratio {omr/p:.3f}")

    results['fits'][f'q{q}'] = {
        'B': float(B), 'A': float(A), 'fit_R2': float(fit_R2),
        'nA_values': [int(x) for x in nA_arr],
        'one_minus_R2': [float(x) for x in omr_arr],
    }
    record(sprint=94, model='sq_potts', q=q, n=0,
           quantity='bw_exp_rate_B', value=B, method='094c_exp_fit')

# Summary
print("\n\n=== EXPONENTIAL RATE B(q) ===")
print(f"{'q':>3} {'B':>8} {'B/ln(q)':>10} {'A':>12} {'fit_R²':>8}")
print("-" * 48)
B_vals = []
for q in sorted(data.keys()):
    f = results['fits'][f'q{q}']
    B_over_lnq = f['B'] / np.log(q)
    print(f"{q:>3} {f['B']:>8.3f} {B_over_lnq:>10.3f} {f['A']:>12.2e} {f['fit_R2']:>8.4f}")
    B_vals.append(f['B'])

# Does B scale with q or ln(q)?
q_arr = np.array(sorted(data.keys()), dtype=float)
B_arr = np.array(B_vals)

# Linear fit B vs q
coeffs_lin = np.polyfit(q_arr, B_arr, 1)
pred_lin = np.polyval(coeffs_lin, q_arr)
R2_lin = 1 - np.sum((B_arr - pred_lin)**2) / np.sum((B_arr - np.mean(B_arr))**2)

# Linear fit B vs ln(q)
lnq_arr = np.log(q_arr)
coeffs_log = np.polyfit(lnq_arr, B_arr, 1)
pred_log = np.polyval(coeffs_log, lnq_arr)
R2_log = 1 - np.sum((B_arr - pred_log)**2) / np.sum((B_arr - np.mean(B_arr))**2)

print(f"\nB vs q linear fit: B = {coeffs_lin[0]:.3f}*q + {coeffs_lin[1]:.3f}, R² = {R2_lin:.4f}")
print(f"B vs ln(q) fit:   B = {coeffs_log[0]:.3f}*ln(q) + {coeffs_log[1]:.3f}, R² = {R2_log:.4f}")

results['B_vs_q'] = {
    'q_values': [int(x) for x in q_arr],
    'B_values': [float(x) for x in B_arr],
    'linear_slope': float(coeffs_lin[0]),
    'linear_intercept': float(coeffs_lin[1]),
    'linear_R2': float(R2_lin),
    'log_slope': float(coeffs_log[0]),
    'log_intercept': float(coeffs_log[1]),
    'log_R2': float(R2_log),
}

# Cross-validate: predict nA where 1-R²=0.5 (BW essentially fails)
print("\n\n=== Predicted nA* where BW fails (1-R²=0.5) ===")
for q in sorted(data.keys()):
    f = results['fits'][f'q{q}']
    nA_star = (np.log(0.5) - np.log(f['A'])) / f['B']
    print(f"  q={q}: nA* = {nA_star:.1f}")
    results['fits'][f'q{q}']['nA_star_50pct'] = float(nA_star)

# nA* where 1-R²=0.01 (1% BW inaccuracy)
print("\n=== Predicted nA* where 1-R²=0.01 (1% inaccuracy) ===")
for q in sorted(data.keys()):
    f = results['fits'][f'q{q}']
    nA_star = (np.log(0.01) - np.log(f['A'])) / f['B']
    print(f"  q={q}: nA* = {nA_star:.1f}")
    results['fits'][f'q{q}']['nA_star_1pct'] = float(nA_star)

with open("results/sprint_094c3_expfit.json", "w") as f:
    json.dump(results, f, indent=2, default=str)

print("\nDone!")
