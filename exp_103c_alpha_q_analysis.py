#!/usr/bin/env python3
"""Sprint 103c: Full alpha(q) analysis combining Sprint 102 + 103 data.

Complete picture of chi_F scaling exponent across walking boundary.
Tests: alpha(q) functional form, crossover behavior at q=4-5.
"""
import numpy as np
import json, time

results = {
    'experiment': '103c_alpha_q_analysis',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_103c_alpha_q_analysis.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# All chi_F data: Sprint 102 + 103
all_data = {
    2: {6: 0.73, 8: 0.95, 10: 1.19, 12: 1.43, 14: 1.67},
    3: {6: 4.09, 8: 6.12, 10: 8.27},
    4: {6: 13.98, 8: 22.75},
    5: {6: 36.294, 8: 66.179, 9: 84.685, 10: 105.653},
    6: {6: 82.313, 8: 162.805},
    7: {6: 166.877, 7: 250.353, 8: 357.602},
}

print("=" * 60)
print("Sprint 103c: Complete alpha(q) analysis")
print("=" * 60)

# FSS fits for each q
print("\n--- FSS: chi_F ~ N^alpha ---")
print(f"{'q':>3} {'#pts':>5} {'alpha':>8} {'nu_eff':>8} {'alpha_err':>10}")
print("-" * 45)

alpha_vals = {}
nu_vals = {}

for q in sorted(all_data.keys()):
    ns = np.array(sorted(all_data[q].keys()))
    chis = np.array([all_data[q][n] for n in ns])

    log_n = np.log(ns)
    log_chi = np.log(chis)

    # Fit with error estimation (need >2 points for cov)
    if len(ns) > 2:
        coeffs, cov = np.polyfit(log_n, log_chi, 1, cov=True)
        alpha_err = np.sqrt(cov[0, 0])
    else:
        coeffs = np.polyfit(log_n, log_chi, 1)
        alpha_err = float('nan')
    alpha = coeffs[0]
    nu = 2 / (alpha + 1)

    alpha_vals[q] = alpha
    nu_vals[q] = nu

    print(f"{q:>3} {len(ns):>5} {alpha:>8.4f} {nu:>8.4f} {alpha_err:>10.4f}")

    # Pairwise
    pairwise = []
    for i in range(len(ns)-1):
        a = (np.log(chis[i+1]) - np.log(chis[i])) / (np.log(ns[i+1]) - np.log(ns[i]))
        pairwise.append(float(a))

    results['data'][str(q)] = {
        'sizes': ns.tolist(),
        'chi_F': chis.tolist(),
        'alpha': float(alpha),
        'alpha_err': float(alpha_err),
        'nu_eff': float(nu),
        'pairwise_alpha': pairwise,
    }

# Compare with known nu
print("\n--- Comparison with known exponents ---")
print(f"{'q':>3} {'nu_exact':>10} {'alpha_pred':>12} {'alpha_meas':>12} {'deviation':>10}")
print("-" * 55)
known_nu = {2: 1.0, 3: 5/6, 4: 2/3}
for q in [2, 3, 4]:
    nu_ex = known_nu[q]
    alpha_pred = 2/nu_ex - 1
    alpha_meas = alpha_vals[q]
    dev = (alpha_meas - alpha_pred) / alpha_pred * 100
    print(f"{q:>3} {nu_ex:>10.4f} {alpha_pred:>12.4f} {alpha_meas:>12.4f} {dev:>9.1f}%")

# Alpha(q) curve: fit functional form
print("\n--- Alpha(q) functional form ---")
qs = np.array(sorted(alpha_vals.keys()))
alphas = np.array([alpha_vals[q] for q in qs])

# Test: alpha(q) = a + b*q for q>=4 (walking regime)
mask45 = qs >= 4
if sum(mask45) >= 2:
    qs_w = qs[mask45]
    als_w = alphas[mask45]
    coeffs_lin = np.polyfit(qs_w, als_w, 1)
    print(f"Linear fit (q>=4): alpha = {coeffs_lin[0]:.4f} * q + {coeffs_lin[1]:.4f}")
    for q in qs_w:
        pred = coeffs_lin[0]*q + coeffs_lin[1]
        meas = alpha_vals[q]
        print(f"  q={q}: pred={pred:.3f}, meas={meas:.3f}, diff={meas-pred:+.3f}")

# Test: alpha(q) = 2 + c*(q-5) for q>=5 (super-first-order)
mask5 = qs >= 5
if sum(mask5) >= 2:
    qs_w = qs[mask5]
    als_w = alphas[mask5]
    coeffs_5 = np.polyfit(qs_w - 5, als_w - 2, 1)
    print(f"\nSuper-first-order fit (q>=5): alpha = 2 + {coeffs_5[0]:.4f}*(q-5) + {coeffs_5[1]:.4f}")
    results['super_first_order_slope'] = float(coeffs_5[0])

# Chi_F at n=6 vs q: exponential growth
print("\n--- chi_F(n=6) vs q ---")
qs_6 = []
chis_6 = []
for q in sorted(all_data.keys()):
    if 6 in all_data[q]:
        qs_6.append(q)
        chis_6.append(all_data[q][6])
        print(f"  q={q}: chi_F(n=6) = {all_data[q][6]:.4f}")

qs_6 = np.array(qs_6)
chis_6 = np.array(chis_6)
log_chi_6 = np.log(chis_6)
coeffs_exp = np.polyfit(qs_6, log_chi_6, 1)
print(f"  Fit: chi_F(n=6) ~ exp({coeffs_exp[0]:.4f} * q + {coeffs_exp[1]:.4f})")
print(f"  Growth rate: {np.exp(coeffs_exp[0]):.4f}x per unit q")
results['chi_n6_growth_rate'] = float(np.exp(coeffs_exp[0]))

# Alpha(q) summary
print("\n--- Summary: alpha(q) across walking boundary ---")
print(f"{'q':>3} {'alpha':>8} {'nu_eff':>8} {'regime':>20}")
print("-" * 45)
for q in sorted(alpha_vals.keys()):
    a = alpha_vals[q]
    nu = nu_vals[q]
    if q <= 3:
        regime = "real CFT"
    elif q == 4:
        regime = "BKT crossover"
    elif q == 5:
        regime = "walking (sweet spot)"
    else:
        regime = "broken walking"
    print(f"{q:>3} {a:>8.3f} {nu:>8.3f} {regime:>20}")

results['summary'] = {
    'q_values': sorted(alpha_vals.keys()),
    'alpha_values': {str(q): float(alpha_vals[q]) for q in sorted(alpha_vals.keys())},
    'nu_eff_values': {str(q): float(nu_vals[q]) for q in sorted(nu_vals.keys())},
}

save()
print("\nResults saved to results/sprint_103c_alpha_q_analysis.json")
