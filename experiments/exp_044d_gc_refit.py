#!/usr/bin/env python3
"""Sprint 044d: Refit g_c(q) scaling law with new q=7 data point.

Updated data:
  q=2: g_c=1.0, q=3: g_c=1.0, q=4: g_c=0.893, q=5: g_c=0.41,
  q=7: g_c=0.259 (NEW), q=10: g_c=0.246

Focus on q>=4 regime where self-duality is broken.
"""
import numpy as np
from scipy.optimize import curve_fit, minimize_scalar
import json

# Updated data
q_all = np.array([2, 3, 4, 5, 7, 10], dtype=float)
gc_all = np.array([1.0, 1.0, 0.893, 0.41, 0.259, 0.246])

q_nsd = np.array([4, 5, 7, 10], dtype=float)
gc_nsd = np.array([0.893, 0.41, 0.259, 0.246])

results = {}

print("=== q>=4 fits (4 data points) ===")

# 1. Power law: g_c = A * q^(-alpha)
def pl(q, A, alpha): return A * q**(-alpha)
popt, _ = curve_fit(pl, q_nsd, gc_nsd, p0=[5, 1])
pred = pl(q_nsd, *popt)
rmse = np.sqrt(np.mean((gc_nsd - pred)**2))
results['power'] = {'A': float(popt[0]), 'alpha': float(popt[1]), 'rmse': rmse}
print(f"  Power: g_c = {popt[0]:.3f} * q^(-{popt[1]:.3f}), RMSE={rmse:.4f}")
for q, g, p in zip(q_nsd, gc_nsd, pred):
    print(f"    q={int(q)}: data={g:.3f}, fit={p:.3f}")

# 2. Pole: g_c = A / (q - q0)
def pole(q, A, q0): return A / (q - q0)
popt2, _ = curve_fit(pole, q_nsd, gc_nsd, p0=[3, 1])
pred2 = pole(q_nsd, *popt2)
rmse2 = np.sqrt(np.mean((gc_nsd - pred2)**2))
results['pole'] = {'A': float(popt2[0]), 'q0': float(popt2[1]), 'rmse': rmse2}
print(f"  Pole: g_c = {popt2[0]:.3f} / (q - {popt2[1]:.3f}), RMSE={rmse2:.4f}")
for q, g, p in zip(q_nsd, gc_nsd, pred2):
    print(f"    q={int(q)}: data={g:.3f}, fit={p:.3f}")

# 3. Exponential: g_c = A * exp(-B*(q-3))
def exp_d(q, A, B): return A * np.exp(-B * (q - 3))
popt3, _ = curve_fit(exp_d, q_nsd, gc_nsd, p0=[1, 0.3])
pred3 = exp_d(q_nsd, *popt3)
rmse3 = np.sqrt(np.mean((gc_nsd - pred3)**2))
results['exp'] = {'A': float(popt3[0]), 'B': float(popt3[1]), 'rmse': rmse3}
print(f"  Exp: g_c = {popt3[0]:.3f} * exp(-{popt3[1]:.3f}*(q-3)), RMSE={rmse3:.4f}")
for q, g, p in zip(q_nsd, gc_nsd, pred3):
    print(f"    q={int(q)}: data={g:.3f}, fit={p:.3f}")

# 4. Power of (q-3): g_c = A * (q-3)^(-beta) — pole at self-dual boundary
def pole3(q, A, beta): return A * (q - 3)**(-beta)
try:
    popt4, _ = curve_fit(pole3, q_nsd, gc_nsd, p0=[0.5, 0.5])
    pred4 = pole3(q_nsd, *popt4)
    rmse4 = np.sqrt(np.mean((gc_nsd - pred4)**2))
    results['pole3_power'] = {'A': float(popt4[0]), 'beta': float(popt4[1]), 'rmse': rmse4}
    print(f"  (q-3)^(-β): g_c = {popt4[0]:.3f} * (q-3)^(-{popt4[1]:.3f}), RMSE={rmse4:.4f}")
    for q, g, p in zip(q_nsd, gc_nsd, pred4):
        print(f"    q={int(q)}: data={g:.3f}, fit={p:.3f}")
except Exception as e:
    print(f"  (q-3)^(-β) failed: {e}")

# 5. Clipped (q-1)^(-gamma) — unified
def clipped(q, A, gamma):
    return np.minimum(1.0, A * (q - 1)**(-gamma))

best_c = (1e10, 1, 1)
for A in np.linspace(0.3, 5.0, 200):
    for gamma in np.linspace(0.3, 4.0, 200):
        pred = clipped(q_all, A, gamma)
        chi2 = np.sum((gc_all - pred)**2)
        if chi2 < best_c[0]:
            best_c = (chi2, A, gamma)
chi2_bc, A_bc, g_bc = best_c
pred_c = clipped(q_all, A_bc, g_bc)
rmse_c = np.sqrt(np.mean((gc_all - pred_c)**2))
results['clipped'] = {'A': A_bc, 'gamma': g_bc, 'rmse': rmse_c}
print(f"  Clipped: g_c = min(1, {A_bc:.3f}*(q-1)^(-{g_bc:.3f})), RMSE={rmse_c:.4f}")
for q, g, p in zip(q_all, gc_all, pred_c):
    print(f"    q={int(q)}: data={g:.3f}, fit={p:.3f}")

# === Key diagnostic: q=7 prediction accuracy from 044a models (without q=7) ===
print("\n=== q=7 PREDICTION vs MEASUREMENT ===")
print(f"  Measured: g_c(7) = 0.259")
# From 044a (fitted WITHOUT q=7):
pred_044a = {
    'power_q4+': 14.69 * 7**(-2.067),
    'pole_q4+': 0.996 / (7 - 2.868),
    'exp_q4+': 1.149 * np.exp(-0.357 * (7-3)),
    'clipped_044a': min(1.0, 3.0 * (7-1)**(-1.227)),
}
for name, p in pred_044a.items():
    err = abs(p - 0.259) / 0.259 * 100
    print(f"  {name}: predicted {p:.3f}, error {err:.1f}%")

# === Ranking ===
print("\n=== FINAL RANKING (by RMSE) ===")
ranked = sorted(results.items(), key=lambda x: x[1]['rmse'])
for i, (name, r) in enumerate(ranked):
    print(f"  {i+1}. {name}: RMSE={r['rmse']:.4f}")

# === What does the data tell us about g_c → 0 behavior? ===
print("\n=== ASYMPTOTIC BEHAVIOR ===")
# q=7 and q=10 are close: 0.259 vs 0.246
# Δg/Δq = (0.259-0.246)/(10-7) = 0.004/q
print(f"  g_c(7) - g_c(10) = {0.259-0.246:.3f} over Δq=3")
print(f"  Average slope: {(0.259-0.246)/3:.4f} per q")
print(f"  If linear: g_c → 0 at q ≈ {10 + 0.246/0.004:.0f}")

# Best model extrapolation
best = ranked[0]
print(f"\n  Best model: {best[0]}")
for q in [15, 20, 50, 100]:
    r = best[1]
    if best[0] == 'pole3_power':
        val = r['A'] * (q-3)**(-r['beta'])
    elif best[0] == 'pole':
        val = r['A'] / (q - r['q0'])
    elif best[0] == 'power':
        val = r['A'] * q**(-r['alpha'])
    elif best[0] == 'exp':
        val = r['A'] * np.exp(-r['B']*(q-3))
    elif best[0] == 'clipped':
        val = min(1.0, r['A'] * (q-1)**(-r['gamma']))
    else:
        val = float('nan')
    print(f"    g_c({q}) ≈ {val:.4f}")

# Save
output = {
    'sprint': '044d',
    'description': 'g_c refit with q=7 data point',
    'data': {'q': q_all.tolist(), 'gc': gc_all.tolist()},
    'fits': results,
    'ranking': [(n, r['rmse']) for n, r in ranked],
    'q7_prediction_test': pred_044a,
    'q7_measured': 0.259,
}
with open('results/sprint_044d_gc_refit.json', 'w') as f:
    json.dump(output, f, indent=2)

print("\nSaved to results/sprint_044d_gc_refit.json")
