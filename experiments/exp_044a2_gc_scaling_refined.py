#!/usr/bin/env python3
"""Sprint 044a2: Refined g_c scaling — separate self-dual and non-self-dual regimes.

Key insight: q=2,3 have g_c=1.0 exactly (self-duality). For q>=4, self-duality breaks
and g_c drops. Fit the q>=4 regime separately.

Also test: g_c(q) = 1/(q-1)^gamma  (which gives g_c=1 at q=2 automatically)
and piecewise: g_c=1 for q<=3, then decay for q>=4.
"""
import numpy as np
from scipy.optimize import curve_fit, minimize_scalar
import json

# Full data
q_all = np.array([2, 3, 4, 5, 10], dtype=float)
gc_all = np.array([1.0, 1.0, 0.893, 0.41, 0.246])

# q>=4 data only (non-self-dual regime)
q_nsd = np.array([4, 5, 10], dtype=float)
gc_nsd = np.array([0.893, 0.41, 0.246])

results = {}

# === Fits for q>=4 only ===
print("=== q>=4 regime (non-self-dual) ===")

# 1. Power law: g_c = A * q^(-alpha)
def pl(q, A, alpha): return A * q**(-alpha)
popt, _ = curve_fit(pl, q_nsd, gc_nsd, p0=[5, 1])
pred = pl(q_nsd, *popt)
resid = np.sqrt(np.mean((gc_nsd - pred)**2))
results['q4+_power'] = {'A': float(popt[0]), 'alpha': float(popt[1]), 'rmse': resid}
print(f"  Power: g_c = {popt[0]:.3f} * q^(-{popt[1]:.3f}), RMSE={resid:.4f}")
for q, g, p in zip(q_nsd, gc_nsd, pred):
    print(f"    q={int(q)}: data={g:.3f}, fit={p:.3f}, err={abs(g-p):.3f}")

# 2. g_c = A / (q-q0)  (pole form)
def pole(q, A, q0): return A / (q - q0)
popt2, _ = curve_fit(pole, q_nsd, gc_nsd, p0=[3, 1])
pred2 = pole(q_nsd, *popt2)
resid2 = np.sqrt(np.mean((gc_nsd - pred2)**2))
results['q4+_pole'] = {'A': float(popt2[0]), 'q0': float(popt2[1]), 'rmse': resid2}
print(f"  Pole: g_c = {popt2[0]:.3f} / (q - {popt2[1]:.3f}), RMSE={resid2:.4f}")
for q, g, p in zip(q_nsd, gc_nsd, pred2):
    print(f"    q={int(q)}: data={g:.3f}, fit={p:.3f}, err={abs(g-p):.3f}")

# 3. Exponential: g_c = A * exp(-B*(q-3))
def exp_decay(q, A, B): return A * np.exp(-B * (q - 3))
popt3, _ = curve_fit(exp_decay, q_nsd, gc_nsd, p0=[1, 0.3])
pred3 = exp_decay(q_nsd, *popt3)
resid3 = np.sqrt(np.mean((gc_nsd - pred3)**2))
results['q4+_exp'] = {'A': float(popt3[0]), 'B': float(popt3[1]), 'rmse': resid3}
print(f"  Exp: g_c = {popt3[0]:.3f} * exp(-{popt3[1]:.3f}*(q-3)), RMSE={resid3:.4f}")
for q, g, p in zip(q_nsd, gc_nsd, pred3):
    print(f"    q={int(q)}: data={g:.3f}, fit={p:.3f}, err={abs(g-p):.3f}")

# === Unified fits (all q) ===
print("\n=== Unified fits (all q) ===")

# 4. g_c = (q-1)^(-gamma) — automatically gives g_c(2)=1
def potts_unified(q, gamma):
    return (q - 1)**(-gamma)

# Scan gamma to find best
def chi2_gamma(gamma):
    pred = potts_unified(q_all, gamma)
    return np.sum((gc_all - pred)**2)

from scipy.optimize import minimize_scalar
res = minimize_scalar(chi2_gamma, bounds=(0.1, 3.0), method='bounded')
gamma_opt = res.x
pred_u = potts_unified(q_all, gamma_opt)
resid_u = np.sqrt(np.mean((gc_all - pred_u)**2))
results['unified_qm1_power'] = {'gamma': gamma_opt, 'rmse': resid_u}
print(f"  (q-1)^(-gamma): gamma={gamma_opt:.4f}, RMSE={resid_u:.4f}")
for q, g, p in zip(q_all, gc_all, pred_u):
    print(f"    q={int(q)}: data={g:.3f}, fit={p:.3f}, err={abs(g-p):.3f}")

# 5. g_c = min(1, A*(q-1)^(-gamma))  — clipped at 1
def potts_clipped(q, A, gamma):
    return np.minimum(1.0, A * (q - 1)**(-gamma))

# Grid search since clipping makes it non-smooth
best = (1e10, 1, 1)
for A in np.linspace(0.5, 3.0, 100):
    for gamma in np.linspace(0.3, 3.0, 100):
        pred = potts_clipped(q_all, A, gamma)
        chi2 = np.sum((gc_all - pred)**2)
        if chi2 < best[0]:
            best = (chi2, A, gamma)
chi2_best, A_best, gamma_best = best
pred_c = potts_clipped(q_all, A_best, gamma_best)
resid_c = np.sqrt(np.mean((gc_all - pred_c)**2))
results['clipped_power'] = {'A': A_best, 'gamma': gamma_best, 'rmse': resid_c}
print(f"  Clipped: g_c = min(1, {A_best:.3f}*(q-1)^(-{gamma_best:.3f})), RMSE={resid_c:.4f}")
for q, g, p in zip(q_all, gc_all, pred_c):
    print(f"    q={int(q)}: data={g:.3f}, fit={p:.3f}, err={abs(g-p):.3f}")

# 6. Known exact: quantum Potts self-dual at g_c = √(q) / (1 + √(q-1))? Let's check
# Actually, for 1D quantum Potts: H = -J Σ P_ij - g Σ (X + X†)
# The self-dual point is at g/J = 1 for this normalization when q≤4c (second-order)
# For our convention, the exact g_c should satisfy a self-duality relation.
# The key formula: for q-state Potts in 1D quantum, g_c = 1/√(q) * correction
# Let's test g_c = c / sqrt(q)
def sqrt_q(q, c):
    return c / np.sqrt(q)

popt_sq, _ = curve_fit(sqrt_q, q_nsd, gc_nsd, p0=[2.0])
pred_sq = sqrt_q(q_nsd, *popt_sq)
resid_sq = np.sqrt(np.mean((gc_nsd - pred_sq)**2))
results['q4+_inv_sqrt'] = {'c': float(popt_sq[0]), 'rmse': resid_sq}
print(f"\n  1/√q: g_c = {popt_sq[0]:.3f} / √q, RMSE={resid_sq:.4f}")
for q, g, p in zip(q_nsd, gc_nsd, pred_sq):
    print(f"    q={int(q)}: data={g:.3f}, fit={p:.3f}, err={abs(g-p):.3f}")

# === Summary and predictions ===
print("\n=== RANKING ===")
ranked = sorted(results.items(), key=lambda x: x[1]['rmse'])
for i, (name, r) in enumerate(ranked):
    print(f"  {i+1}. {name}: RMSE={r['rmse']:.4f}")

print("\n=== PREDICTIONS for q=6,7,8,15,20 ===")
best_name = ranked[0][0]
q_pred = [6, 7, 8, 15, 20]
predictions = {}
for q in q_pred:
    preds = {}
    # q4+ fits
    if 'q4+_power' in results:
        r = results['q4+_power']
        preds['power'] = r['A'] * q**(-r['alpha'])
    if 'q4+_pole' in results:
        r = results['q4+_pole']
        preds['pole'] = r['A'] / (q - r['q0'])
    if 'q4+_exp' in results:
        r = results['q4+_exp']
        preds['exp'] = r['A'] * np.exp(-r['B'] * (q - 3))
    if 'clipped_power' in results:
        r = results['clipped_power']
        preds['clipped'] = min(1.0, r['A'] * (q-1)**(-r['gamma']))
    if 'q4+_inv_sqrt' in results:
        r = results['q4+_inv_sqrt']
        preds['inv_sqrt'] = r['c'] / np.sqrt(q)

    vals = list(preds.values())
    predictions[q] = preds
    print(f"  q={q}: " + ", ".join(f"{k}={v:.3f}" for k, v in sorted(preds.items())))
    print(f"         mean={np.mean(vals):.3f} ± {np.std(vals):.3f}")

# Best prediction for q=7 (our verification target)
print(f"\n=== q=7 PREDICTION SUMMARY ===")
p7 = predictions[7]
vals7 = list(p7.values())
print(f"  Mean: {np.mean(vals7):.3f}")
print(f"  Range: [{min(vals7):.3f}, {max(vals7):.3f}]")
print(f"  Std: {np.std(vals7):.3f}")
print(f"  Best model ({ranked[0][0]}): {p7.get('clipped', p7.get('power', 'N/A'))}")

# Save
output = {
    'sprint': '044a2',
    'description': 'Refined g_c scaling fits (q>=4 separate + unified)',
    'data': {'q': q_all.tolist(), 'gc': gc_all.tolist()},
    'fits': results,
    'ranking': [(n, r['rmse']) for n, r in ranked],
    'predictions': {str(q): preds for q, preds in predictions.items()},
}
with open('results/sprint_044a_gc_scaling.json', 'w') as f:
    json.dump(output, f, indent=2)

print("\nSaved to results/sprint_044a_gc_scaling.json")
