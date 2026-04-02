#!/usr/bin/env python3
"""Sprint 073c: Cross-q Ly convergence analysis.

Compare how q=2 and q=3 critical couplings converge to 2D as Ly increases.
Fit exponential and power-law models. Predict 2D g_c(q=3) from convergence rate.

Also compute: effective coordination number, universal Ly scaling.
"""
import numpy as np
import json
import time

results = {
    'experiment': '073c_ly_convergence',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

print("=" * 60, flush=True)
print("Sprint 073c: Cross-q Ly convergence analysis", flush=True)
print("=" * 60, flush=True)

# Collected data
# q=2: g_c(1D)=0.250, Ly=2: 0.451, Ly=3: 0.655, Ly=4: 0.688, 2D: 0.771
# q=3: g_c(1D)=0.333, Ly=2: 0.565, Ly=3: 0.797, 2D: 1.267

q2_data = {
    'gc_1d': 0.250,
    'gc_ly': {1: 0.250, 2: 0.451, 3: 0.655, 4: 0.688},
    'gc_2d': 0.771,
}
q3_data = {
    'gc_1d': 0.333,
    'gc_ly': {1: 0.333, 2: 0.565, 3: 0.797},
    'gc_2d': 1.267,
}

# ---- Analysis 1: Normalized convergence ----
print("\n--- Normalized convergence: f(Ly) = [g_c(Ly) - g_c(1D)] / [g_c(2D) - g_c(1D)] ---", flush=True)

for label, data in [('q=2', q2_data), ('q=3', q3_data)]:
    gc_1d = data['gc_1d']
    gc_2d = data['gc_2d']
    span = gc_2d - gc_1d
    print(f"\n{label}: span = {gc_2d} - {gc_1d} = {span:.3f}", flush=True)
    for Ly in sorted(data['gc_ly'].keys()):
        gc = data['gc_ly'][Ly]
        f = (gc - gc_1d) / span
        print(f"  Ly={Ly}: g_c={gc:.4f}, f={f:.4f} ({f*100:.1f}%)", flush=True)

# ---- Analysis 2: Exponential fit for q=2 ----
print("\n--- Exponential fit g_c(Ly) = g_c(2D) - A*exp(-Ly/B) ---", flush=True)
from scipy.optimize import curve_fit

def exp_conv(Ly, A, B, gc_2d):
    return gc_2d - A * np.exp(-np.array(Ly) / B)

# q=2: fit with known gc_2d
Lys_q2 = [1, 2, 3, 4]
gcs_q2 = [q2_data['gc_ly'][ly] for ly in Lys_q2]

def exp_conv_q2(Ly, A, B):
    return 0.771 - A * np.exp(-np.array(Ly) / B)

popt2, pcov2 = curve_fit(exp_conv_q2, Lys_q2, gcs_q2, p0=[0.5, 1.5])
A2, B2 = popt2
print(f"\nq=2: g_c(Ly) = 0.771 - {A2:.4f}*exp(-Ly/{B2:.3f})", flush=True)
for Ly in range(1, 7):
    pred = exp_conv_q2(Ly, A2, B2)
    actual = q2_data['gc_ly'].get(Ly)
    err_str = f"  (actual: {actual:.4f}, err: {abs(pred-actual):.4f})" if actual else ""
    print(f"  Ly={Ly}: g_c_pred = {pred:.4f}{err_str}", flush=True)

# q=3: fit with 3 points, gc_2d = 1.267
Lys_q3 = [1, 2, 3]
gcs_q3 = [q3_data['gc_ly'][ly] for ly in Lys_q3]

def exp_conv_q3(Ly, A, B):
    return 1.267 - A * np.exp(-np.array(Ly) / B)

popt3, pcov3 = curve_fit(exp_conv_q3, Lys_q3, gcs_q3, p0=[0.5, 1.5])
A3, B3 = popt3
print(f"\nq=3: g_c(Ly) = 1.267 - {A3:.4f}*exp(-Ly/{B3:.3f})", flush=True)
for Ly in range(1, 7):
    pred = exp_conv_q3(Ly, A3, B3)
    actual = q3_data['gc_ly'].get(Ly)
    err_str = f"  (actual: {actual:.4f}, err: {abs(pred-actual):.4f})" if actual else ""
    print(f"  Ly={Ly}: g_c_pred = {pred:.4f}{err_str}", flush=True)

results['exp_fit_q2'] = {'A': float(A2), 'B': float(B2), 'gc_2d': 0.771}
results['exp_fit_q3'] = {'A': float(A3), 'B': float(B3), 'gc_2d': 1.267}

# ---- Analysis 3: Alternative — fit gc_2d as free parameter for q=3 ----
print("\n--- Free gc_2d fit for q=3 ---", flush=True)
try:
    popt3f, pcov3f = curve_fit(exp_conv, Lys_q3, gcs_q3, p0=[0.5, 1.5, 1.3], maxfev=5000)
    A3f, B3f, gc2d_3f = popt3f
    perr = np.sqrt(np.diag(pcov3f))
    print(f"q=3 free fit: g_c(Ly) = {gc2d_3f:.4f} - {A3f:.4f}*exp(-Ly/{B3f:.3f})", flush=True)
    print(f"  Predicted g_c(2D) = {gc2d_3f:.4f} ± {perr[2]:.4f}", flush=True)
    print(f"  Known g_c(2D) = 1.267", flush=True)
    print(f"  Error = {abs(gc2d_3f - 1.267):.4f} ({abs(gc2d_3f-1.267)/1.267*100:.1f}%)", flush=True)
    results['free_fit_q3'] = {'A': float(A3f), 'B': float(B3f),
                               'gc_2d_pred': float(gc2d_3f), 'gc_2d_err': float(perr[2])}
except Exception as e:
    print(f"  Free fit failed: {e}", flush=True)

# ---- Analysis 4: Convergence ratio ----
print("\n--- Convergence ratio r(Ly) = [g_c(2D) - g_c(Ly)] / [g_c(2D) - g_c(Ly-1)] ---", flush=True)
for label, data in [('q=2', q2_data), ('q=3', q3_data)]:
    gc_2d = data['gc_2d']
    print(f"\n{label}:", flush=True)
    prev_deficit = None
    for Ly in sorted(data['gc_ly'].keys()):
        gc = data['gc_ly'][Ly]
        deficit = gc_2d - gc
        if prev_deficit is not None and prev_deficit > 0:
            ratio = deficit / prev_deficit
            print(f"  Ly={Ly}: deficit={deficit:.4f}, ratio={ratio:.4f}", flush=True)
        else:
            print(f"  Ly={Ly}: deficit={deficit:.4f}", flush=True)
        prev_deficit = deficit

# ---- Analysis 5: Effective coordination ----
print("\n--- Effective coordination number z_eff(Ly) ---", flush=True)
print("z = 2 + 2·(1 - 1/Ly) for Ly≥2 (periodic y, open x interior)", flush=True)
for Ly in range(1, 7):
    if Ly == 1:
        z = 2
    else:
        z = 2 + 2 * (1 - 1/Ly)
    print(f"  Ly={Ly}: z_eff = {z:.3f}", flush=True)

# ---- Analysis 6: Compare convergence exponents ----
print("\n--- Convergence exponent: does q=3 decay length differ from q=2? ---", flush=True)
print(f"q=2 decay length B = {B2:.3f} (in Ly units)", flush=True)
print(f"q=3 decay length B = {B3:.3f} (in Ly units)", flush=True)
if B3 > B2:
    print(f"q=3 converges SLOWER: B3/B2 = {B3/B2:.2f}", flush=True)
else:
    print(f"q=3 converges FASTER: B3/B2 = {B3/B2:.2f}", flush=True)

# ---- Analysis 7: Rescaled plot data ----
print("\n--- Rescaled: g_c(Ly)/g_c(1D) vs Ly ---", flush=True)
for label, data in [('q=2', q2_data), ('q=3', q3_data)]:
    gc_1d = data['gc_1d']
    gc_2d = data['gc_2d']
    print(f"\n{label}: (2D/1D ratio = {gc_2d/gc_1d:.3f})", flush=True)
    for Ly in sorted(data['gc_ly'].keys()):
        gc = data['gc_ly'][Ly]
        print(f"  Ly={Ly}: g_c/g_c(1D) = {gc/gc_1d:.4f}", flush=True)
    print(f"  2D:  g_c/g_c(1D) = {gc_2d/gc_1d:.4f}", flush=True)

results['convergence_comparison'] = {
    'q2': {
        'decay_length': float(B2),
        'progress': {str(ly): float((gc - q2_data['gc_1d']) / (q2_data['gc_2d'] - q2_data['gc_1d']))
                     for ly, gc in q2_data['gc_ly'].items()},
        'ratio_2d_1d': q2_data['gc_2d'] / q2_data['gc_1d'],
    },
    'q3': {
        'decay_length': float(B3),
        'progress': {str(ly): float((gc - q3_data['gc_1d']) / (q3_data['gc_2d'] - q3_data['gc_1d']))
                     for ly, gc in q3_data['gc_ly'].items()},
        'ratio_2d_1d': q3_data['gc_2d'] / q3_data['gc_1d'],
    },
}

# ---- Summary ----
print(f"\n{'='*60}", flush=True)
print("SUMMARY", flush=True)
print(f"{'='*60}", flush=True)
print(f"\nLy convergence comparison:", flush=True)
print(f"{'Ly':>4} | {'q=2 g_c':>8} | {'q=2 %':>6} | {'q=3 g_c':>8} | {'q=3 %':>6}", flush=True)
print(f"{'-'*4}-+-{'-'*8}-+-{'-'*6}-+-{'-'*8}-+-{'-'*6}", flush=True)

all_Lys = sorted(set(list(q2_data['gc_ly'].keys()) + list(q3_data['gc_ly'].keys())))
for Ly in all_Lys:
    gc2 = q2_data['gc_ly'].get(Ly)
    gc3 = q3_data['gc_ly'].get(Ly)
    f2 = f"{(gc2 - q2_data['gc_1d']) / (q2_data['gc_2d'] - q2_data['gc_1d'])*100:.1f}" if gc2 else "—"
    f3 = f"{(gc3 - q3_data['gc_1d']) / (q3_data['gc_2d'] - q3_data['gc_1d'])*100:.1f}" if gc3 else "—"
    gc2_s = f"{gc2:.4f}" if gc2 else "—"
    gc3_s = f"{gc3:.4f}" if gc3 else "—"
    print(f"{Ly:>4} | {gc2_s:>8} | {f2:>6} | {gc3_s:>8} | {f3:>6}", flush=True)

print(f"  2D | {q2_data['gc_2d']:>8.4f} | {'100.0':>6} | {q3_data['gc_2d']:>8.4f} | {'100.0':>6}", flush=True)
print(f"\nq=2 exp decay length B = {B2:.3f}", flush=True)
print(f"q=3 exp decay length B = {B3:.3f}", flush=True)
print(f"q=3 converges {'SLOWER' if B3 > B2 else 'FASTER'} by factor {B3/B2:.2f}", flush=True)
print(f"\n2D/1D g_c ratio: q=2 ({q2_data['gc_2d']/q2_data['gc_1d']:.3f}), q=3 ({q3_data['gc_2d']/q3_data['gc_1d']:.3f})", flush=True)
print(f"q=3 has {q3_data['gc_2d']/q3_data['gc_1d']/(q2_data['gc_2d']/q2_data['gc_1d']):.2f}x larger 2D/1D ratio — more room to close", flush=True)

# Save
with open("results/sprint_073c_ly_convergence.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_073c_ly_convergence.json", flush=True)

from db_utils import record
record(sprint=73, model='hybrid_cyl', q=2, n=4,
       quantity='ly_decay_length', value=float(B2),
       method='exp_fit_Ly1-4',
       notes=f'g_c(Ly) = 0.771 - {A2:.3f}*exp(-Ly/{B2:.3f})')
record(sprint=73, model='hybrid_cyl', q=3, n=3,
       quantity='ly_decay_length', value=float(B3),
       method='exp_fit_Ly1-3',
       notes=f'g_c(Ly) = 1.267 - {A3:.3f}*exp(-Ly/{B3:.3f})')
