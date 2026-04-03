#!/usr/bin/env python3
"""Sprint 075c: Final cylinder convergence analysis for q=2.

Combine all g_c(Ly) data with entropy c_eff(Ly) for unified picture.
Fit updated exponential model. Analyze convergence rate.
Compare g_c convergence with c_eff convergence.
"""
import numpy as np
import json
import time

results = {
    'experiment': '075c_convergence_analysis',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

# ---- g_c convergence data for q=2 ----
print("=" * 60, flush=True)
print("q=2 Cylinder g_c Convergence (Full Dataset)", flush=True)
print("=" * 60, flush=True)

gc_data = {
    1: 0.250,  # 1D (exact)
    2: 0.451,  # Ly=2 (3 crossing pairs)
    3: 0.655,  # Ly=3 (4 crossing pairs)
    4: 0.688,  # Ly=4 (2 crossing pairs)
    5: 0.701,  # Ly=5 (1 crossing pair, 075a)
}
gc_2d = 0.771

Ly_arr = np.array(sorted(gc_data.keys()))
gc_arr = np.array([gc_data[ly] for ly in Ly_arr])
progress = (gc_arr - gc_data[1]) / (gc_2d - gc_data[1]) * 100

print(f"{'Ly':>4} {'g_c':>8} {'Progress':>10} {'Increment':>10} {'Ratio':>8}")
print("-" * 45)
for i, ly in enumerate(Ly_arr):
    inc = gc_arr[i] - gc_arr[i-1] if i > 0 else 0
    ratio = inc / (gc_arr[i-1] - gc_arr[i-2]) if i > 1 else float('nan')
    print(f"{ly:>4d} {gc_arr[i]:>8.3f} {progress[i]:>9.1f}% {inc:>10.4f} {ratio:>8.3f}" if i > 1 else
          f"{ly:>4d} {gc_arr[i]:>8.3f} {progress[i]:>9.1f}% {inc:>10.4f}        —")

# Increment analysis
increments = np.diff(gc_arr)
print(f"\nIncrements: {[f'{d:.4f}' for d in increments]}")
print(f"Increment ratios: {[f'{increments[i+1]/increments[i]:.3f}' for i in range(len(increments)-1)]}")

results['gc_data'] = {str(k): v for k, v in gc_data.items()}
results['gc_2d'] = gc_2d
results['increments'] = [float(d) for d in increments]

# ---- Exponential fit: g_c(Ly) = g_c(2D) - A*exp(-Ly/B) ----
print(f"\n{'='*60}", flush=True)
print("Exponential Fit: g_c(Ly) = 0.771 - A*exp(-Ly/B)", flush=True)
print("=" * 60, flush=True)

from scipy.optimize import curve_fit

def exp_model(Ly, A, B):
    return gc_2d - A * np.exp(-np.array(Ly) / B)

# Fit using Ly=2-5 (skip Ly=1 which is exact 1D)
Ly_fit = Ly_arr[1:]  # Ly=2,3,4,5
gc_fit = gc_arr[1:]

try:
    popt, pcov = curve_fit(exp_model, Ly_fit, gc_fit, p0=[0.5, 2.0])
    A_fit, B_fit = popt
    print(f"  A = {A_fit:.4f}, B = {B_fit:.4f}")
    print(f"  g_c(Ly) = {gc_2d:.3f} - {A_fit:.4f}*exp(-Ly/{B_fit:.4f})")

    # Predictions
    for ly_pred in [6, 7, 8, 10]:
        gc_pred = exp_model(ly_pred, A_fit, B_fit)
        pct = (gc_pred - 0.250) / (gc_2d - 0.250) * 100
        print(f"  Ly={ly_pred}: g_c = {gc_pred:.4f} ({pct:.1f}% to 2D)")

    # Residuals
    gc_model = exp_model(Ly_fit, A_fit, B_fit)
    residuals = gc_fit - gc_model
    print(f"\n  Residuals: {[f'{r:.4f}' for r in residuals]}")
    print(f"  RMS residual: {np.sqrt(np.mean(residuals**2)):.5f}")

    results['exp_fit'] = {
        'A': float(A_fit), 'B': float(B_fit),
        'residuals': [float(r) for r in residuals],
        'rms': float(np.sqrt(np.mean(residuals**2))),
    }
except Exception as e:
    print(f"  Fit failed: {e}")

# ---- Power-law fit: g_c(Ly) = g_c(2D) - C/Ly^alpha ----
print(f"\n{'='*60}", flush=True)
print("Power-Law Fit: g_c(Ly) = 0.771 - C/Ly^alpha", flush=True)
print("=" * 60, flush=True)

def power_model(Ly, C, alpha):
    return gc_2d - C / np.array(Ly) ** alpha

try:
    popt2, pcov2 = curve_fit(power_model, Ly_fit, gc_fit, p0=[0.5, 1.0])
    C_fit, alpha_fit = popt2
    print(f"  C = {C_fit:.4f}, alpha = {alpha_fit:.4f}")
    print(f"  g_c(Ly) = {gc_2d:.3f} - {C_fit:.4f}/Ly^{alpha_fit:.4f}")

    gc_model2 = power_model(Ly_fit, C_fit, alpha_fit)
    residuals2 = gc_fit - gc_model2
    print(f"\n  Residuals: {[f'{r:.4f}' for r in residuals2]}")
    print(f"  RMS residual: {np.sqrt(np.mean(residuals2**2)):.5f}")

    for ly_pred in [6, 7, 8, 10]:
        gc_pred = power_model(ly_pred, C_fit, alpha_fit)
        pct = (gc_pred - 0.250) / (gc_2d - 0.250) * 100
        print(f"  Ly={ly_pred}: g_c = {gc_pred:.4f} ({pct:.1f}% to 2D)")

    results['power_fit'] = {
        'C': float(C_fit), 'alpha': float(alpha_fit),
        'residuals': [float(r) for r in residuals2],
        'rms': float(np.sqrt(np.mean(residuals2**2))),
    }
except Exception as e:
    print(f"  Fit failed: {e}")

# ---- c_eff convergence ----
print(f"\n{'='*60}", flush=True)
print("c_eff Convergence on Cylinders", flush=True)
print("=" * 60, flush=True)

c_eff_data = {
    1: 0.610,  # 1D (overshoots c=0.5 at small Lx)
    2: 0.752,
    3: 0.817,
}
c_exact_1d = 0.500

print(f"{'Ly':>4} {'c_eff':>8} {'Δc_eff':>8}")
print("-" * 25)
for i, ly in enumerate(sorted(c_eff_data.keys())):
    dc = c_eff_data[ly] - c_eff_data.get(ly-1, 0) if ly > 1 else 0
    print(f"{ly:>4d} {c_eff_data[ly]:>8.3f} {dc:>8.3f}" if ly > 1 else
          f"{ly:>4d} {c_eff_data[ly]:>8.3f}       —")

# c_eff grows: 0.61 → 0.75 → 0.82
# Increments: 0.14, 0.065 — decelerating
# If pure 2D Ising c is different... 2D Ising has c=0.5 too!
# So the growth of c_eff is a finite-size effect, NOT an approach to higher c.
print(f"\nNote: 2D Ising also has c=0.5. c_eff growth is a finite-size effect.")
print(f"The Lx values used (4-12) overshoot c. Overshooting increases with Ly")
print(f"because more bonds are cut (area-law contribution masks the log).")

# ---- Convergence comparison: g_c vs c_eff ----
print(f"\n{'='*60}", flush=True)
print("Cross-Analysis: g_c Convergence vs Entropy", flush=True)
print("=" * 60, flush=True)

# At Lx=4: S values across Ly
# From 075b data
S_at_Lx4 = {
    1: 0.2849,
    2: 0.3146,
    3: 0.3218,
    4: 0.3564,
}

print(f"{'Ly':>4} {'g_c':>8} {'S(Lx=4)':>10} {'Cuts':>6}")
print("-" * 35)
for ly in sorted(S_at_Lx4.keys()):
    # Number of bonds cut at half-way (Ly*1 vertical bonds + Ly horizontal wraps)
    # Actually at x=Lx/2 cut: Ly vertical bonds cross the cut
    n_cuts = ly  # Ly vertical bonds cross
    print(f"{ly:>4d} {gc_data[ly]:>8.3f} {S_at_Lx4[ly]:>10.4f} {n_cuts:>6d}")

# S per cut
print(f"\nS per bond cut at Lx=4:")
for ly in sorted(S_at_Lx4.keys()):
    s_per_cut = S_at_Lx4[ly] / ly
    print(f"  Ly={ly}: S/Ly = {s_per_cut:.4f}")

results['S_at_Lx4'] = {str(k): v for k, v in S_at_Lx4.items()}
results['c_eff_data'] = {str(k): v for k, v in c_eff_data.items()}

# ---- Key finding: convergence deceleration ----
print(f"\n{'='*60}", flush=True)
print("KEY FINDINGS", flush=True)
print("=" * 60, flush=True)

inc = np.diff(gc_arr)
print(f"1. g_c increments DECELERATING:")
print(f"   Ly=1→2: +{inc[0]:.3f}")
print(f"   Ly=2→3: +{inc[1]:.3f}")
print(f"   Ly=3→4: +{inc[2]:.3f}")
print(f"   Ly=4→5: +{inc[3]:.3f}")
print(f"   Ratios: {inc[1]/inc[0]:.2f}, {inc[2]/inc[1]:.2f}, {inc[3]/inc[2]:.2f}")

print(f"\n2. Exponential prediction OVERSHOOTS:")
print(f"   Predicted g_c(Ly=5) = 0.745, Measured = 0.701")
print(f"   Convergence is SLOWER than exponential near 2D")

print(f"\n3. c_eff grows with Ly (0.61→0.75→0.82) but this is FSS overshoot,")
print(f"   not approach to higher c. 2D Ising has c=0.5 same as 1D.")

print(f"\n4. S per bond cut DECREASING: {S_at_Lx4[1]/1:.4f} → {S_at_Lx4[2]/2:.4f} → {S_at_Lx4[3]/3:.4f}")
print(f"   Bonds become less entangled as Ly grows — approaching area law.")

# Save
with open("results/sprint_075c_convergence.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to results/sprint_075c_convergence.json", flush=True)

from db_utils import record
record(sprint=75, model='hybrid_cyl', q=2, n=5,
       quantity='gc_gap', value=0.701,
       method='convergence_analysis_Ly5',
       notes='Ly=5 confirmed. Exp fit: A={:.3f}, B={:.3f}'.format(
           results.get('exp_fit', {}).get('A', 0),
           results.get('exp_fit', {}).get('B', 0)))
