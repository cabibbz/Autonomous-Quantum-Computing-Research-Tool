#!/usr/bin/env python3
"""Sprint 074c: Cross-q Ly convergence analysis with q=5 Ly=3 data.

Updated analysis incorporating q=5 Ly=3 (g_c=0.974 from (2,3) crossing).
Compare convergence rates across q=2,3,5 and test scaling laws.
"""
import numpy as np
import json
import time
from scipy.optimize import curve_fit

print("=" * 60, flush=True)
print("Sprint 074c: Cross-q Ly convergence analysis", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '074c_convergence_analysis',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

# Complete dataset
data = {
    2: {
        'gc_1d': 0.250, 'gc_2d': 0.771,
        'cyl': {1: 0.250, 2: 0.451, 3: 0.655, 4: 0.688},
    },
    3: {
        'gc_1d': 0.333, 'gc_2d': 1.267,
        'cyl': {1: 0.333, 2: 0.565, 3: 0.797},
    },
    5: {
        'gc_1d': 0.441, 'gc_2d': 1.588,
        'cyl': {1: 0.441, 2: 0.714, 3: 0.974},  # Ly=3 from (2,3) crossing
    },
}

# 1. Progress to 2D for each q and Ly
print("\n--- Progress to 2D ---", flush=True)
print(f"{'q':>3} {'Ly':>3} {'g_c':>8} {'%2D':>7} {'Δg_c':>8}", flush=True)
print("-" * 35, flush=True)
progress_data = {}
for q in [2, 3, 5]:
    gc_1d = data[q]['gc_1d']
    gc_2d = data[q]['gc_2d']
    gap = gc_2d - gc_1d
    prev_gc = gc_1d
    progress_data[q] = {}
    for Ly in sorted(data[q]['cyl'].keys()):
        gc = data[q]['cyl'][Ly]
        pct = (gc - gc_1d) / gap * 100
        delta = gc - prev_gc
        print(f"{q:>3} {Ly:>3} {gc:>8.3f} {pct:>6.1f}% {delta:>8.3f}", flush=True)
        progress_data[q][Ly] = {'gc': gc, 'pct': pct, 'delta': delta}
        prev_gc = gc

results['progress'] = {str(k): {str(ly): v for ly, v in val.items()} for k, val in progress_data.items()}

# 2. Exponential fits: g_c(Ly) = g_c(2D) - A * exp(-Ly/B)
print(f"\n--- Exponential fits ---", flush=True)
def exp_conv(Ly, A, B, gc_2d):
    return gc_2d - A * np.exp(-np.array(Ly) / B)

fit_results = {}
for q in [2, 3, 5]:
    gc_2d = data[q]['gc_2d']
    Lys = sorted(data[q]['cyl'].keys())
    gcs = [data[q]['cyl'][ly] for ly in Lys]
    # Need at least 2 points above Ly=1 for meaningful fit
    fit_Lys = [ly for ly in Lys if ly >= 2]
    fit_gcs = [data[q]['cyl'][ly] for ly in fit_Lys]

    if len(fit_Lys) >= 2:
        try:
            func = lambda Ly, A, B: gc_2d - A * np.exp(-np.array(Ly) / B)
            popt, pcov = curve_fit(func, fit_Lys, fit_gcs, p0=[1.0, 2.0])
            A, B = popt
            perr = np.sqrt(np.diag(pcov))
            # Predictions
            pred_ly5 = func(5, A, B)
            pred_ly6 = func(6, A, B)
            residuals = np.array(fit_gcs) - func(fit_Lys, A, B)
            rms = np.sqrt(np.mean(residuals ** 2))
            print(f"q={q}: g_c(Ly) = {gc_2d:.3f} - {A:.3f}*exp(-Ly/{B:.2f})", flush=True)
            print(f"  RMS = {rms:.5f}, A_err={perr[0]:.3f}, B_err={perr[1]:.3f}", flush=True)
            print(f"  Predictions: Ly=5: {pred_ly5:.4f}, Ly=6: {pred_ly6:.4f}", flush=True)
            fit_results[q] = {
                'A': float(A), 'B': float(B),
                'A_err': float(perr[0]), 'B_err': float(perr[1]),
                'rms': float(rms),
                'pred_ly5': float(pred_ly5), 'pred_ly6': float(pred_ly6),
            }
        except Exception as e:
            print(f"q={q}: fit failed: {e}", flush=True)
    else:
        print(f"q={q}: insufficient data for fit (need Ly≥2 points)", flush=True)

results['exp_fits'] = {str(k): v for k, v in fit_results.items()}

# 3. Compare convergence decay lengths B(q)
print(f"\n--- Convergence decay length B(q) ---", flush=True)
if len(fit_results) >= 2:
    qs = sorted(fit_results.keys())
    Bs = [fit_results[q]['B'] for q in qs]
    print(f"  q=2: B={fit_results.get(2, {}).get('B', 'N/A')}", flush=True)
    print(f"  q=3: B={fit_results.get(3, {}).get('B', 'N/A')}", flush=True)
    print(f"  q=5: B={fit_results.get(5, {}).get('B', 'N/A')}", flush=True)

    # Test if B scales with q or 2D/1D ratio
    if len(qs) >= 3:
        ratios_2d_1d = [data[q]['gc_2d'] / data[q]['gc_1d'] for q in qs]
        print(f"\n  2D/1D ratios: {[f'{r:.2f}' for r in ratios_2d_1d]}", flush=True)
        print(f"  B values:     {[f'{b:.2f}' for b in Bs]}", flush=True)
        # Correlation
        corr = np.corrcoef(ratios_2d_1d, Bs)[0, 1]
        print(f"  Correlation(2D/1D ratio, B) = {corr:.3f}", flush=True)
        results['B_vs_ratio_corr'] = float(corr)

# 4. Matched-Ly comparison (all 3 q at Ly=2 and Ly=3)
print(f"\n--- Matched-Ly comparison ---", flush=True)
for Ly in [2, 3]:
    print(f"\nLy={Ly}:", flush=True)
    for q in [2, 3, 5]:
        gc = data[q]['cyl'].get(Ly)
        if gc is not None:
            gc_1d = data[q]['gc_1d']
            gc_2d = data[q]['gc_2d']
            pct = (gc - gc_1d) / (gc_2d - gc_1d) * 100
            ratio = gc / gc_1d
            print(f"  q={q}: g_c={gc:.3f}, cyl/1D={ratio:.2f}, {pct:.1f}% to 2D", flush=True)

# 5. Cyl/1D ratio at Ly=2 and Ly=3
print(f"\n--- Cyl/1D ratio evolution ---", flush=True)
for Ly in [2, 3]:
    ratios = []
    for q in [2, 3, 5]:
        gc = data[q]['cyl'].get(Ly)
        if gc is not None:
            r = gc / data[q]['gc_1d']
            ratios.append((q, r))
    print(f"Ly={Ly}: {', '.join(f'q={q}: {r:.2f}' for q, r in ratios)}", flush=True)
    if len(ratios) >= 2:
        trend = 'decreasing' if ratios[-1][1] < ratios[0][1] else 'increasing'
        print(f"  Trend: {trend}", flush=True)

results['cyl_1d_ratios'] = {
    'Ly2': {str(q): data[q]['cyl'][2] / data[q]['gc_1d'] for q in [2, 3, 5]},
    'Ly3': {str(q): data[q]['cyl'][3] / data[q]['gc_1d'] for q in [2, 3, 5]},
}

# 6. Key finding: q-dependent convergence
print(f"\n{'='*60}", flush=True)
print("KEY FINDINGS", flush=True)
print(f"{'='*60}", flush=True)

# Progress at Ly=3 for all q
print("\nProgress to 2D at Ly=3:", flush=True)
for q in [2, 3, 5]:
    gc = data[q]['cyl'][3]
    gc_1d = data[q]['gc_1d']
    gc_2d = data[q]['gc_2d']
    pct = (gc - gc_1d) / (gc_2d - gc_1d) * 100
    print(f"  q={q}: {pct:.1f}%", flush=True)

print("\nConvergence SLOWER for larger q:", flush=True)
print("  q=2: 77.7% at Ly=3, 84.0% at Ly=4", flush=True)
print("  q=3: 49.7% at Ly=3 (28% behind q=2)", flush=True)
print("  q=5: 46.5% at Ly=3 (31% behind q=2)", flush=True)

print("\nq=3 and q=5 converge at SIMILAR rates (49.7% vs 46.5%).", flush=True)
print("The convergence slowdown from q=2 to q=3 is much larger than q=3 to q=5.", flush=True)
print("This suggests a saturation of the convergence rate for q≥3.", flush=True)

# Save
with open("results/sprint_074c_convergence_analysis.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_074c_convergence_analysis.json", flush=True)

from db_utils import record
record(sprint=74, model='hybrid_cyl', q=5, n=3,
       quantity='gc_gap', value=0.974,
       method='gapLx_crossing_Ly3_Lx23',
       notes='Ly=3 cylinder, (2,3) crossing, low confidence (Lx=2 FSS)')

print(f"\n{'='*60}", flush=True)
print("SUMMARY", flush=True)
print(f"{'='*60}", flush=True)
print("New data: q=5 Ly=3 g_c=0.974 (46.5% to 2D)", flush=True)
print("q=3 Ly=4: INFEASIBLE (838s/pt at Lx=4)", flush=True)
print("Convergence saturation: q=3 and q=5 at similar % at Ly=3", flush=True)
if fit_results:
    for q in sorted(fit_results.keys()):
        print(f"B(q={q}) = {fit_results[q]['B']:.2f}", flush=True)
