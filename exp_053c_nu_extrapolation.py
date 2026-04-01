#!/usr/bin/env python3
"""Sprint 053c: Extrapolate ν(q=3) from finite-size slope ratios.

The slope d(Δ·N)/dg at g_c scales as N^{1/ν}.
Ratios of consecutive sizes give effective ν(N_eff).
These should converge to exact 5/6 as N→∞.

Also: pseudo-critical drift analysis from crossing points.
Also: data collapse with larger dataset (n=4,6,8,10).
"""
import numpy as np, json, time
from scipy.optimize import minimize_scalar, minimize
from scipy.interpolate import interp1d

# Load data
with open('results/sprint_053b_gap_larger.json') as f:
    data = json.load(f)['data']

g_c = 1.0 / 3.0

# === Extract slopes at g_c for all sizes ===
print("=== Slopes d(Δ·N)/dg at g_c ===", flush=True)
slopes = {}
for nk in ['n4', 'n6', 'n8', 'n10']:
    pts = data[nk]
    n = int(nk[1:])
    g_arr = np.array([p['g'] for p in pts])
    y_arr = np.array([p['gap_x_n'] for p in pts])
    f = interp1d(g_arr, y_arr, kind='cubic')
    dg = 0.005
    slope = float((f(g_c + dg) - f(g_c - dg)) / (2 * dg))
    slopes[n] = slope
    print(f"  n={n}: slope = {slope:.4f}, log(slope) = {np.log(slope):.4f}", flush=True)

# === Direct power-law fit: slope = a * N^{1/ν} ===
print("\n=== Direct power-law fit: slope ~ N^{1/ν} ===", flush=True)
sizes = np.array(sorted(slopes.keys()))
s_arr = np.array([slopes[n] for n in sizes])
# log-log fit: ln(slope) = ln(a) + (1/ν)*ln(N)
log_n = np.log(sizes)
log_s = np.log(s_arr)

# Simple linear regression
A = np.vstack([log_n, np.ones_like(log_n)]).T
result = np.linalg.lstsq(A, log_s, rcond=None)
inv_nu, ln_a = result[0]
nu_direct = 1.0 / inv_nu
print(f"  1/ν = {inv_nu:.4f}, ν = {nu_direct:.4f}", flush=True)
print(f"  Exact ν = 5/6 = {5/6:.4f}", flush=True)
print(f"  Error: {abs(nu_direct - 5/6)/(5/6)*100:.1f}%", flush=True)

# With corrections: slope = a * N^{1/ν} * (1 + b/N)
print("\n=== Corrected fit: slope ~ N^{1/ν} * (1 + b/N) ===", flush=True)
def fit_corrected(params):
    a, inv_nu, b = params
    pred = a * sizes**inv_nu * (1 + b/sizes)
    return np.sum((pred - s_arr)**2 / s_arr**2)

best = minimize(fit_corrected, [1.0, 1.2, 1.0], method='Nelder-Mead')
a_fit, inv_nu_fit, b_fit = best.x
nu_corrected = 1.0 / inv_nu_fit
print(f"  1/ν = {inv_nu_fit:.4f}, ν = {nu_corrected:.4f}, b = {b_fit:.4f}", flush=True)
print(f"  Exact ν = {5/6:.4f}", flush=True)
print(f"  Error: {abs(nu_corrected - 5/6)/(5/6)*100:.1f}%", flush=True)

# === Pairwise ν and extrapolation ===
print("\n=== Pairwise ν estimates ===", flush=True)
pairs = []
for i, n1 in enumerate(sizes):
    for n2 in sizes[i+1:]:
        ratio = slopes[n2] / slopes[n1]
        nu_est = np.log(n2 / n1) / np.log(ratio)
        n_eff = np.sqrt(n1 * n2)
        pairs.append((n1, n2, n_eff, nu_est))
        print(f"  ({n1},{n2}): ratio={ratio:.4f}, ν={nu_est:.4f}, N_eff={n_eff:.1f}", flush=True)

# Extrapolate pairwise ν to N→∞: ν(N) = ν_∞ + c/N
# Use consecutive pairs only
consec = [(p[2], p[3]) for p in pairs if p[1] == p[0] + 2]  # consecutive
print("\nConsecutive pairs:", flush=True)
for n_eff, nu in consec:
    print(f"  N_eff={n_eff:.1f}: ν={nu:.4f}", flush=True)

if len(consec) >= 2:
    n_effs = np.array([c[0] for c in consec])
    nus = np.array([c[1] for c in consec])

    # Linear fit: ν = ν_∞ + c/N_eff
    A2 = np.vstack([1.0/n_effs, np.ones_like(n_effs)]).T
    res2 = np.linalg.lstsq(A2, nus, rcond=None)
    c_fit, nu_inf = res2[0]
    print(f"\n  Extrapolated: ν(∞) = {nu_inf:.4f} + {c_fit:.2f}/N", flush=True)
    print(f"  Exact: {5/6:.4f}", flush=True)

    # Also try ν = ν_∞ + c/N^2
    A3 = np.vstack([1.0/n_effs**2, np.ones_like(n_effs)]).T
    res3 = np.linalg.lstsq(A3, nus, rcond=None)
    c_fit2, nu_inf2 = res3[0]
    print(f"  Extrapolated (1/N²): ν(∞) = {nu_inf2:.4f} + {c_fit2:.2f}/N²", flush=True)

# === All-pairs ν extrapolation ===
print("\n=== All-pairs extrapolation ===", flush=True)
all_n_effs = np.array([p[2] for p in pairs])
all_nus = np.array([p[3] for p in pairs])
A4 = np.vstack([1.0/all_n_effs, np.ones_like(all_n_effs)]).T
res4 = np.linalg.lstsq(A4, all_nus, rcond=None)
print(f"  All-pairs 1/N: ν(∞) = {res4[0][1]:.4f}", flush=True)

A5 = np.vstack([1.0/all_n_effs**2, np.ones_like(all_n_effs)]).T
res5 = np.linalg.lstsq(A5, all_nus, rcond=None)
print(f"  All-pairs 1/N²: ν(∞) = {res5[0][1]:.4f}", flush=True)

# === Pseudo-critical point drift ===
print("\n=== Pseudo-critical point drift ===", flush=True)

# Find crossings from raw data
crossings = {}
for n1k, n2k in [('n4','n6'), ('n4','n8'), ('n6','n8'), ('n6','n10'), ('n8','n10')]:
    d1, d2 = data[n1k], data[n2k]
    g1 = np.array([p['g'] for p in d1])
    y1 = np.array([p['gap_x_n'] for p in d1])
    g2 = np.array([p['g'] for p in d2])
    y2 = np.array([p['gap_x_n'] for p in d2])

    g_min = max(g1.min(), g2.min())
    g_max = min(g1.max(), g2.max())
    g_common = np.linspace(g_min, g_max, 500)
    f1 = interp1d(g1, y1, kind='cubic')
    f2 = interp1d(g2, y2, kind='cubic')
    diff = f1(g_common) - f2(g_common)
    for j in range(len(diff) - 1):
        if diff[j] * diff[j+1] < 0:
            gc = g_common[j] - diff[j] * (g_common[j+1] - g_common[j]) / (diff[j+1] - diff[j])
            n1, n2 = int(n1k[1:]), int(n2k[1:])
            n_eff = np.sqrt(n1 * n2)
            crossings[(n1, n2)] = gc
            shift = g_c - gc
            print(f"  ({n1},{n2}): g_cross = {gc:.5f}, shift = {shift:.5f}, N_eff={n_eff:.1f}", flush=True)
            break

# Fit: g_c - g_cross(N) ~ N^{-1/ν - ω}
# For q=3 Potts, ω = 4/5 (correction exponent from irrelevant operator)
# So shift ~ N^{-1/ν - 4/5} = N^{-6/5 - 4/5} = N^{-2} for ν=5/6
cross_neffs = []
cross_shifts = []
for (n1, n2), gc in sorted(crossings.items()):
    n_eff = np.sqrt(n1 * n2)
    shift = g_c - gc
    if shift > 0:
        cross_neffs.append(n_eff)
        cross_shifts.append(shift)

if len(cross_neffs) >= 2:
    cross_neffs = np.array(cross_neffs)
    cross_shifts = np.array(cross_shifts)
    # log-log fit
    log_n_cross = np.log(cross_neffs)
    log_s_cross = np.log(cross_shifts)
    A6 = np.vstack([log_n_cross, np.ones_like(log_n_cross)]).T
    res6 = np.linalg.lstsq(A6, log_s_cross, rcond=None)
    exponent = res6[0][0]
    print(f"\n  Shift exponent: {exponent:.4f}", flush=True)
    print(f"  Expected for ν=5/6, ω=4/5: -(1/ν + ω) = -{6/5 + 4/5:.4f} = -2.000", flush=True)
    # Extract ν assuming ω=4/5
    nu_from_shift = 1.0 / (-exponent - 4/5)
    print(f"  ν from shift (assuming ω=4/5): {nu_from_shift:.4f}", flush=True)

# === Data collapse with n=4,6,8,10, fixed g_c ===
print("\n=== Data collapse (4 sizes, fixed g_c=1/3) ===", flush=True)

def collapse_quality(nu, g_c, size_keys, data):
    curves = {}
    for nk in size_keys:
        n = int(nk[1:])
        pts = data[nk]
        x = np.array([(p['g'] - g_c) * n**(1.0/nu) for p in pts])
        y = np.array([p['gap_x_n'] for p in pts])
        curves[n] = (x, y)

    total_err = 0.0
    n_pairs = 0
    sizes_list = sorted(curves.keys())
    for i, n1 in enumerate(sizes_list):
        for n2 in sizes_list[i+1:]:
            x1, y1 = curves[n1]
            x2, y2 = curves[n2]
            x_min = max(x1.min(), x2.min())
            x_max = min(x1.max(), x2.max())
            if x_min >= x_max:
                continue
            x_common = np.linspace(x_min + 0.01*(x_max-x_min), x_max - 0.01*(x_max-x_min), 50)
            f1 = interp1d(x1, y1, kind='cubic', bounds_error=False, fill_value=np.nan)
            f2 = interp1d(x2, y2, kind='cubic', bounds_error=False, fill_value=np.nan)
            y1i, y2i = f1(x_common), f2(x_common)
            mask = np.isfinite(y1i) & np.isfinite(y2i)
            if mask.sum() < 5:
                continue
            y_mean = (np.abs(y1i[mask]) + np.abs(y2i[mask])) / 2
            y_mean = np.maximum(y_mean, 1e-10)
            err = np.mean(((y1i[mask] - y2i[mask]) / y_mean)**2)
            total_err += err
            n_pairs += 1

    return total_err / max(n_pairs, 1)

# Full scan
size_keys = ['n4', 'n6', 'n8', 'n10']
nu_scan = np.linspace(0.4, 2.0, 161)
q_scan = [collapse_quality(nu, g_c, size_keys, data) for nu in nu_scan]
q_scan = np.array(q_scan)
best_idx = np.argmin(q_scan)
print(f"  All 4 sizes: best ν = {nu_scan[best_idx]:.4f} (quality = {q_scan[best_idx]:.6f})", flush=True)

# Excluding smallest size
size_keys2 = ['n6', 'n8', 'n10']
q_scan2 = [collapse_quality(nu, g_c, size_keys2, data) for nu in nu_scan]
q_scan2 = np.array(q_scan2)
best_idx2 = np.argmin(q_scan2)
print(f"  Excluding n=4: best ν = {nu_scan[best_idx2]:.4f} (quality = {q_scan2[best_idx2]:.6f})", flush=True)

# At exact ν=5/6
q_exact = collapse_quality(5/6, g_c, size_keys, data)
q_exact2 = collapse_quality(5/6, g_c, size_keys2, data)
print(f"  Quality at ν=5/6: all={q_exact:.6f}, excl n=4={q_exact2:.6f}", flush=True)

# Relative quality
print(f"\n  Ratio best/exact (all sizes): {q_scan[best_idx] / q_exact:.4f}", flush=True)
print(f"  Ratio best/exact (excl n=4): {q_scan2[best_idx2] / q_exact2:.4f}", flush=True)

# Save results
results_summary = {
    'sprint': '053c', 'q': 3, 'g_c_exact': g_c, 'nu_exact': 5/6,
    'direct_powerlaw_nu': float(nu_direct),
    'corrected_powerlaw_nu': float(nu_corrected),
    'pairwise_nu': {f"({p[0]},{p[1]})": float(p[3]) for p in pairs},
    'extrapolated_nu_1_over_N': float(nu_inf) if len(consec) >= 2 else None,
    'collapse_nu_all4': float(nu_scan[best_idx]),
    'collapse_nu_excl_n4': float(nu_scan[best_idx2]),
    'crossing_points': {f"({k[0]},{k[1]})": float(v) for k, v in crossings.items()},
}
with open('results/sprint_053c_nu_extrapolation.json', 'w') as f:
    json.dump(results_summary, f, indent=2)

print("\n=== SUMMARY ===", flush=True)
print(f"  Exact ν(q=3) = 5/6 = {5/6:.4f}", flush=True)
print(f"  Direct power-law: ν = {nu_direct:.4f} ({abs(nu_direct-5/6)/(5/6)*100:.1f}% error)", flush=True)
print(f"  Corrected power-law: ν = {nu_corrected:.4f} ({abs(nu_corrected-5/6)/(5/6)*100:.1f}% error)", flush=True)
if len(consec) >= 2:
    print(f"  1/N extrapolation: ν = {nu_inf:.4f} ({abs(nu_inf-5/6)/(5/6)*100:.1f}% error)", flush=True)
print(f"  Data collapse (4 sizes): ν = {nu_scan[best_idx]:.4f} ({abs(nu_scan[best_idx]-5/6)/(5/6)*100:.1f}% error)", flush=True)
print(f"  Data collapse (excl n=4): ν = {nu_scan[best_idx2]:.4f} ({abs(nu_scan[best_idx2]-5/6)/(5/6)*100:.1f}% error)", flush=True)

print("\nDone!", flush=True)
