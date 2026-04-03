#!/usr/bin/env python3
"""Sprint 099b: Test oscillatory vs monotonic Casimir corrections.

Compare three models for E₀/N:
  M1: ε + A/N² + B/N⁴  (3-param, monotonic)
  M2: ε + A/N² + B/N⁴ + D/N⁶  (4-param, monotonic)
  M3: ε + A/N² + C·cos(ω·ln(N)+φ)/N⁴  (4-param, oscillatory with ω FIXED from complex CFT)

If Im(c) ≠ 0 causes oscillations, M3 should beat M2 for q=5,7 but NOT for q=2.
"""
import numpy as np
import json, time
from scipy.optimize import curve_fit

# Load data
with open("results/sprint_099a_casimir_dense.json") as f:
    data = json.load(f)

results = {
    'experiment': '099b_oscillatory_fit',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_099b_oscillatory_fit.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# Complex CFT parameters
# alpha = arccosh(sqrt(q)/2), omega = 2*alpha
def get_omega(q):
    if q <= 4:
        return 0.0  # real CFT
    return 2 * np.arccosh(np.sqrt(q) / 2)

def get_imc(q):
    """Im(c) from complex CFT."""
    if q <= 4:
        return 0.0
    alpha = np.arccosh(np.sqrt(q) / 2)
    p_real = np.pi**2 / (np.pi**2 + alpha**2)  # Re part of p
    p_imag = np.pi * alpha / (np.pi**2 + alpha**2)
    # c = 1 - 6/(p(p-1)), p complex
    p_complex = complex(p_real * np.pi / alpha, np.pi**2 / (np.pi**2 + alpha**2) * np.pi / alpha)
    # Actually use the exact formula
    alpha_val = alpha
    p = complex(0, np.pi / alpha_val)
    c = 1 - 6 / (p * (p - 1))
    return c.imag

print("Sprint 099b: Oscillatory vs Monotonic Casimir Corrections")
print("=" * 70, flush=True)

for q_str in ['2', '5', '7']:
    q = int(q_str)
    d = data['data'][q_str]
    N_arr = np.array(d['sizes'], dtype=float)
    y = np.array(d['E0_per_N'])
    n_pts = len(N_arr)
    omega = get_omega(q)

    print(f"\n{'='*70}")
    print(f"q={q}, {n_pts} data points, omega={omega:.4f}")
    print(f"{'='*70}")

    # ---- Model 1: ε + A/N² + B/N⁴ ----
    x2 = 1.0 / N_arr**2
    x4 = 1.0 / N_arr**4
    x6 = 1.0 / N_arr**6
    A1 = np.vstack([x2, x4, np.ones(n_pts)]).T
    c1, _, _, _ = np.linalg.lstsq(A1, y, rcond=None)
    y_pred1 = A1 @ c1
    resid1 = y - y_pred1
    rss1 = np.sum(resid1**2)
    bic1 = n_pts * np.log(rss1 / n_pts) + 3 * np.log(n_pts)

    # ---- Model 2: ε + A/N² + B/N⁴ + D/N⁶ ----
    A2 = np.vstack([x2, x4, x6, np.ones(n_pts)]).T
    c2, _, _, _ = np.linalg.lstsq(A2, y, rcond=None)
    y_pred2 = A2 @ c2
    resid2 = y - y_pred2
    rss2 = np.sum(resid2**2)
    bic2 = n_pts * np.log(rss2 / n_pts + 1e-300) + 4 * np.log(n_pts)

    # ---- Model 3: ε + A/N² + C·cos(ω·ln(N)+φ)/N⁴ ----
    # With omega FIXED from complex CFT prediction
    if omega > 0:
        def model3(N, eps, A, C, phi):
            return eps + A / N**2 + C * np.cos(omega * np.log(N) + phi) / N**4

        try:
            # Initial guess from M1
            p0 = [c1[2], c1[0], 0.0, 0.0]
            popt3, pcov3 = curve_fit(model3, N_arr, y, p0=p0, maxfev=10000)
            y_pred3 = model3(N_arr, *popt3)
            resid3 = y - y_pred3
            rss3 = np.sum(resid3**2)
            bic3 = n_pts * np.log(rss3 / n_pts + 1e-300) + 4 * np.log(n_pts)
            C_amp = popt3[2]
            phi_fit = popt3[3]
        except Exception as e:
            print(f"  M3 fit failed: {e}")
            rss3, bic3, C_amp, phi_fit = rss1, bic1, 0, 0
            resid3 = resid1
    else:
        # For q=2 (real CFT), use a "fake" oscillatory model with omega as free param
        def model3_free(N, eps, A, C, phi, w):
            return eps + A / N**2 + C * np.cos(w * np.log(N) + phi) / N**4

        try:
            p0 = [c1[2], c1[0], 0.0, 0.0, 1.0]
            popt3f, pcov3f = curve_fit(model3_free, N_arr, y, p0=p0, maxfev=10000)
            y_pred3f = model3_free(N_arr, *popt3f)
            resid3f = y - y_pred3f
            rss3 = np.sum(resid3f**2)
            bic3 = n_pts * np.log(rss3 / n_pts + 1e-300) + 5 * np.log(n_pts)
            C_amp = popt3f[2]
            phi_fit = popt3f[3]
            omega_fit = popt3f[4]
            resid3 = resid3f
            print(f"  q=2 free omega fit: omega={omega_fit:.3f}, C={C_amp:.4e}")
        except:
            rss3, bic3, C_amp, phi_fit = rss1, bic1, 0, 0
            resid3 = resid1

    # ---- Model 4: ε + A/N² + C·cos(ω·ln(N)+φ)/N⁴ + B/N⁴ (5-param, mixed) ----
    if omega > 0 and n_pts >= 6:
        def model4(N, eps, A, B, C, phi):
            return eps + A / N**2 + B / N**4 + C * np.cos(omega * np.log(N) + phi) / N**4

        try:
            p0 = [c1[2], c1[0], c1[1], 0.0, 0.0]
            popt4, pcov4 = curve_fit(model4, N_arr, y, p0=p0, maxfev=10000)
            y_pred4 = model4(N_arr, *popt4)
            resid4 = y - y_pred4
            rss4 = np.sum(resid4**2)
            bic4 = n_pts * np.log(rss4 / n_pts + 1e-300) + 5 * np.log(n_pts)
            C_amp4 = popt4[3]
            phi_fit4 = popt4[4]
        except:
            rss4, bic4 = rss2, bic2
            C_amp4 = 0
    else:
        rss4, bic4, C_amp4 = rss2, bic2, 0

    # ---- Summary ----
    print(f"\n  Model comparison (RSS = residual sum of squares):")
    print(f"  M1 (3p, ε+A/N²+B/N⁴):             RSS = {rss1:.4e}, BIC = {bic1:.2f}")
    print(f"  M2 (4p, ε+A/N²+B/N⁴+D/N⁶):        RSS = {rss2:.4e}, BIC = {bic2:.2f}")
    if omega > 0:
        print(f"  M3 (4p, ε+A/N²+C·cos/N⁴):         RSS = {rss3:.4e}, BIC = {bic3:.2f}, C = {C_amp:.4e}")
    else:
        print(f"  M3 (5p, free ω):                    RSS = {rss3:.4e}, BIC = {bic3:.2f}, C = {C_amp:.4e}")

    if omega > 0 and n_pts >= 6:
        print(f"  M4 (5p, ε+A/N²+B/N⁴+C·cos/N⁴):   RSS = {rss4:.4e}, BIC = {bic4:.2f}, C = {C_amp4:.4e}")

    # F-test: M2 vs M1
    if rss1 > 0 and rss2 > 0:
        dof1 = n_pts - 3
        dof2 = n_pts - 4
        if dof2 > 0 and rss2 > 0:
            F_21 = ((rss1 - rss2) / 1) / (rss2 / dof2)
            print(f"\n  F-test M2 vs M1: F = {F_21:.2f} (dof {dof1},{dof2})")

    # Compare oscillatory amplitude
    if omega > 0:
        print(f"\n  Oscillatory amplitude C = {C_amp:.4e}")
        print(f"  At N=6: oscillatory term = {abs(C_amp)/6**4:.4e}")
        print(f"  Ratio oscillatory/systematic (|C|/|B_M1|): {abs(C_amp)/abs(c1[1]):.4f}" if c1[1] != 0 else "")

    # Residual sign pattern
    sign_pattern = ''.join(['+' if r > 0 else '-' for r in resid1])
    print(f"\n  M1 residual pattern: {sign_pattern}")
    sign_pattern2 = ''.join(['+' if r > 0 else '-' for r in resid2])
    print(f"  M2 residual pattern: {sign_pattern2}")
    sign_pattern3 = ''.join(['+' if r > 0 else '-' for r in resid3])
    print(f"  M3 residual pattern: {sign_pattern3}")

    # Pairwise c_implied and its DERIVATIVE (looking for oscillation in dc/dN)
    print(f"\n  Pairwise c_implied/Re(c) and curvature:")
    rec = d['rec']
    x_sig = {2: 0.123, 5: 0.136, 7: 0.132}[q]
    pw_ratios = []
    for i in range(len(N_arr) - 1):
        N1, N2 = N_arr[i], N_arr[i+1]
        dE = y[i+1] - y[i]
        dx = 1/N2**2 - 1/N1**2
        vc_pair = -dE / dx * 6 / np.pi
        v_pair = d['gap_N'][i+1] / (2*np.pi*x_sig)
        c_pair = vc_pair / v_pair
        ratio = c_pair / rec
        pw_ratios.append(ratio)
        print(f"    ({int(N1)},{int(N2)}): c/Re(c) = {ratio:.6f}")

    # Check monotonicity of pairwise ratios
    pw_arr = np.array(pw_ratios)
    pw_diffs = np.diff(pw_arr)
    pw_signs = ''.join(['+' if d > 0 else '-' for d in pw_diffs])
    print(f"  Pairwise differences sign: {pw_signs}")
    monotonic = all(d <= 0 for d in pw_diffs) or all(d >= 0 for d in pw_diffs)
    print(f"  Monotonic: {'YES' if monotonic else 'NO — possible oscillation signal!'}")

    if not monotonic and omega > 0:
        # Count reversals
        reversals = sum(1 for i in range(len(pw_diffs)-1) if pw_diffs[i]*pw_diffs[i+1] < 0)
        print(f"  Reversals in pairwise convergence: {reversals}")

    # Store results
    q_results = {
        'n_points': n_pts,
        'omega': float(omega),
        'M1_rss': float(rss1), 'M1_bic': float(bic1),
        'M2_rss': float(rss2), 'M2_bic': float(bic2),
        'M3_rss': float(rss3), 'M3_bic': float(bic3),
        'M3_C_amplitude': float(C_amp),
        'pairwise_ratios': [float(r) for r in pw_ratios],
        'pairwise_monotonic': bool(monotonic),
        'M1_resid_signs': sign_pattern,
    }
    if omega > 0 and n_pts >= 6:
        q_results['M4_rss'] = float(rss4)
        q_results['M4_bic'] = float(bic4)
        q_results['M4_C_amplitude'] = float(C_amp4)
    results['data'][q_str] = q_results

save()

# === Grand summary ===
print(f"\n{'='*70}")
print("GRAND SUMMARY: Oscillation Evidence")
print(f"{'='*70}")
print(f"{'q':>3} {'Im(c)':>7} {'omega':>7} {'M2_RSS':>12} {'M3_RSS':>12} {'M3/M2':>7} {'C_osc':>12} {'monotonic':>10}")
for q_str in ['2', '5', '7']:
    d = results['data'][q_str]
    ratio = d['M3_rss'] / d['M2_rss'] if d['M2_rss'] > 0 else 999
    imc = {2: 0.0, 5: 0.021, 7: 0.029}[int(q_str)]
    print(f"{q_str:>3} {imc:7.3f} {d['omega']:7.3f} {d['M2_rss']:12.4e} {d['M3_rss']:12.4e} "
          f"{ratio:7.3f} {d['M3_C_amplitude']:12.4e} {'YES' if d['pairwise_monotonic'] else 'NO':>10}")

print(f"\nIf oscillations from Im(c):")
print(f"  - M3/M2 < 1 for q>4 (oscillatory fits better than monotonic)")
print(f"  - M3/M2 ≈ 1 for q=2 (no advantage)")
print(f"  - pairwise non-monotonic for q>4")
print(f"  - C_osc increases with Im(c)")

save()
print(f"\nSaved to results/sprint_099b_oscillatory_fit.json")
