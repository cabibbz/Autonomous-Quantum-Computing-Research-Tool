#!/usr/bin/env python3
"""Sprint 099c: Deep pairwise analysis — separate Casimir from velocity.

The key signal from 099b: q=5 pairwise c/Re(c) is non-monotonic.
Question: is the non-monotonicity in vc (pure Casimir) or in v (velocity)?
And is the reversal quantitatively consistent with complex CFT omega?
"""
import numpy as np
import json, time

with open("results/sprint_099a_casimir_dense.json") as f:
    data = json.load(f)

results = {
    'experiment': '099c_pairwise_analysis',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_099c_pairwise_analysis.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

x_sigma = {2: 0.123, 5: 0.136, 7: 0.132}
rec_values = {2: 0.500, 5: 1.138, 7: 1.351}

print("Sprint 099c: Pairwise Casimir Analysis — Separating vc from v")
print("=" * 70)

for q_str in ['2', '5', '7']:
    q = int(q_str)
    d = data['data'][q_str]
    N_arr = np.array(d['sizes'], dtype=float)
    y = np.array(d['E0_per_N'])
    gapN = np.array(d['gap_N'])
    rec = rec_values[q]
    xs = x_sigma[q]

    print(f"\n{'='*70}")
    print(f"q={q}, Re(c)={rec}")
    print(f"{'='*70}")

    # Pairwise vc (pure Casimir, no velocity needed)
    print(f"\n  {'pair':>8} {'vc':>10} {'v_geo':>10} {'c_imp':>10} {'c/Rec':>10} {'Δ(c/Rec)':>10}")
    pw_vc = []
    pw_v = []
    pw_c = []
    pw_ratio = []

    for i in range(len(N_arr) - 1):
        N1, N2 = N_arr[i], N_arr[i+1]
        dE = y[i+1] - y[i]
        dx = 1/N2**2 - 1/N1**2
        vc_pair = -dE / dx * 6 / np.pi

        # Velocity from geometric mean of the two sizes
        v1 = gapN[i] / (2*np.pi*xs)
        v2 = gapN[i+1] / (2*np.pi*xs)
        v_geo = np.sqrt(v1 * v2)  # geometric mean
        c_pair = vc_pair / v_geo
        ratio = c_pair / rec

        pw_vc.append(float(vc_pair))
        pw_v.append(float(v_geo))
        pw_c.append(float(c_pair))
        pw_ratio.append(float(ratio))

        delta = ratio - pw_ratio[-2] if len(pw_ratio) > 1 else 0.0
        print(f"  ({int(N1):2d},{int(N2):2d}) {vc_pair:10.6f} {v_geo:10.6f} "
              f"{c_pair:10.6f} {ratio:10.6f} {delta:+10.6f}")

    # Check monotonicity of vc itself
    vc_arr = np.array(pw_vc)
    vc_diffs = np.diff(vc_arr)
    vc_mono = 'MONO-' if all(d < 0 for d in vc_diffs) else ('MONO+' if all(d > 0 for d in vc_diffs) else 'NON-MONO')

    # Check monotonicity of v
    v_arr = np.array(pw_v)
    v_diffs = np.diff(v_arr)
    v_mono = 'MONO-' if all(d < 0 for d in v_diffs) else ('MONO+' if all(d > 0 for d in v_diffs) else 'NON-MONO')

    # Check monotonicity of c/Re(c)
    r_arr = np.array(pw_ratio)
    r_diffs = np.diff(r_arr)
    r_mono = 'MONO-' if all(d < 0 for d in r_diffs) else ('MONO+' if all(d > 0 for d in r_diffs) else 'NON-MONO')

    print(f"\n  vc trend: {vc_mono}")
    print(f"  v trend:  {v_mono}")
    print(f"  c/Rec:    {r_mono}")

    # ---- KEY TEST: vc_pair convergence without velocity ----
    # The product vc should converge to v_∞ * Re(c).
    # Check vc directly — no velocity correction needed
    print(f"\n  Pairwise vc (product, no velocity separation):")
    for i, vc in enumerate(pw_vc):
        delta = pw_vc[i] - pw_vc[i-1] if i > 0 else 0
        N_mid = np.sqrt(N_arr[i] * N_arr[i+1])
        print(f"    N_mid={N_mid:.1f}: vc={vc:.8f} Δvc={delta:+.4e}")

    vc_diffs_signs = ''.join(['+' if d > 0 else '-' for d in vc_diffs])
    print(f"  vc diff signs: {vc_diffs_signs}")

    # ---- Richardson extrapolation of vc (higher order) ----
    # If vc(N) ~ vc_∞ + a/N_mid^2 + b/N_mid^4 + ...
    # Use consecutive vc pairs to eliminate leading correction
    if len(pw_vc) >= 3:
        print(f"\n  Richardson-extrapolated vc (eliminating leading correction):")
        N_mids = [np.sqrt(N_arr[i]*N_arr[i+1]) for i in range(len(N_arr)-1)]
        for i in range(len(pw_vc) - 1):
            N1, N2 = N_mids[i], N_mids[i+1]
            h1, h2 = 1/N1**2, 1/N2**2
            vc_rich = (pw_vc[i+1]*h1 - pw_vc[i]*h2) / (h1 - h2)
            print(f"    ({N1:.1f},{N2:.1f}): vc_rich = {vc_rich:.8f}")

    # ---- Detrended vc: remove power-law and look for wiggles ----
    # Fit vc(N_mid) = a + b/N_mid^2
    N_mids_arr = np.array([np.sqrt(N_arr[i]*N_arr[i+1]) for i in range(len(N_arr)-1)])
    vc_arr_fit = np.array(pw_vc)
    A_fit = np.vstack([1/N_mids_arr**2, np.ones(len(N_mids_arr))]).T
    (slope, intercept), _, _, _ = np.linalg.lstsq(A_fit, vc_arr_fit, rcond=None)
    vc_pred = slope / N_mids_arr**2 + intercept
    vc_resid = vc_arr_fit - vc_pred

    print(f"\n  Detrended vc (after removing linear 1/N² trend):")
    print(f"  vc_∞ = {intercept:.8f}, slope = {slope:.8f}")
    for i, N_mid in enumerate(N_mids_arr):
        print(f"    N_mid={N_mid:.1f}: vc_resid = {vc_resid[i]:+.4e}")

    vc_resid_signs = ''.join(['+' if r > 0 else '-' for r in vc_resid])
    print(f"  Detrended sign pattern: {vc_resid_signs}")

    vc_sign_changes = sum(1 for i in range(len(vc_resid)-1) if vc_resid[i]*vc_resid[i+1] < 0)
    print(f"  Sign changes: {vc_sign_changes}/{len(vc_resid)-1}")

    # For q>4: compare vc non-monotonicity to omega prediction
    if q > 4:
        alpha = np.arccosh(np.sqrt(q)/2)
        omega = 2 * alpha
        half_period = np.pi / omega  # in ln(N) units
        print(f"\n  Complex CFT: omega={omega:.4f}, half-period={half_period:.4f} in ln(N)")
        print(f"  Our ln(N) range: {np.log(N_arr[0]):.3f} to {np.log(N_arr[-1]):.3f} "
              f"(span = {np.log(N_arr[-1]/N_arr[0]):.3f})")
        print(f"  Fraction of period covered: {np.log(N_arr[-1]/N_arr[0]) / (2*np.pi/omega):.3f}")

    results['data'][q_str] = {
        'pairwise_vc': pw_vc,
        'pairwise_v': pw_v,
        'pairwise_c_over_rec': pw_ratio,
        'vc_monotonic': vc_mono,
        'v_monotonic': v_mono,
        'c_ratio_monotonic': r_mono,
        'vc_detrended': [float(r) for r in vc_resid],
        'vc_detrended_sign_changes': vc_sign_changes,
    }

save()

# Grand summary
print(f"\n{'='*70}")
print("SUMMARY: Where is the non-monotonicity?")
print(f"{'='*70}")
print(f"{'q':>3} {'vc':>10} {'v':>10} {'c/Rec':>10} {'vc_detrend_sign_chg':>20}")
for q_str in ['2', '5', '7']:
    d = results['data'][q_str]
    print(f"{q_str:>3} {d['vc_monotonic']:>10} {d['v_monotonic']:>10} "
          f"{d['c_ratio_monotonic']:>10} {d['vc_detrended_sign_changes']:>20}")

print(f"\nKey question: Does detrended vc show MORE oscillatory behavior for q>4?")
for q_str in ['2', '5', '7']:
    d = results['data'][q_str]
    resids = d['vc_detrended']
    if len(resids) > 2:
        # RMS of detrended vc
        rms = np.sqrt(np.mean(np.array(resids)**2))
        print(f"  q={q_str}: vc_detrend RMS = {rms:.4e}, "
              f"sign changes = {d['vc_detrended_sign_changes']}/{len(resids)-1}")

save()
print(f"\nSaved to results/sprint_099c_pairwise_analysis.json")
