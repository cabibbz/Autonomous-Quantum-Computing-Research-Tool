#!/usr/bin/env python3
"""Sprint 099a: Dense Casimir energy scan at ALL sizes for q=2,5,7.

Generate E₀(N) at every integer N to look for oscillatory corrections
predicted by complex CFT (Im(c) ≠ 0 for q>4).

q=2: N=4-14 (control — real CFT, no oscillations expected)
q=5: N=4-10 (walking — Im(c)≈0.021, possible oscillations)
q=7: N=4-8 (broken walking — Im(c)≈0.029, possible oscillations)
"""
import numpy as np
import json, time
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh

results = {
    'experiment': '099a_casimir_dense',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'description': 'Dense Casimir scan at all integer N for oscillation detection',
    'data': {},
}

def save():
    with open("results/sprint_099a_casimir_dense.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_periodic(n, q, g):
    """Build S_q Potts Hamiltonian on periodic chain — vectorized."""
    dim = q**n
    all_idx = np.arange(dim, dtype=np.int64)
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q
        tmp //= q
    powers = q ** np.arange(n, dtype=np.int64)

    # Diagonal: coupling
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n):
        nxt = (site + 1) % n
        diag_vals -= (digits[:, site] == digits[:, nxt]).astype(np.float64)

    rows_list = [all_idx]
    cols_list = [all_idx]
    vals_list = [diag_vals]

    # Off-diagonal: field
    for site in range(n):
        pw = powers[site]
        old_digit = digits[:, site]
        for k in range(1, q):
            new_digit = (old_digit + k) % q
            delta = (new_digit.astype(np.int64) - old_digit.astype(np.int64)) * pw
            new_idx = all_idx + delta
            rows_list.append(new_idx)
            cols_list.append(all_idx)
            vals_list.append(np.full(dim, -g, dtype=np.float64))

    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    vals = np.concatenate(vals_list)
    return csr_matrix(coo_matrix((vals, (rows, cols)), shape=(dim, dim)))

# Complex CFT values
rec_values = {2: 0.500, 5: 1.138, 7: 1.351}
imc_values = {2: 0.000, 5: 0.021, 7: 0.029}  # Im(c) from complex CFT
x_sigma = {2: 0.123, 5: 0.136, 7: 0.132}

# Size lists — ALL integers
size_lists = {
    2: list(range(4, 15)),   # N=4-14, max dim=16384
    5: list(range(4, 11)),   # N=4-10, max dim=9.8M
    7: list(range(4, 9)),    # N=4-8, max dim=5.7M
}

print("Sprint 099a: Dense Casimir Energy Scan for Oscillation Detection")
print("=" * 70, flush=True)

for q in [2, 5, 7]:
    gc = 1.0 / q
    sizes = size_lists[q]
    rec = rec_values[q]
    imc = imc_values[q]
    print(f"\n{'='*70}")
    print(f"q={q}, g_c={gc:.4f}, Re(c)={rec}, Im(c)={imc}")
    print(f"Sizes: {sizes} (max dim={q**max(sizes):,})")
    print(f"{'='*70}", flush=True)

    q_data = {'q': q, 'gc': gc, 'rec': rec, 'imc': imc,
              'sizes': [], 'E0_per_N': [], 'E0_raw': [],
              'gaps': [], 'gap_N': []}

    for n in sizes:
        dim = q**n
        print(f"  n={n}, dim={dim:,} ... ", end='', flush=True)
        t0 = time.time()
        H = build_sq_potts_periodic(n, q, gc)
        t_build = time.time() - t0
        t0 = time.time()
        evals, _ = eigsh(H, k=4, which='SA')
        t_eig = time.time() - t0
        evals = np.sort(np.real(evals))
        E0 = float(evals[0])
        gap = float(evals[1] - evals[0])
        E0_per_N = E0 / n
        gapN = gap * n

        q_data['sizes'].append(n)
        q_data['E0_per_N'].append(E0_per_N)
        q_data['E0_raw'].append(E0)
        q_data['gaps'].append(gap)
        q_data['gap_N'].append(gapN)

        print(f"E0/N={E0_per_N:.12f}, gap*N={gapN:.6f}, "
              f"build={t_build:.1f}s, eig={t_eig:.1f}s", flush=True)

    # === Casimir fit: E0/N = eps_inf + A/N^2 + B/N^4 ===
    N_arr = np.array(q_data['sizes'], dtype=float)
    y = np.array(q_data['E0_per_N'])
    x2 = 1.0 / N_arr**2
    x4 = 1.0 / N_arr**4

    # 3-param fit
    A3 = np.vstack([x2, x4, np.ones_like(x2)]).T
    coeffs, _, _, _ = np.linalg.lstsq(A3, y, rcond=None)
    slope2, slope4, eps_inf = coeffs
    vc_cas = -slope2 * 6 / np.pi
    y_pred3 = slope2 * x2 + slope4 * x4 + eps_inf
    resid3 = y - y_pred3
    ss_res3 = np.sum(resid3**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    R2_3 = 1 - ss_res3 / ss_tot

    # Velocity from largest size
    v_last = q_data['gap_N'][-1] / (2 * np.pi * x_sigma[q])
    c_implied = vc_cas / v_last

    print(f"\n  3-param fit: E0/N = {eps_inf:.12f} + ({slope2:.8f})/N² + ({slope4:.8f})/N⁴")
    print(f"  R² = {R2_3:.12f}")
    print(f"  vc = {vc_cas:.8f}, v = {v_last:.6f}")
    print(f"  c_implied = {c_implied:.6f}, c/Re(c) = {c_implied/rec:.6f}")

    # === Residuals from 3-param fit ===
    print(f"\n  Residuals from 3-param fit (looking for oscillations):")
    print(f"  {'N':>4} {'ln(N)':>8} {'resid':>14} {'resid*N^3':>14} {'resid*N^5':>14}")
    for i, n_val in enumerate(q_data['sizes']):
        r = resid3[i]
        print(f"  {n_val:4d} {np.log(n_val):8.4f} {r:+14.4e} {r*n_val**3:+14.4e} {r*n_val**5:+14.4e}")

    # === Test oscillatory model: resid ~ A*cos(omega*ln(N) + phi) / N^p ===
    # For complex CFT, omega = 2*alpha where alpha = arccosh(sqrt(q)/2)
    if q > 4:
        alpha = np.arccosh(np.sqrt(q) / 2)
        omega_pred = 2 * alpha
        print(f"\n  Complex CFT prediction: alpha = {alpha:.6f}, omega = 2*alpha = {omega_pred:.6f}")
        print(f"  Period in ln(N) = 2*pi/omega = {2*np.pi/omega_pred:.4f}")
        print(f"  Full oscillation needs N ratio = exp(2pi/omega) = {np.exp(2*np.pi/omega_pred):.1f}")

    # Phase of residuals: sign changes?
    sign_changes = 0
    for i in range(len(resid3) - 1):
        if resid3[i] * resid3[i+1] < 0:
            sign_changes += 1
    print(f"\n  Residual sign changes: {sign_changes} out of {len(resid3)-1} consecutive pairs")
    print(f"  Residual range: [{min(resid3):.4e}, {max(resid3):.4e}]")

    # Pairwise c_implied
    print(f"\n  Pairwise c_implied:")
    pairwise = []
    for i in range(len(N_arr) - 1):
        N1, N2 = N_arr[i], N_arr[i+1]
        dE = y[i+1] - y[i]
        dx = 1/N2**2 - 1/N1**2
        vc_pair = -dE / dx * 6 / np.pi
        v_pair = q_data['gap_N'][i+1] / (2*np.pi*x_sigma[q])
        c_pair = vc_pair / v_pair
        pairwise.append({
            'pair': f'({int(N1)},{int(N2)})',
            'vc': float(vc_pair), 'c_implied': float(c_pair),
            'c_over_rec': float(c_pair / rec),
        })
        print(f"    ({int(N1)},{int(N2)}): c/Re(c) = {c_pair/rec:.6f}")

    q_data['fit_3param'] = {
        'eps_inf': float(eps_inf), 'slope2': float(slope2), 'slope4': float(slope4),
        'vc': float(vc_cas), 'v_last': float(v_last),
        'c_implied': float(c_implied), 'c_over_rec': float(c_implied/rec),
        'R2': float(R2_3),
    }
    q_data['residuals_3param'] = [float(r) for r in resid3]
    q_data['pairwise'] = pairwise
    q_data['sign_changes'] = sign_changes
    results['data'][str(q)] = q_data
    save()

# === Grand comparison: residual structure ===
print(f"\n{'='*70}")
print("RESIDUAL STRUCTURE COMPARISON")
print(f"{'='*70}")
for q in [2, 5, 7]:
    d = results['data'][str(q)]
    resids = d['residuals_3param']
    sizes = d['sizes']
    print(f"\nq={q} (Im(c)={d['imc']}): {d['sign_changes']} sign changes in {len(resids)-1} pairs")
    print(f"  Max |resid|: {max(abs(r) for r in resids):.4e}")
    print(f"  Residual pattern: {' '.join(['+' if r > 0 else '-' for r in resids])}")

    # Autocorrelation of residuals (simple)
    r = np.array(resids)
    if len(r) > 3:
        r_centered = r - np.mean(r)
        var = np.sum(r_centered**2)
        if var > 0:
            ac1 = np.sum(r_centered[:-1] * r_centered[1:]) / var
            ac2 = np.sum(r_centered[:-2] * r_centered[2:]) / var if len(r) > 4 else 0
            print(f"  Autocorrelation: lag-1={ac1:.3f}, lag-2={ac2:.3f}")
            print(f"  (oscillatory if lag-1 < 0 and lag-2 > 0)")

save()
print(f"\nSaved to results/sprint_099a_casimir_dense.json")
