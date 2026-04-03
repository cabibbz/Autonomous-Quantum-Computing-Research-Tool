#!/usr/bin/env python3
"""Sprint 099a (continued): Get q=7 n=8 and combine with saved data."""
import numpy as np
import json, time
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh

def build_sq_potts_periodic(n, q, g):
    dim = q**n
    all_idx = np.arange(dim, dtype=np.int64)
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q
        tmp //= q
    powers = q ** np.arange(n, dtype=np.int64)
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n):
        nxt = (site + 1) % n
        diag_vals -= (digits[:, site] == digits[:, nxt]).astype(np.float64)
    rows_list = [all_idx]
    cols_list = [all_idx]
    vals_list = [diag_vals]
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

# q=7 n=8
q, n = 7, 8
gc = 1.0/q
print(f"Computing q={q} n={n}, dim={q**n:,} ...", flush=True)
t0 = time.time()
H = build_sq_potts_periodic(n, q, gc)
t_build = time.time() - t0
print(f"  Build: {t_build:.1f}s", flush=True)

t0 = time.time()
evals, _ = eigsh(H, k=4, which='SA')
t_eig = time.time() - t0
evals = np.sort(np.real(evals))
E0 = float(evals[0])
gap = float(evals[1] - evals[0])
print(f"  Eigsh: {t_eig:.1f}s")
print(f"  E0/N = {E0/n:.12f}, gap*N = {gap*n:.6f}")

# Load and update saved results
with open("results/sprint_099a_casimir_dense.json") as f:
    results = json.load(f)

# q=7 data from the partial run (n=4-7)
q7_partial = results['data'].get('7', None)
if q7_partial is None:
    # Build from Sprint 098a data
    with open("results/sprint_098a_casimir_gpu_extended.json") as f:
        old = json.load(f)
    q7_old = old['data']['7']
    q7 = {'q': 7, 'gc': gc, 'rec': 1.351, 'imc': 0.029,
           'sizes': q7_old['sizes'][:4],  # n=4,5,6,7
           'E0_per_N': q7_old['E0_per_N'][:4],
           'E0_raw': q7_old['E0_raw'][:4],
           'gaps': q7_old['gaps'][:4],
           'gap_N': q7_old['gap_N'][:4]}
else:
    q7 = q7_partial

# Append n=8
q7['sizes'].append(n)
q7['E0_per_N'].append(E0/n)
q7['E0_raw'].append(E0)
q7['gaps'].append(gap)
q7['gap_N'].append(gap*n)

# Fit and analyze
N_arr = np.array(q7['sizes'], dtype=float)
y = np.array(q7['E0_per_N'])
x2 = 1.0 / N_arr**2
x4 = 1.0 / N_arr**4
x_sigma_q7 = 0.132

A3 = np.vstack([x2, x4, np.ones_like(x2)]).T
coeffs, _, _, _ = np.linalg.lstsq(A3, y, rcond=None)
slope2, slope4, eps_inf = coeffs
vc_cas = -slope2 * 6 / np.pi
y_pred3 = slope2 * x2 + slope4 * x4 + eps_inf
resid3 = y - y_pred3
ss_res3 = np.sum(resid3**2)
ss_tot = np.sum((y - np.mean(y))**2)
R2_3 = 1 - ss_res3/ss_tot

v_last = q7['gap_N'][-1] / (2*np.pi*x_sigma_q7)
c_implied = vc_cas / v_last

print(f"\n3-param fit: eps_inf={eps_inf:.12f}, slope2={slope2:.8f}, slope4={slope4:.8f}")
print(f"R² = {R2_3:.12f}")
print(f"vc = {vc_cas:.8f}, v = {v_last:.6f}, c = {c_implied:.6f}, c/Re(c) = {c_implied/1.351:.6f}")

print(f"\nResiduals:")
for i, n_val in enumerate(q7['sizes']):
    r = resid3[i]
    print(f"  N={n_val}: resid={r:+.4e}, resid*N³={r*n_val**3:+.4e}")

sign_changes = sum(1 for i in range(len(resid3)-1) if resid3[i]*resid3[i+1] < 0)
print(f"\nSign changes: {sign_changes}/{len(resid3)-1}")

alpha = np.arccosh(np.sqrt(7)/2)
omega = 2*alpha
print(f"Complex CFT: alpha={alpha:.6f}, omega={omega:.6f}, period={2*np.pi/omega:.4f}")

# Pairwise c_implied
print(f"\nPairwise c/Re(c):")
pairwise = []
for i in range(len(N_arr)-1):
    N1, N2 = N_arr[i], N_arr[i+1]
    dE = y[i+1] - y[i]
    dx = 1/N2**2 - 1/N1**2
    vc_pair = -dE/dx*6/np.pi
    v_pair = q7['gap_N'][i+1]/(2*np.pi*x_sigma_q7)
    c_pair = vc_pair/v_pair
    pairwise.append({'pair': f'({int(N1)},{int(N2)})', 'c_over_rec': float(c_pair/1.351)})
    print(f"  ({int(N1)},{int(N2)}): {c_pair/1.351:.6f}")

q7['fit_3param'] = {
    'eps_inf': float(eps_inf), 'slope2': float(slope2), 'slope4': float(slope4),
    'vc': float(vc_cas), 'v_last': float(v_last),
    'c_implied': float(c_implied), 'c_over_rec': float(c_implied/1.351),
    'R2': float(R2_3),
}
q7['residuals_3param'] = [float(r) for r in resid3]
q7['pairwise'] = pairwise
q7['sign_changes'] = sign_changes

results['data']['7'] = q7
with open("results/sprint_099a_casimir_dense.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved updated results.")
