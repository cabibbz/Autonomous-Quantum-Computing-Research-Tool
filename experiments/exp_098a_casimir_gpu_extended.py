#!/usr/bin/env python3
"""Sprint 098a: GPU-extended Casimir energy extraction for q=3,5,7,8.

Extends Sprint 083 data with larger system sizes enabled by GPU eigsh:
  q=3: add n=12 (dim=531k)
  q=5: add n=10 (dim=9.8M)
  q=7: add n=8 (dim=5.7M)
  q=8: add n=7 (dim=5.7M)

Uses vectorized Hamiltonian builder for large matrices.
"""
import numpy as np
import json, time
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh

results = {
    'experiment': '098a_casimir_gpu_extended',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_098a_casimir_gpu_extended.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_periodic_fast(n, q, g):
    """Build S_q Potts Hamiltonian on periodic chain — vectorized."""
    dim = q**n
    all_idx = np.arange(dim, dtype=np.int64)

    # Decode all states: digits[i, site] = spin at site for state i
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q
        tmp //= q

    # Powers of q for index encoding
    powers = q ** np.arange(n, dtype=np.int64)

    # --- Diagonal: coupling = -Σ δ(s_i, s_{i+1}) ---
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n):
        nxt = (site + 1) % n
        diag_vals -= (digits[:, site] == digits[:, nxt]).astype(np.float64)

    rows_list = [all_idx]
    cols_list = [all_idx]
    vals_list = [diag_vals]

    # --- Off-diagonal: field = -g Σ_{site} Σ_{k=1}^{q-1} |shifted><original| ---
    for site in range(n):
        pw = powers[site]
        old_digit = digits[:, site]  # shape (dim,)
        for k in range(1, q):
            new_digit = (old_digit + k) % q
            # New index = old index + (new_digit - old_digit) * q^site
            delta = (new_digit.astype(np.int64) - old_digit.astype(np.int64)) * pw
            new_idx = all_idx + delta
            rows_list.append(new_idx)
            cols_list.append(all_idx)
            vals_list.append(np.full(dim, -g, dtype=np.float64))

    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    vals = np.concatenate(vals_list)

    H = coo_matrix((vals, (rows, cols)), shape=(dim, dim))
    return csr_matrix(H)

# Complete size lists including new GPU-enabled sizes
size_lists = {
    2: [6, 8, 10, 12, 14],        # 5 points (from 083a)
    3: [4, 6, 8, 10, 12],          # 5 points (add n=12, dim=531k)
    5: [4, 6, 8, 10],              # 4 points (add n=10, dim=9.8M)
    7: [4, 5, 6, 7, 8],            # 5 points (add n=8, dim=5.7M)
    8: [4, 5, 6, 7],               # 4 points (add n=7, dim=5.7M)
}

# Complex CFT Re(c) values
rec_values = {2: 0.500, 3: 0.800, 5: 1.138, 7: 1.351, 8: 1.438}

# c_eff values (from entropy, Sprint 079-080)
c_eff_values = {2: 0.500, 3: 0.893, 5: 1.152, 7: 1.059, 8: 1.062}

# x_sigma from Sprint 082
x_sigma_values = {2: 0.123, 3: 0.132, 5: 0.136, 7: 0.132, 8: 0.130}

print("Sprint 098a: GPU-Extended Casimir Energy Extraction")
print("=" * 70, flush=True)

for q in [2, 3, 5, 7, 8]:
    gc = 1.0 / q
    sizes = size_lists[q]
    rec = rec_values[q]
    print(f"\n{'='*70}")
    print(f"q={q}, g_c=1/{q}={gc:.4f}, Re(c)={rec}")
    print(f"Sizes: {sizes} (max dim={q**max(sizes):,})")
    print(f"{'='*70}", flush=True)

    q_data = {'q': q, 'gc': gc, 'sizes': [], 'E0_per_N': [], 'E0_raw': [],
              'gaps': [], 'gap_N': []}

    for n in sizes:
        dim = q**n
        print(f"  n={n}, dim={dim:,} ... ", end='', flush=True)

        t0 = time.time()
        H = build_sq_potts_periodic_fast(n, q, gc)
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

        print(f"E0/N={E0_per_N:.10f}, gap={gap:.6f}, gap*N={gapN:.4f}, "
              f"build={t_build:.1f}s, eig={t_eig:.1f}s", flush=True)

    # Fit Casimir formula: E0/N = eps_inf - pi*v*c/(6*N^2)
    N_arr = np.array(q_data['sizes'], dtype=float)
    y = np.array(q_data['E0_per_N'])
    x = 1.0 / N_arr**2

    # 2-parameter fit
    A = np.vstack([x, np.ones_like(x)]).T
    (slope, intercept), _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    eps_inf = intercept
    vc_casimir = -slope * 6 / np.pi

    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    # 3-parameter fit: + 1/N^4 term
    if len(N_arr) >= 4:
        x4 = 1.0 / N_arr**4
        A3 = np.vstack([x, x4, np.ones_like(x)]).T
        (slope3, slope4, intercept3), _, _, _ = np.linalg.lstsq(A3, y, rcond=None)
        vc3 = -slope3 * 6 / np.pi
        y_pred3 = slope3 * x + slope4 * x4 + intercept3
        ss_res3 = np.sum((y - y_pred3)**2)
        R2_3 = 1 - ss_res3 / ss_tot if ss_tot > 0 else 0
    else:
        vc3, R2_3, slope4 = None, None, None

    # Pairwise vc from consecutive sizes
    pairwise = []
    for i in range(len(N_arr) - 1):
        N1, N2 = N_arr[i], N_arr[i+1]
        dE = y[i+1] - y[i]
        dx = 1/N2**2 - 1/N1**2
        vc_pair = -dE / dx * 6 / np.pi
        # Use velocity from the LARGER size for c_implied
        v_pair = q_data['gap_N'][i+1] / (2*np.pi*x_sigma_values[q])
        c_implied_pair = vc_pair / v_pair
        pairwise.append({
            'pair': f'({int(N1)},{int(N2)})',
            'vc': float(vc_pair),
            'v_pair': float(v_pair),
            'c_implied': float(c_implied_pair),
            'c_imp_over_rec': float(c_implied_pair / rec),
        })

    # c_implied from full fit
    v_gap = q_data['gap_N'][-1] / (2 * np.pi * x_sigma_values[q])
    c_implied = vc_casimir / v_gap
    c_imp_over_rec = c_implied / rec

    c_eff = c_eff_values.get(q, None)
    c_eff_over_rec = c_eff / rec if c_eff else None

    print(f"\n  --- 2-parameter Casimir fit ---")
    print(f"  E0/N = {eps_inf:.10f} + ({slope:.8f})/N²")
    print(f"  R² = {R2:.10f}")
    print(f"  v*c = {vc_casimir:.6f}")
    print(f"  v (from gap) = {v_gap:.5f}")
    print(f"  c_implied = vc/v = {c_implied:.4f}")
    print(f"  c_implied/Re(c) = {c_imp_over_rec:.4f}")
    if c_eff:
        print(f"  c_eff/Re(c)     = {c_eff_over_rec:.4f}  (from entropy)")
        print(f"  ** Casimir sees Re(c) {abs(c_imp_over_rec-1)*100:.1f}% off, "
              f"entropy {abs(c_eff_over_rec-1)*100:.1f}% off **")

    if vc3 is not None:
        c_implied3 = vc3 / v_gap
        print(f"\n  --- 3-parameter fit (+ 1/N⁴ term) ---")
        print(f"  v*c = {vc3:.6f}, c_implied = {c_implied3:.4f}, "
              f"c_imp/Re(c) = {c_implied3/rec:.4f}")
        print(f"  R² = {R2_3:.10f}")
        print(f"  1/N⁴ coefficient = {slope4:.6f}")

    print(f"\n  Pairwise v*c from consecutive sizes:")
    for p in pairwise:
        print(f"    {p['pair']}: v*c={p['vc']:.6f}, c_implied={p['c_implied']:.4f}, "
              f"c/Re(c)={p['c_imp_over_rec']:.4f}")

    if len(pairwise) >= 2:
        vc_first = pairwise[0]['vc']
        vc_last = pairwise[-1]['vc']
        drift_pct = (vc_last - vc_first) / vc_first * 100
        print(f"  Pairwise vc drift: {drift_pct:+.2f}% (first→last)")

    print(f"\n  Residuals:")
    for i, n_val in enumerate(q_data['sizes']):
        resid = y[i] - y_pred[i]
        print(f"    n={n_val}: resid={resid:+.4e}")

    q_data['fit_2param'] = {
        'eps_inf': float(eps_inf), 'slope': float(slope),
        'vc_casimir': float(vc_casimir), 'v_gap': float(v_gap),
        'c_implied': float(c_implied), 'c_imp_over_rec': float(c_imp_over_rec),
        'R2': float(R2), 'rec': rec,
    }
    if vc3 is not None:
        q_data['fit_3param'] = {
            'vc': float(vc3), 'c_implied': float(vc3/v_gap),
            'c_imp_over_rec': float(vc3/v_gap/rec),
            'R2': float(R2_3), 'coeff_N4': float(slope4),
        }
    q_data['pairwise'] = pairwise
    if c_eff:
        q_data['c_eff'] = c_eff
        q_data['c_eff_over_rec'] = c_eff_over_rec

    results['data'][str(q)] = q_data
    save()

# Grand summary
print(f"\n{'='*70}")
print("GRAND SUMMARY: c_Casimir vs c_entropy vs Re(c)")
print(f"{'='*70}")
print(f"{'q':>3} {'#pts':>4} {'vc':>8} {'v_gap':>7} {'c_Cas':>6} {'c_Cas/Rec':>9} "
      f"{'c_eff/Rec':>9} {'R²(2p)':>12}")
for q in [2, 3, 5, 7, 8]:
    d = results['data'][str(q)]
    f2 = d['fit_2param']
    npts = len(d['sizes'])
    c_eff_str = f"{d.get('c_eff_over_rec', 0):.4f}" if 'c_eff_over_rec' in d else "N/A"
    print(f"{q:3d} {npts:4d} {f2['vc_casimir']:8.4f} {f2['v_gap']:7.4f} "
          f"{f2['c_implied']:6.4f} {f2['c_imp_over_rec']:9.4f} {c_eff_str:>9} "
          f"{f2['R2']:12.10f}")

# Pairwise convergence summary
print(f"\n{'='*70}")
print("PAIRWISE CONVERGENCE: c_implied/Re(c) from consecutive pairs")
print(f"{'='*70}")
for q in [2, 3, 5, 7, 8]:
    d = results['data'][str(q)]
    pairs = d['pairwise']
    print(f"  q={q}: ", end='')
    for p in pairs:
        print(f"{p['pair']}:{p['c_imp_over_rec']:.3f}  ", end='')
    print()

save()
print(f"\nSaved to results/sprint_098a_casimir_gpu_extended.json")

from db_utils import record
for q in [2, 3, 5, 7, 8]:
    d = results['data'][str(q)]
    f2 = d['fit_2param']
    record(sprint=98, model='sq_potts', q=q, n=max(d['sizes']),
           quantity='c_casimir', value=f2['c_implied'],
           method='casimir_energy_gpu',
           notes=f"vc={f2['vc_casimir']:.6f}, v_gap={f2['v_gap']:.5f}, "
                 f"c/Rec={f2['c_imp_over_rec']:.4f}, R2={f2['R2']:.8f}, "
                 f"{len(d['sizes'])}pts")
    record(sprint=98, model='sq_potts', q=q, n=max(d['sizes']),
           quantity='c_imp_over_rec', value=f2['c_imp_over_rec'],
           method='casimir_energy_gpu',
           notes=f"Casimir c / Re(c), {len(d['sizes'])} sizes")
print("Recorded to DB.")
