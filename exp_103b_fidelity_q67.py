#!/usr/bin/env python3
"""Sprint 103b: Fidelity susceptibility at q=6 (n=6,8) and q=7 (n=6,7,8) to map alpha(q).

Maps alpha(q) across the walking boundary: q=5 (walking), q=6 (marginal), q=7 (broken).
q=7 n=8 = 5.76M dim — GPU needed, ~300s.
"""
import numpy as np
import json, time
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh

results = {
    'experiment': '103b_fidelity_q67',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'method': 'chi_F at g_c=1/q, dg=1e-4, S_q Potts periodic',
    'data': {},
}

def save():
    with open("results/sprint_103b_fidelity_q67.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts(n, q):
    """Vectorized coupling/field split: H(g) = H_coup + g * H_field."""
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
    H_coup = csr_matrix((diag_vals, (all_idx, all_idx)), shape=(dim, dim))

    # Off-diagonal: field
    rows_list, cols_list, vals_list = [], [], []
    for site in range(n):
        pw = powers[site]
        old_digit = digits[:, site]
        for k in range(1, q):
            new_digit = (old_digit + k) % q
            delta = (new_digit.astype(np.int64) - old_digit.astype(np.int64)) * pw
            new_idx = all_idx + delta
            rows_list.append(new_idx)
            cols_list.append(all_idx)
            vals_list.append(np.full(dim, -1.0, dtype=np.float64))

    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    vals = np.concatenate(vals_list)
    H_field = csr_matrix(coo_matrix((vals, (rows, cols)), shape=(dim, dim)))

    return H_coup, H_field

dg = 1e-4

# Configurations: q -> list of n values
configs = {
    6: {'sizes': [6, 8], 'g_c': 1/6},
    7: {'sizes': [6, 7, 8], 'g_c': 1/7},
}

print("=" * 60)
print("Sprint 103b: chi_F at q=6,7 — mapping alpha(q)")
print("=" * 60)

for q, cfg in configs.items():
    g_c = cfg['g_c']
    results['data'][str(q)] = {'g_c': g_c, 'sizes': {}}

    for n in cfg['sizes']:
        dim = q**n
        print(f"\nq={q}, n={n} (dim={dim:,})...", flush=True)
        t0 = time.time()

        H_coup, H_field = build_sq_potts_parts(n, q)
        t_build = time.time() - t0
        print(f"  Built H parts in {t_build:.1f}s", flush=True)

        # Ground state at g_c
        H = H_coup + g_c * H_field
        t1 = time.time()
        _, vecs1 = eigsh(H, k=1, which='SA')
        psi1 = vecs1[:, 0]
        t_eig1 = time.time() - t1
        print(f"  eigsh #1: {t_eig1:.1f}s", flush=True)

        # Ground state at g_c + dg
        H2 = H_coup + (g_c + dg) * H_field
        t2 = time.time()
        _, vecs2 = eigsh(H2, k=1, which='SA')
        psi2 = vecs2[:, 0]
        t_eig2 = time.time() - t2
        print(f"  eigsh #2: {t_eig2:.1f}s", flush=True)

        overlap = abs(np.vdot(psi1, psi2))
        overlap = min(overlap, 1.0 - 1e-15)
        chi_F = 2.0 * (1.0 - overlap) / (dg**2 * n)

        dt = time.time() - t0
        print(f"  chi_F = {chi_F:.4f}, overlap = {overlap:.15f}")
        print(f"  Total: {dt:.1f}s")

        results['data'][str(q)]['sizes'][str(n)] = {
            'dim': dim, 'chi_F': float(chi_F),
            'overlap': float(overlap), 'time_s': dt,
        }
        save()

        # Time guard for q=7 n=8
        if dt > 250:
            print(f"  WARNING: {dt:.0f}s — skipping larger sizes")
            break

# Include Sprint 102 data for combined analysis
sprint102_data = {
    2: {6: 0.73, 8: 0.95, 10: 1.19, 12: 1.43, 14: 1.67},
    3: {6: 4.09, 8: 6.12, 10: 8.27},
    4: {6: 13.98, 8: 22.75},
    5: {6: 36.29, 8: 66.18},
    7: {6: 166.88},
}

print("\n" + "=" * 60)
print("Combined FSS analysis (Sprint 102 + 103b)")
print("=" * 60)

# Merge q=7 n=6 from Sprint 102
for q_str, qdata in results['data'].items():
    q = int(q_str)
    if q in sprint102_data:
        for n_str in list(qdata['sizes'].keys()):
            n = int(n_str)
            if n in sprint102_data[q]:
                print(f"  q={q} n={n}: 102={sprint102_data[q][n]:.2f}, 103={qdata['sizes'][n_str]['chi_F']:.2f}")

# Alpha fits for q=6 and q=7
for q_str, qdata in results['data'].items():
    q = int(q_str)
    # Combine with Sprint 102 data
    all_n = []
    all_chi = []
    if q in sprint102_data:
        for n, chi in sprint102_data[q].items():
            if str(n) not in qdata['sizes']:
                all_n.append(n)
                all_chi.append(chi)
    for n_str, d in qdata['sizes'].items():
        all_n.append(int(n_str))
        all_chi.append(d['chi_F'])

    all_n = np.array(sorted(all_n))
    all_chi = np.array([c for _, c in sorted(zip(all_n, all_chi))])

    if len(all_n) >= 2:
        log_n = np.log(all_n)
        log_chi = np.log(all_chi)
        coeffs = np.polyfit(log_n, log_chi, 1)
        alpha = coeffs[0]
        nu = 2 / (alpha + 1)
        print(f"\n  q={q}: alpha = {alpha:.4f}, nu_eff = {nu:.4f}")
        for i in range(len(all_n)):
            print(f"    n={all_n[i]}: chi_F = {all_chi[i]:.4f}")
        for i in range(len(all_n)-1):
            a = (np.log(all_chi[i+1]) - np.log(all_chi[i])) / (np.log(all_n[i+1]) - np.log(all_n[i]))
            print(f"    pairwise ({all_n[i]},{all_n[i+1]}): alpha = {a:.4f}")

        results['data'][q_str]['alpha'] = float(alpha)
        results['data'][q_str]['nu_eff'] = float(nu)

save()
print("\nResults saved to results/sprint_103b_fidelity_q67.json")
