#!/usr/bin/env python3
"""Sprint 103a: Fidelity susceptibility at q=5 for n=6,8,9,10 (GPU for large sizes).

Uses vectorized Hamiltonian builder (from Sprint 098) for n=9,10.
Single-point evaluation at g_c=0.2 (zero FSS shift for q>=5, Sprint 102).
"""
import numpy as np
import json, time
from scipy.sparse import coo_matrix, csr_matrix, lil_matrix
from gpu_utils import eigsh

results = {
    'experiment': '103a_fidelity_q5_gpu',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'method': 'chi_F at g_c=1/q, dg=1e-4, S_q Potts periodic',
    'data': {},
}

def save():
    with open("results/sprint_103a_fidelity_q5_gpu.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts_fast(n, q):
    """Vectorized coupling/field split: H(g) = H_coup + g * H_field."""
    dim = q**n
    all_idx = np.arange(dim, dtype=np.int64)

    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q
        tmp //= q

    powers = q ** np.arange(n, dtype=np.int64)

    # Diagonal: coupling = -sum delta(s_i, s_{i+1})
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n):
        nxt = (site + 1) % n
        diag_vals -= (digits[:, site] == digits[:, nxt]).astype(np.float64)

    H_coup = csr_matrix((diag_vals, (all_idx, all_idx)), shape=(dim, dim))

    # Off-diagonal: field = -sum_site sum_{k=1}^{q-1} |shifted><original|
    rows_list = []
    cols_list = []
    vals_list = []
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

def build_sq_potts_parts_slow(n, q):
    """Loop-based builder for small sizes."""
    dim = q**n
    H_coup = lil_matrix((dim, dim), dtype=float)
    H_field = lil_matrix((dim, dim), dtype=float)

    for idx in range(dim):
        c = []
        tmp = idx
        for _ in range(n):
            c.append(tmp % q)
            tmp //= q

        diag = 0.0
        for site in range(n):
            nxt = (site + 1) % n
            if c[site] == c[nxt]:
                diag -= 1.0
        H_coup[idx, idx] = diag

        for site in range(n):
            for k in range(1, q):
                c2 = c.copy()
                c2[site] = (c[site] + k) % q
                idx2 = 0
                for i in range(n - 1, -1, -1):
                    idx2 = idx2 * q + c2[i]
                H_field[idx, idx2] += -1.0

    return csr_matrix(H_coup), csr_matrix(H_field)

q = 5
g_c = 0.2
dg = 1e-4

print("=" * 60)
print("Sprint 103a: chi_F at q=5, extending to n=9,10 (GPU)")
print("=" * 60)

for n in [6, 8, 9, 10]:
    dim = q**n
    print(f"\nq={q}, n={n} (dim={dim:,})...", flush=True)
    t0 = time.time()

    # Use fast builder for large sizes
    if dim > 100000:
        H_coup, H_field = build_sq_potts_parts_fast(n, q)
    else:
        H_coup, H_field = build_sq_potts_parts_slow(n, q)
    t_build = time.time() - t0
    print(f"  Built H parts in {t_build:.1f}s", flush=True)

    # Ground state at g_c
    H = H_coup + g_c * H_field
    t1 = time.time()
    _, vecs1 = eigsh(H, k=1, which='SA')
    psi1 = vecs1[:, 0]
    t_eig1 = time.time() - t1
    print(f"  eigsh #1 (g_c): {t_eig1:.1f}s", flush=True)

    # Ground state at g_c + dg
    H2 = H_coup + (g_c + dg) * H_field
    t2 = time.time()
    _, vecs2 = eigsh(H2, k=1, which='SA')
    psi2 = vecs2[:, 0]
    t_eig2 = time.time() - t2
    print(f"  eigsh #2 (g_c+dg): {t_eig2:.1f}s", flush=True)

    overlap = abs(np.vdot(psi1, psi2))
    overlap = min(overlap, 1.0 - 1e-15)
    chi_F = 2.0 * (1.0 - overlap) / (dg**2 * n)

    dt = time.time() - t0
    print(f"  chi_F = {chi_F:.4f}")
    print(f"  overlap = {overlap:.15f}")
    print(f"  Total time: {dt:.1f}s")

    results['data'][str(n)] = {
        'q': q, 'n': n, 'dim': dim, 'g_c': g_c, 'dg': dg,
        'chi_F': float(chi_F), 'overlap': float(overlap),
        'time_build_s': t_build, 'time_eig1_s': t_eig1,
        'time_eig2_s': t_eig2, 'time_total_s': dt,
    }
    save()
    print("  Saved.", flush=True)

# FSS analysis
print("\n" + "=" * 60)
print("FSS: chi_F ~ N^alpha (log-log fit)")
print("=" * 60)

ns = []
chis = []
for n_str, d in sorted(results['data'].items(), key=lambda x: int(x[0])):
    n_val = int(n_str)
    ns.append(n_val)
    chis.append(d['chi_F'])
    print(f"  n={n_val}: chi_F = {d['chi_F']:.4f}")

ns = np.array(ns)
chis = np.array(chis)

if len(ns) >= 2:
    log_n = np.log(ns)
    log_chi = np.log(chis)
    coeffs = np.polyfit(log_n, log_chi, 1)
    alpha_full = coeffs[0]
    print(f"\n  Full fit ({len(ns)} points): alpha = {alpha_full:.4f}")
    print(f"  -> nu_eff = 2/(alpha+1) = {2/(alpha_full+1):.4f}")

    print("\n  Pairwise alpha:")
    pairwise = []
    for i in range(len(ns)-1):
        a = (np.log(chis[i+1]) - np.log(chis[i])) / (np.log(ns[i+1]) - np.log(ns[i]))
        pairwise.append(float(a))
        print(f"    ({ns[i]},{ns[i+1]}): alpha = {a:.4f} -> nu = {2/(a+1):.4f}")

    results['fss'] = {
        'alpha_full': float(alpha_full),
        'nu_eff': float(2/(alpha_full+1)),
        'pairwise_alpha': pairwise,
        'sizes': ns.tolist(),
        'chi_F_values': chis.tolist(),
    }
    save()

print("\nResults saved to results/sprint_103a_fidelity_q5_gpu.json")
