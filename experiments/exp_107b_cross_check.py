"""Sprint 107b: Cross-check spectral chi_F vs finite-difference chi_F.

Two independent methods must agree:
1. Spectral: chi_F = (1/N) sum_n |<n|H_field|0>|^2 / (E_n-E0)^2
2. Finite-diff: chi_F = 2(1 - |<psi(g)|psi(g+dg)>|) / (dg^2 * N)

Compare at q=2,3,5 for sizes where both are fast (n<=8).
Also verify Sprint 103 finite-diff values match Sprint 106 spectral values.
"""
import numpy as np
import json, time, os
from scipy.sparse import lil_matrix, csr_matrix, coo_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '107b_cross_check',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_107b_cross_check.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts(n, q):
    """Vectorized S_q Potts Hamiltonian builder."""
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
    H_coup = csr_matrix((diag_vals, (all_idx, all_idx)), shape=(dim, dim))

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

configs = [
    (2, [6, 8, 10, 12]),
    (3, [6, 8]),
    (5, [6, 8]),
]

dg = 1e-4  # Same as Sprint 103

print("=" * 70)
print("Sprint 107b: Cross-check spectral vs finite-difference chi_F")
print("=" * 70)

for q, sizes in configs:
    g_c = 1.0 / q
    results['data'][f'q{q}'] = {'q': q, 'g_c': g_c, 'sizes': {}}

    for n in sizes:
        dim = q**n
        print(f"\nq={q}, n={n} (dim={dim:,})...", flush=True)
        t0 = time.time()

        H_coup, H_field = build_sq_potts_parts(n, q)

        # === Method 1: Spectral decomposition ===
        H = H_coup + g_c * H_field
        k_need = min(q + 2, dim - 2)  # enough for (q-1)-fold multiplet
        evals, evecs = eigsh(H, k=k_need, which='SA')
        order = np.argsort(evals)
        evals = evals[order]
        evecs = evecs[:, order]

        E0 = evals[0]
        psi0 = evecs[:, 0]
        H_field_psi0 = H_field.dot(psi0)

        total_chi_spec = 0
        for i in range(1, k_need):
            gap = evals[i] - E0
            if gap < 1e-12:
                continue
            me = np.dot(evecs[:, i], H_field_psi0)
            total_chi_spec += me**2 / gap**2
        chi_F_spec = total_chi_spec / n

        # === Method 2: Finite difference ===
        # Only need ground states at g and g+dg
        _, vecs_g = eigsh(H, k=1, which='SA')
        psi_g = vecs_g[:, 0]

        H2 = H_coup + (g_c + dg) * H_field
        _, vecs_g2 = eigsh(H2, k=1, which='SA')
        psi_g2 = vecs_g2[:, 0]

        overlap = abs(np.vdot(psi_g, psi_g2))
        overlap = min(overlap, 1.0 - 1e-15)
        chi_F_fd = 2.0 * (1.0 - overlap) / (dg**2 * n)

        # Relative difference
        rel_diff = abs(chi_F_spec - chi_F_fd) / chi_F_fd if chi_F_fd > 0 else 0

        dt = time.time() - t0

        results['data'][f'q{q}']['sizes'][str(n)] = {
            'dim': dim,
            'chi_F_spectral': float(chi_F_spec),
            'chi_F_finite_diff': float(chi_F_fd),
            'relative_diff': float(rel_diff),
            'k_states': k_need,
            'time_s': dt,
        }

        record(sprint=107, model='S_q_Potts', q=q, n=n,
               quantity='chi_F_fd', value=chi_F_fd, method='finite_diff')
        record(sprint=107, model='S_q_Potts', q=q, n=n,
               quantity='chi_F_spec', value=chi_F_spec, method='spectral_k_states',
               notes=f'k={k_need}')

        status = "MATCH" if rel_diff < 0.01 else "MISMATCH" if rel_diff < 0.1 else "FAIL"
        print(f"  Spectral:  chi_F = {chi_F_spec:.6f}")
        print(f"  Finite-diff: chi_F = {chi_F_fd:.6f}")
        print(f"  Rel diff: {rel_diff:.2e}  [{status}]  ({dt:.1f}s)")

        save()

# Summary
print(f"\n{'=' * 70}")
print("SUMMARY: Spectral vs Finite-Difference chi_F")
print(f"{'q':>3} {'n':>4} {'chi_spec':>12} {'chi_fd':>12} {'rel_diff':>10} {'status':>8}")
print("-" * 55)
for q, sizes in configs:
    for n in sizes:
        d = results['data'][f'q{q}']['sizes'].get(str(n), {})
        if d:
            cs = d['chi_F_spectral']
            cf = d['chi_F_finite_diff']
            rd = d['relative_diff']
            st = "MATCH" if rd < 0.01 else "MISMATCH" if rd < 0.1 else "FAIL"
            print(f"{q:>3} {n:>4} {cs:>12.6f} {cf:>12.6f} {rd:>10.2e} {st:>8}")

save()
print("\nResults saved.")
