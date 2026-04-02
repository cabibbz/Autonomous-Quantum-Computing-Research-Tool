"""Sprint 107a (part 2): q=5 n=10 and q=7 n=8 with vectorized builder + GPU.

The loop-based builder is too slow for dim>1M. Use vectorized builder from Sprint 103.
q=5 n=10: dim=9,765,625 — build ~20s vectorized, eigsh ~160s GPU
q=7 n=8: dim=5,764,801 — build ~12s vectorized, eigsh ~90s GPU
"""
import numpy as np
import json, time, os
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '107a2_gpu_extend',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_107a2_gpu_extend.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts_fast(n, q):
    """Vectorized S_q Potts Hamiltonian builder — fast for large dim."""
    dim = q**n
    all_idx = np.arange(dim, dtype=np.int64)
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q
        tmp //= q
    powers = q ** np.arange(n, dtype=np.int64)

    # Diagonal coupling
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n):
        nxt = (site + 1) % n
        diag_vals -= (digits[:, site] == digits[:, nxt]).astype(np.float64)
    H_coup = csr_matrix((diag_vals, (all_idx, all_idx)), shape=(dim, dim))

    # Off-diagonal field
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
    (5, 10),  # dim=9,765,625
    (7, 8),   # dim=5,764,801
]

print("=" * 70)
print("Sprint 107a2: GPU spectral decomposition for large sizes")
print("=" * 70)

for q, n in configs:
    dim = q**n
    g_c = 1.0 / q
    k_need = q + 2  # enough for multiplet at level q-1

    print(f"\nq={q}, n={n} (dim={dim:,})...", flush=True)
    t0 = time.time()

    H_coup, H_field = build_sq_potts_parts_fast(n, q)
    t_build = time.time() - t0
    print(f"  Vectorized build: {t_build:.1f}s", flush=True)

    H = H_coup + g_c * H_field
    k_use = min(k_need, dim - 2)
    evals, evecs = eigsh(H, k=k_use, which='SA')
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]
    t_eig = time.time() - t0 - t_build
    print(f"  Eigsh (GPU): {t_eig:.1f}s (k={k_use})", flush=True)

    E0 = evals[0]
    psi0 = evecs[:, 0]
    H_field_psi0 = H_field.dot(psi0)
    gap_spectral = evals[1] - E0

    best_chi = 0
    best_gap = 0
    best_me_sq = 0
    best_level = -1
    total_chi = 0

    for i in range(1, k_use):
        gap = evals[i] - E0
        if gap < 1e-12:
            continue
        me = np.dot(evecs[:, i], H_field_psi0)
        chi_n = me**2 / gap**2
        total_chi += chi_n
        if chi_n > best_chi:
            best_chi = chi_n
            best_gap = gap
            best_me_sq = me**2
            best_level = i

    chi_F = total_chi / n
    dom_frac = best_chi / total_chi if total_chi > 0 else 0
    dt = time.time() - t0

    results['data'][f'q{q}_n{n}'] = {
        'q': q, 'n': n, 'dim': dim,
        'chi_F': float(chi_F),
        'gap_spectral': float(gap_spectral),
        'gap_multiplet': float(best_gap),
        'me_sq': float(best_me_sq),
        'dominant_level': best_level,
        'dominant_fraction': float(dom_frac),
        'time_build_s': t_build,
        'time_eigsh_s': t_eig,
        'time_total_s': dt,
    }

    record(sprint=107, model='S_q_Potts', q=q, n=n,
           quantity='gap_multiplet', value=best_gap, method='spectral_decomp_gpu')
    record(sprint=107, model='S_q_Potts', q=q, n=n,
           quantity='me_sq_multiplet', value=best_me_sq, method='spectral_decomp_gpu')
    record(sprint=107, model='S_q_Potts', q=q, n=n,
           quantity='chi_F', value=chi_F, method='spectral_decomp_gpu')

    print(f"  gap_spec={gap_spectral:.6f}  gap_mult={best_gap:.6f}")
    print(f"  |me|²={best_me_sq:.4f}  chi_F={chi_F:.4f}")
    print(f"  frac={dom_frac:.4f}  level={best_level}  total={dt:.1f}s")

    save()

save()
print("\nAll done. Results saved.")
