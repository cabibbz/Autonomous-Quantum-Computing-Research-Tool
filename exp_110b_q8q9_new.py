"""Sprint 110b: chi_F spectral decomposition for q=8 (n=6,7) and q=9 (n=6,7).

Extends alpha(q) linear formula to higher q. Currently tested at q=5,6,7.
Prediction: alpha(q=8) = 0.315*8 + 0.469 = 2.99, alpha(q=9) = 0.315*9 + 0.469 = 3.30.

Dimensions: q=8 n=7 = 2.1M (GPU), q=9 n=7 = 4.8M (GPU).
Uses vectorized builder (Sprint 098 pattern).
"""
import numpy as np
import json, time, os
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '110b_q8q9_new',
    'sprint': 110,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_110b_q8q9_new.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts_fast(n, q):
    """Build S_q Potts coupling and field parts separately — vectorized."""
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
    H_field = coo_matrix((vals, (rows, cols)), shape=(dim, dim)).tocsr()

    return H_coup, H_field

configs = [
    (8, [6, 7], 10),   # q=8: k=10 (7-fold multiplet + singlet + margin)
    (9, [6, 7], 12),   # q=9: k=12 (8-fold multiplet + singlet + margin)
]

print("=" * 70)
print("Sprint 110b: q=8,9 chi_F spectral decomposition")
print("=" * 70)

for q, sizes, k_need in configs:
    g_c = 1.0 / q
    print(f"\n--- q={q}, g_c={g_c:.6f}, k={k_need} ---")

    for n in sizes:
        dim = q**n
        print(f"\n  n={n} (dim={dim:,})...", end="", flush=True)
        t0 = time.time()

        H_coup, H_field = build_sq_potts_parts_fast(n, q)
        t_build = time.time() - t0
        print(f" [build {t_build:.1f}s]", end="", flush=True)

        H = H_coup + g_c * H_field
        k_use = min(k_need, dim - 2)
        evals, evecs = eigsh(H, k=k_use, which='SA')
        order = np.argsort(evals)
        evals = evals[order]
        evecs = evecs[:, order]

        E0 = evals[0]
        psi0 = evecs[:, 0]
        H_field_psi0 = H_field.dot(psi0)

        best_chi = 0
        best_gap = 0
        best_me_sq = 0
        total_chi = 0
        spectral_gap = None
        contributions = []
        for i in range(1, k_use):
            gap = evals[i] - E0
            if gap < 1e-12:
                continue
            if spectral_gap is None:
                spectral_gap = gap
            me = np.dot(evecs[:, i], H_field_psi0)
            chi_n = me**2 / gap**2
            total_chi += chi_n
            contributions.append({'level': i, 'gap': float(gap), 'me_sq': float(me**2),
                                  'chi_contrib': float(chi_n)})
            if chi_n > best_chi:
                best_chi = chi_n
                best_gap = gap
                best_me_sq = float(me**2)

        chi_F = total_chi / n
        dominant_frac = best_chi / total_chi if total_chi > 0 else 0
        dt = time.time() - t0
        print(f" gap_m={best_gap:.6f}, |me|²={best_me_sq:.4f}, chi_F={chi_F:.4f}, "
              f"frac={dominant_frac:.4f} ({dt:.1f}s)")

        key = f"q{q}_n{n}"
        results['data'][key] = {
            'n': n, 'dim': dim, 'q': q,
            'gap_multiplet': float(best_gap),
            'gap_spectral': float(spectral_gap) if spectral_gap else None,
            'me_sq': best_me_sq,
            'chi_F': float(chi_F),
            'dominant_fraction': float(dominant_frac),
            'E0': float(E0),
            'time_s': round(dt, 1),
            'contributions': contributions,
        }
        save()

        record(sprint=110, model='sq_potts', q=q, n=n,
               quantity='chi_F', value=float(chi_F), method='spectral_exact_gpu')
        record(sprint=110, model='sq_potts', q=q, n=n,
               quantity='gap_multiplet', value=float(best_gap), method='spectral_exact_gpu')
        record(sprint=110, model='sq_potts', q=q, n=n,
               quantity='me_sq', value=best_me_sq, method='spectral_exact_gpu')

# Pairwise for each q
print(f"\n{'=' * 70}")
print("PAIRWISE EXPONENTS")
print("=" * 70)

for q, sizes, _ in configs:
    Ns, gaps, mes, chis = [], [], [], []
    for n in sizes:
        key = f"q{q}_n{n}"
        d = results['data'].get(key)
        if d:
            Ns.append(d['n'])
            gaps.append(d['gap_multiplet'])
            mes.append(d['me_sq'])
            chis.append(d['chi_F'])

    if len(Ns) < 2:
        continue

    Ns = np.array(Ns, dtype=float)
    log_N = np.log(Ns)
    log_gap = np.log(np.array(gaps))
    log_me = np.log(np.array(mes))
    log_chi = np.log(np.array(chis))

    dln = log_N[1] - log_N[0]
    z_m = -(log_gap[1] - log_gap[0]) / dln
    beta_me = (log_me[1] - log_me[0]) / dln
    alpha = (log_chi[1] - log_chi[0]) / dln
    predicted_alpha = beta_me + 2 * z_m - 1
    formula_alpha = 0.315 * q + 0.469

    print(f"\n  q={q}: n={sizes}")
    print(f"    z_m = {z_m:.4f}")
    print(f"    beta_me = {beta_me:.4f}")
    print(f"    alpha = {alpha:.4f}")
    print(f"    beta+2z-1 = {predicted_alpha:.4f} (diff={alpha-predicted_alpha:.6f})")
    print(f"    Formula prediction: {formula_alpha:.3f}")
    print(f"    Measured/predicted: {alpha/formula_alpha:.4f}")

    results[f'pairwise_q{q}'] = {
        'z_m': float(z_m),
        'beta_me': float(beta_me),
        'alpha': float(alpha),
        'predicted_alpha': float(predicted_alpha),
        'formula_alpha': float(formula_alpha),
    }

save()
print("\nResults saved.")
