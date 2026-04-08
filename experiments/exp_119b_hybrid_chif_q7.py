"""Sprint 119b: chi_F spectral decomposition on HYBRID Potts-clock model at q=7.

Second data point for hybrid alpha(q). Also tests q=3 (where hybrid = S_q)
as a sanity check, and q=10 for wider coverage.

Hybrid g_c values (from Sprints 076-077):
  q=3: 0.333 (= S_q g_c = 1/3, models identical)
  q=7: 0.535
  q=10: 0.684

Size limits (q^n <= 10M):
  q=3: n up to 14 (4.8M), but 14 is ~250s. Use n=4-12.
  q=7: n up to 8 (5.8M). Use n=4-8.
  q=10: n up to 7 (10M). Use n=4-7.
"""
import numpy as np
import json, time, os
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '119b_hybrid_chif_multi_q',
    'sprint': 119,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'model': 'hybrid',
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results', 'sprint_119b_hybrid_chif_multi_q.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)


def build_hybrid_parts_fast(n, q):
    """Build hybrid Potts-clock Hamiltonian: Potts delta coupling + (X+X†) field."""
    dim = q**n
    all_idx = np.arange(dim, dtype=np.int64)
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q
        tmp //= q
    powers = q ** np.arange(n, dtype=np.int64)

    # Coupling: -δ(s_i, s_{i+1}) periodic
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n):
        nxt = (site + 1) % n
        diag_vals -= (digits[:, site] == digits[:, nxt]).astype(np.float64)
    H_coup = csr_matrix((diag_vals, (all_idx, all_idx)), shape=(dim, dim))

    # Field: -(X + X†) at each site
    rows_list, cols_list, vals_list = [], [], []
    for site in range(n):
        pw = powers[site]
        old_digit = digits[:, site]
        for k in [1, q - 1]:  # X (+1) and X† (-1 mod q)
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


def run_chif_spectral(q, g_c, sizes, k_need=12, label=''):
    """Run chi_F spectral decomposition for given q and sizes."""
    print(f"\n{'=' * 70}")
    print(f"HYBRID q={q}, g_c={g_c}, sizes={sizes} {label}")
    print("=" * 70)

    for n in sizes:
        dim = q**n
        if dim > 12_000_000:
            print(f"  n={n} (dim={dim:,}) -- SKIP (too large)")
            continue
        print(f"\n  n={n} (dim={dim:,})...", end="", flush=True)
        t0 = time.time()

        H_coup, H_field = build_hybrid_parts_fast(n, q)
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
        best_level = -1
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
            contributions.append({
                'level': i, 'gap': float(gap),
                'me_sq': float(me**2), 'chi_contrib': float(chi_n),
            })
            if chi_n > best_chi:
                best_chi = chi_n
                best_gap = gap
                best_me_sq = float(me**2)
                best_level = i

        chi_F = total_chi / n
        dominant_frac = best_chi / total_chi if total_chi > 0 else 0
        dt = time.time() - t0

        degen_count = sum(1 for c in contributions
                          if abs(c['gap'] - best_gap) < 0.01 * best_gap)

        print(f" gap_m={best_gap:.6f} (lvl {best_level}, deg~{degen_count}), "
              f"|me|²={best_me_sq:.4e}, chi_F={chi_F:.4f}, "
              f"frac={dominant_frac:.4f} ({dt:.1f}s)")

        contributions.sort(key=lambda x: -x['chi_contrib'])
        for j, c in enumerate(contributions[:3]):
            frac_j = c['chi_contrib'] / total_chi if total_chi > 0 else 0
            print(f"    #{j+1}: lvl {c['level']}, gap={c['gap']:.6f}, "
                  f"|me|²={c['me_sq']:.4e}, frac={frac_j:.4f}")

        key = f"hybrid_q{q}_n{n}"
        results['data'][key] = {
            'n': n, 'dim': dim, 'q': q, 'model': 'hybrid',
            'g_c': g_c,
            'gap_multiplet': float(best_gap),
            'gap_spectral': float(spectral_gap) if spectral_gap else None,
            'dominant_level': best_level,
            'dominant_degeneracy': degen_count,
            'me_sq': best_me_sq,
            'chi_F': float(chi_F),
            'dominant_fraction': float(dominant_frac),
            'E0': float(E0),
            'time_s': round(dt, 1),
            'top3': contributions[:3],
        }
        save()
        record(sprint=119, model='hybrid', q=q, n=n,
               quantity='chi_F', value=float(chi_F), method='spectral_exact_gpu')
        record(sprint=119, model='hybrid', q=q, n=n,
               quantity='gap_multiplet', value=float(best_gap), method='spectral_exact_gpu')


def analyze_scaling(q):
    """Compute pairwise and global alpha for a given q."""
    keys = sorted([k for k in results['data'] if k.startswith(f'hybrid_q{q}_')],
                  key=lambda k: results['data'][k]['n'])
    if len(keys) < 2:
        return

    Ns = np.array([results['data'][k]['n'] for k in keys], dtype=float)
    chis = np.array([results['data'][k]['chi_F'] for k in keys])
    gaps = np.array([results['data'][k]['gap_multiplet'] for k in keys])
    mes = np.array([results['data'][k]['me_sq'] for k in keys])

    log_N = np.log(Ns)
    log_chi = np.log(chis)
    log_gap = np.log(gaps)
    log_me = np.log(mes)

    print(f"\n  --- q={q} pairwise exponents ---")
    print(f"  {'Pair':>8s}  {'alpha':>8s}  {'z_m':>8s}  {'beta_me':>8s}  {'recon':>8s}")
    for i in range(len(Ns)-1):
        dln = log_N[i+1] - log_N[i]
        pw_a = (log_chi[i+1] - log_chi[i]) / dln
        pw_z = -(log_gap[i+1] - log_gap[i]) / dln
        pw_b = (log_me[i+1] - log_me[i]) / dln
        recon = pw_b + 2*pw_z - 1
        print(f"  ({int(Ns[i]):2d},{int(Ns[i+1]):2d})  {pw_a:8.4f}  {pw_z:8.4f}  {pw_b:8.4f}  {recon:8.4f}")

    if len(Ns) >= 3:
        p_chi = np.polyfit(log_N, log_chi, 1)
        p_gap = np.polyfit(log_N, log_gap, 1)
        p_me = np.polyfit(log_N, log_me, 1)
        alpha = p_chi[0]
        z_m = -p_gap[0]
        beta_me = p_me[0]
        nu_eff = 1.0 / (0.5 * (alpha + 1))
        print(f"  Global: alpha={alpha:.4f}, z_m={z_m:.4f}, beta_me={beta_me:.4f}, "
              f"recon={beta_me + 2*z_m - 1:.4f}, nu_eff={nu_eff:.4f}")

        results[f'q{q}_global_alpha'] = float(alpha)
        results[f'q{q}_z_m'] = float(z_m)
        results[f'q{q}_beta_me'] = float(beta_me)
        results[f'q{q}_nu_eff'] = float(nu_eff)


# Run experiments
# q=3 sanity check (hybrid = S_q at q=3)
run_chif_spectral(3, 0.333, [4, 6, 8, 10, 12], k_need=10, label='(sanity: hybrid=S_q)')

# q=7 main target
run_chif_spectral(7, 0.535, [4, 5, 6, 7, 8], k_need=14)

# q=10 extension
run_chif_spectral(10, 0.684, [4, 5, 6, 7], k_need=16)

# Analysis
print(f"\n{'=' * 70}")
print("SCALING ANALYSIS")
print("=" * 70)
analyze_scaling(3)
analyze_scaling(7)
analyze_scaling(10)

# Summary comparison with S_q
print(f"\n{'=' * 70}")
print("MODEL COMPARISON: HYBRID vs S_q")
print("=" * 70)
sq_alpha = {3: 1.40, 5: 2.09, 7: 2.65}  # Known S_q values
for q_val in [3, 7, 10]:
    key = f'q{q_val}_global_alpha'
    if key in results:
        hybrid_a = results[key]
        sq_a = sq_alpha.get(q_val, None)
        sq_str = f"{sq_a:.2f}" if sq_a else "N/A"
        print(f"  q={q_val:2d}: hybrid alpha={hybrid_a:.4f}, S_q alpha={sq_str}")

save()
print("\nAll results saved.")
