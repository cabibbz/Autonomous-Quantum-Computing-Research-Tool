"""Sprint 118a: q=4 chi_F extended to more system sizes.

Sprint 102c had only n=6,8 -> alpha=1.69, nu=0.743.
Exact q=4 Potts: nu=2/3, so alpha=2/nu-1=2.0.
Can we resolve the discrepancy with more sizes?

q=4: n=4 (256), n=5 (1024), n=6 (4096), n=7 (16384),
      n=8 (65536), n=9 (262144), n=10 (1048576), n=11 (4194304).
All feasible on GPU. Use spectral decomposition like q=5+ experiments.
"""
import numpy as np
import json, time, os
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '118a_q4_chif_extend',
    'sprint': 118,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results', 'sprint_118a_q4_chif.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts_fast(n, q):
    """Build S_q Potts coupling and field parts separately -- vectorized."""
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
    H_field = coo_matrix((vals, (rows, cols)), shape=(dim, dim)).tocsr()
    return H_coup, H_field

q = 4
g_c = 1.0 / q  # = 0.25, self-dual point
# n=4 to n=11: dims 256 to 4.2M (all GPU feasible)
all_sizes = [4, 5, 6, 7, 8, 9, 10, 11]
k_need = 10  # (q-1)=3-fold multiplet + margin

print("=" * 70)
print(f"Sprint 118a: q={q} chi_F spectral decomposition, n={all_sizes}")
print(f"  g_c = {g_c:.6f}")
print(f"  Exact q=4 Potts: nu=2/3, predicted alpha=2/nu-1=2.0")
print(f"  Sprint 102c: alpha=1.69 (n=6,8 only)")
print("=" * 70)

for n in all_sizes:
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
    print(f" gap_m={best_gap:.6f}, |me|^2={best_me_sq:.4f}, chi_F={chi_F:.4f}, "
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
    }
    save()
    record(sprint=118, model='sq_potts', q=q, n=n,
           quantity='chi_F', value=float(chi_F), method='spectral_exact_gpu')
    record(sprint=118, model='sq_potts', q=q, n=n,
           quantity='gap_multiplet', value=float(best_gap), method='spectral_exact_gpu')

# Pairwise and global analysis
print(f"\n{'=' * 70}")
print("PAIRWISE AND GLOBAL EXPONENTS")
print("=" * 70)

Ns, gaps, mes, chis = [], [], [], []
for n in all_sizes:
    key = f"q{q}_n{n}"
    d = results['data'].get(key)
    if d:
        Ns.append(d['n'])
        gaps.append(d['gap_multiplet'])
        mes.append(d['me_sq'])
        chis.append(d['chi_F'])

Ns = np.array(Ns, dtype=float)
gaps = np.array(gaps)
mes = np.array(mes)
chis = np.array(chis)

log_N = np.log(Ns)
log_chi = np.log(chis)
log_gap = np.log(gaps)
log_me = np.log(mes)

# Pairwise exponents
print(f"\n  {'Pair':>8s}  {'alpha':>8s}  {'z_m':>8s}  {'beta_me':>8s}  {'recon':>8s}")
for i in range(len(Ns)-1):
    dln = log_N[i+1] - log_N[i]
    pw_a = (log_chi[i+1] - log_chi[i]) / dln
    pw_z = -(log_gap[i+1] - log_gap[i]) / dln
    pw_b = (log_me[i+1] - log_me[i]) / dln
    print(f"  ({int(Ns[i]):2d},{int(Ns[i+1]):2d})  {pw_a:8.4f}  {pw_z:8.4f}  {pw_b:8.4f}  {pw_b + 2*pw_z - 1:8.4f}")

# Progressive global fits
print(f"\n  Progressive global alpha (using sizes from n_min to n_max):")
for start in range(len(Ns)-1):
    for end in range(start+2, len(Ns)+1):
        subset_N = Ns[start:end]
        subset_chi = chis[start:end]
        if len(subset_N) >= 2:
            p = np.polyfit(np.log(subset_N), np.log(subset_chi), 1)
            print(f"    n={int(subset_N[0])}-{int(subset_N[-1])} ({len(subset_N)} pts): alpha={p[0]:.4f}")

# Overall fit
p_all = np.polyfit(log_N, log_chi, 1)
alpha_global = p_all[0]
nu_eff = 1.0 / (0.5 * (alpha_global + 1))

print(f"\n  GLOBAL alpha (all {len(Ns)} sizes): {alpha_global:.4f}")
print(f"  Effective nu = 1/(0.5*(alpha+1)) = {nu_eff:.4f}")
print(f"  Exact q=4: alpha=2.0, nu=2/3=0.667")
print(f"  Our logarithmic formula: alpha = 1.87*ln(4) - 0.97 = {1.87*np.log(4) - 0.97:.3f}")

# Check if BKT-like (alpha increasing with N)
print(f"\n  Is alpha INCREASING (BKT) or CONVERGING (power-law)?")
if len(Ns) >= 4:
    first_half = np.polyfit(log_N[:len(Ns)//2+1], log_chi[:len(Ns)//2+1], 1)[0]
    second_half = np.polyfit(log_N[len(Ns)//2:], log_chi[len(Ns)//2:], 1)[0]
    print(f"    First half (n={int(Ns[0])}-{int(Ns[len(Ns)//2])}):  alpha={first_half:.4f}")
    print(f"    Second half (n={int(Ns[len(Ns)//2])}-{int(Ns[-1])}): alpha={second_half:.4f}")
    if second_half > first_half + 0.05:
        print(f"    -> INCREASING: suggests BKT (alpha->infinity) or walking")
    elif abs(second_half - first_half) < 0.05:
        print(f"    -> STABLE: suggests true power-law")
    else:
        print(f"    -> DECREASING: suggests convergence to smaller value")

results['global_alpha'] = float(alpha_global)
results['nu_eff'] = float(nu_eff)
save()
print("\nResults saved.")
