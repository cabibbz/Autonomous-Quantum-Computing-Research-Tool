"""Sprint 119a: chi_F spectral decomposition on the HYBRID Potts-clock model at q=5.

The hybrid model uses Potts delta coupling + Z_q clock field (X + X†).
This differs from S_q Potts which uses Σ_{k=1}^{q-1} X^k as the field.

Hybrid q=5: g_c = 0.438 (from Sprint 076), Z_5 symmetry, continuous transition (ν≈0.83).
S_q q=5:   g_c = 0.200 (= 1/q), S_5 symmetry, walking (ν_eff≈0.65).

Key question: does the hybrid show super-scaling (alpha > 2) like S_q,
or ordinary scaling (alpha ≈ 2/ν - 1 ≈ 1.41)?

Sizes: n=4 (625), n=5 (3125), n=6 (15625), n=7 (78125), n=8 (390625),
       n=9 (1953125), n=10 (9765625).
GPU for n>=8. All within 10M dim limit.
"""
import numpy as np
import json, time, os
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '119a_hybrid_chif_q5',
    'sprint': 119,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'model': 'hybrid',
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results', 'sprint_119a_hybrid_chif_q5.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)


def build_hybrid_parts_fast(n, q):
    """Build hybrid Potts-clock Hamiltonian parts: H_coup (Potts delta) and H_field (X+X†).

    H_coup = -Σ_<ij> δ(s_i, s_j)   [Potts coupling, same as S_q]
    H_field = -Σ_i (X_i + X_i†)     [clock field, NOT S_q field]

    X|s> = |(s+1) mod q>, so X†|s> = |(s-1) mod q>.
    X + X† shifts by +1 and -1 only (nearest-neighbor in Z_q space).
    """
    dim = q**n
    all_idx = np.arange(dim, dtype=np.int64)
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q
        tmp //= q
    powers = q ** np.arange(n, dtype=np.int64)

    # Coupling: -δ(s_i, s_{i+1}) for periodic chain
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n):
        nxt = (site + 1) % n
        diag_vals -= (digits[:, site] == digits[:, nxt]).astype(np.float64)
    H_coup = csr_matrix((diag_vals, (all_idx, all_idx)), shape=(dim, dim))

    # Field: -(X + X†) at each site = shift by +1 and shift by -1
    rows_list, cols_list, vals_list = [], [], []
    for site in range(n):
        pw = powers[site]
        old_digit = digits[:, site]
        for k in [1, q - 1]:  # X (shift +1) and X† (shift -1, equiv to +q-1 mod q)
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


q = 5
g_c = 0.438  # Hybrid q=5 critical point (Sprint 076)
all_sizes = [4, 5, 6, 7, 8, 9, 10]
k_need = 12  # need enough states to find dominant multiplet

print("=" * 70)
print(f"Sprint 119a: HYBRID Potts-clock chi_F spectral, q={q}")
print(f"  Field: X + X† (clock, Z_{q} symmetry)")
print(f"  g_c = {g_c} (Sprint 076)")
print(f"  Prediction: alpha ≈ 2/0.83 - 1 ≈ 1.41 (continuous, ν≈0.83)")
print(f"  Comparison: S_q q=5 alpha = 2.09 (walking)")
print("=" * 70)

for n in all_sizes:
    dim = q**n
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

    # Count near-degenerate states around dominant level
    degen_count = sum(1 for c in contributions
                      if abs(c['gap'] - best_gap) < 0.01 * best_gap)

    print(f" gap_m={best_gap:.6f} (level {best_level}, degen~{degen_count}), "
          f"|me|^2={best_me_sq:.4f}, chi_F={chi_F:.4f}, "
          f"frac={dominant_frac:.4f} ({dt:.1f}s)")

    # Show top 5 contributions
    contributions.sort(key=lambda x: -x['chi_contrib'])
    for j, c in enumerate(contributions[:5]):
        frac_j = c['chi_contrib'] / total_chi if total_chi > 0 else 0
        print(f"    #{j+1}: level {c['level']}, gap={c['gap']:.6f}, "
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
        'top5': contributions[:5],
    }
    save()
    record(sprint=119, model='hybrid', q=q, n=n,
           quantity='chi_F', value=float(chi_F), method='spectral_exact_gpu')
    record(sprint=119, model='hybrid', q=q, n=n,
           quantity='gap_multiplet', value=float(best_gap), method='spectral_exact_gpu')

# Pairwise and global analysis
print(f"\n{'=' * 70}")
print("PAIRWISE AND GLOBAL EXPONENTS")
print("=" * 70)

Ns, gaps, mes, chis = [], [], [], []
for n in all_sizes:
    key = f"hybrid_q{q}_n{n}"
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
    recon = pw_b + 2*pw_z - 1
    print(f"  ({int(Ns[i]):2d},{int(Ns[i+1]):2d})  {pw_a:8.4f}  {pw_z:8.4f}  {pw_b:8.4f}  {recon:8.4f}")

# Global fits
if len(Ns) >= 3:
    p_chi = np.polyfit(log_N, log_chi, 1)
    p_gap = np.polyfit(log_N, log_gap, 1)
    p_me = np.polyfit(log_N, log_me, 1)
    alpha_global = p_chi[0]
    z_m_global = -p_gap[0]
    beta_me_global = p_me[0]
    nu_eff = 1.0 / (0.5 * (alpha_global + 1))
    recon_global = beta_me_global + 2*z_m_global - 1

    print(f"\n  GLOBAL (all {len(Ns)} sizes):")
    print(f"    alpha  = {alpha_global:.4f}")
    print(f"    z_m    = {z_m_global:.4f}")
    print(f"    beta_me= {beta_me_global:.4f}")
    print(f"    recon  = {recon_global:.4f} (should = alpha)")
    print(f"    nu_eff = {nu_eff:.4f}")
    print(f"\n  COMPARISON:")
    print(f"    Hybrid prediction: alpha ≈ 1.41 (ν≈0.83)")
    print(f"    S_q q=5 measured:  alpha = 2.09 (walking)")
    print(f"    Ratio S_q/hybrid:  {2.09/alpha_global:.2f}x")

    results['global_alpha'] = float(alpha_global)
    results['global_z_m'] = float(z_m_global)
    results['global_beta_me'] = float(beta_me_global)
    results['nu_eff'] = float(nu_eff)
    save()

print("\nResults saved.")
