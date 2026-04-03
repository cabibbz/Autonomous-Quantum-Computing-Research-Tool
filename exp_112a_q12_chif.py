"""Sprint 112a: chi_F spectral decomposition for q=12 at n=4,5,6.

Tests power-law vs linear alpha(q) at q=12 where models diverge:
  Linear  (0.260q+0.827):    alpha = 3.947
  Power-law (0.69*q^0.69):   alpha = 4.502
  Difference: 0.55

Dimensions: n=4 (21k), n=5 (249k), n=6 (3.0M GPU).
"""
import numpy as np
import json, time, os
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '112a_q12_chif',
    'sprint': 112,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_112a_q12_chif.json')
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

    # Diagonal coupling: -Σ δ(s_i, s_{i+1})
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n):
        nxt = (site + 1) % n
        diag_vals -= (digits[:, site] == digits[:, nxt]).astype(np.float64)

    H_coup = csr_matrix((diag_vals, (all_idx, all_idx)), shape=(dim, dim))

    # Off-diagonal field: -Σ_{site} Σ_{k=1}^{q-1} |shifted><original|
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

q = 12
g_c = 1.0 / q
all_sizes = [4, 5, 6]
k_need = 14  # (q-1)=11-fold multiplet + singlet + margin

print("=" * 70)
print(f"Sprint 112a: q={q} chi_F spectral decomposition, n={all_sizes}")
print(f"  g_c = {g_c:.6f}, k = {k_need}")
print(f"  Linear prediction:    alpha = {0.260*q + 0.827:.3f}")
print(f"  Power-law prediction: alpha = {0.69 * q**0.69:.3f}")
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

    # Spectral decomposition of chi_F
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

    record(sprint=112, model='sq_potts', q=q, n=n,
           quantity='chi_F', value=float(chi_F), method='spectral_exact_gpu')
    record(sprint=112, model='sq_potts', q=q, n=n,
           quantity='gap_multiplet', value=float(best_gap), method='spectral_exact_gpu')
    record(sprint=112, model='sq_potts', q=q, n=n,
           quantity='me_sq', value=best_me_sq, method='spectral_exact_gpu')

# Pairwise analysis
print(f"\n{'=' * 70}")
print("PAIRWISE EXPONENTS")
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
log_gap = np.log(gaps)
log_me = np.log(mes)
log_chi = np.log(chis)

pw_z = []; pw_b = []; pw_a = []; pw_N = []
for i in range(len(Ns)-1):
    dln = log_N[i+1] - log_N[i]
    pw_z.append(-(log_gap[i+1] - log_gap[i]) / dln)
    pw_b.append((log_me[i+1] - log_me[i]) / dln)
    pw_a.append((log_chi[i+1] - log_chi[i]) / dln)
    pw_N.append(0.5 * (Ns[i] + Ns[i+1]))

print(f"  Sizes: {[int(n) for n in Ns]}")
print(f"  Pairwise z_m:     {['%.4f' % z for z in pw_z]}")
print(f"  Pairwise beta_me: {['%.4f' % b for b in pw_b]}")
print(f"  Pairwise alpha:   {['%.4f' % a for a in pw_a]}")

# Global log-log fit
if len(Ns) >= 2:
    p_a = np.polyfit(log_N, log_chi, 1)
    alpha_global = p_a[0]
    print(f"\n  Global alpha (all {len(Ns)} sizes): {alpha_global:.4f}")

    # Decomposition check
    for i in range(len(pw_z)):
        predicted = pw_b[i] + 2 * pw_z[i] - 1
        print(f"  Pair {int(Ns[i])},{int(Ns[i+1])}: alpha={pw_a[i]:.4f}, "
              f"beta+2z-1={predicted:.4f}, diff={pw_a[i]-predicted:.6f}")

    # Component global fits
    p_z = np.polyfit(log_N, log_gap, 1)
    z_m_global = -p_z[0]
    p_b = np.polyfit(log_N, log_me, 1)
    beta_me_global = p_b[0]
    print(f"\n  Global z_m:     {z_m_global:.4f}")
    print(f"  Global beta_me: {beta_me_global:.4f}")
    print(f"  Reconstructed alpha = beta + 2z - 1 = {beta_me_global + 2*z_m_global - 1:.4f}")

    # Compare to predictions
    alpha_linear = 0.260 * q + 0.827
    alpha_power = 0.69 * q**0.69
    alpha_quad = -0.009 * q**2 + 0.389 * q + 0.386
    print(f"\n  Predictions:")
    print(f"    Linear  (0.260q+0.827):        {alpha_linear:.3f}")
    print(f"    Power-law (0.69*q^0.69):       {alpha_power:.3f}")
    print(f"    Quadratic (-0.009q²+0.389q+0.386): {alpha_quad:.3f}")
    print(f"    Measured (global):             {alpha_global:.3f}")
    print(f"    Residual linear:               {alpha_global - alpha_linear:+.3f}")
    print(f"    Residual power-law:            {alpha_global - alpha_power:+.3f}")
    print(f"    Residual quadratic:            {alpha_global - alpha_quad:+.3f}")

    results['pairwise'] = {
        'z_m': [float(z) for z in pw_z],
        'beta_me': [float(b) for b in pw_b],
        'alpha': [float(a) for a in pw_a],
        'N_mid': [float(n) for n in pw_N],
        'alpha_global': float(alpha_global),
        'z_m_global': float(z_m_global),
        'beta_me_global': float(beta_me_global),
        'alpha_pred_linear': float(alpha_linear),
        'alpha_pred_power': float(alpha_power),
        'alpha_pred_quadratic': float(alpha_quad),
    }

save()
print("\nResults saved.")
