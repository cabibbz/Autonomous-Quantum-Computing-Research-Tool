"""Sprint 109a: Extend q=3 chi_F spectral decomposition to n=12, 14 via GPU.

q=3 exact diag had n=4-10 (Sprint 108b). GPU enables n=12 (dim=531k) and n=14 (dim=4.8M).
These two extra sizes push 1/ln(N) range from 0.434 to 0.379, giving 13% more leverage
for distinguishing power-law from log corrections.
"""
import numpy as np
import json, time, os
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '109a_q3_gpu_extend',
    'sprint': 109,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_109a_q3_gpu_extend.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts(n, q):
    """Build coupling and field parts separately for S_q Potts."""
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

q = 3
g_c = 1.0 / q
# All sizes: existing n=4-10 + new n=12, 14
all_sizes = [4, 5, 6, 7, 8, 9, 10, 12, 14]
k_need = 6  # q=3: spectral gap = 2-fold, dominant singlet at level 3

print("=" * 70)
print(f"Sprint 109a: q={q} chi_F spectral decomposition, n=4-14")
print("=" * 70)

for n in all_sizes:
    dim = q**n
    print(f"\n  n={n} (dim={dim:,})...", end="", flush=True)
    t0 = time.time()

    H_coup, H_field = build_sq_potts_parts(n, q)
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

    results['data'][str(n)] = {
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

    # Record to DB
    record(sprint=109, model='sq_potts', q=q, n=n,
           quantity='chi_F', value=float(chi_F), method='spectral_exact_gpu')
    record(sprint=109, model='sq_potts', q=q, n=n,
           quantity='gap_multiplet', value=float(best_gap), method='spectral_exact_gpu')
    record(sprint=109, model='sq_potts', q=q, n=n,
           quantity='me_sq', value=best_me_sq, method='spectral_exact_gpu')

# Pairwise analysis
print(f"\n{'=' * 70}")
print("PAIRWISE EXPONENTS")
print("=" * 70)

Ns = []
gaps = []
mes = []
chis = []
for n in all_sizes:
    d = results['data'].get(str(n))
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
from numpy.polynomial.polynomial import polyfit as npfit
p_a = np.polyfit(log_N, log_chi, 1)
alpha_global = p_a[0]
print(f"\n  Global alpha (all sizes): {alpha_global:.4f}")
p_a_last5 = np.polyfit(log_N[-5:], log_chi[-5:], 1)
alpha_last5 = p_a_last5[0]
print(f"  Global alpha (last 5):    {alpha_last5:.4f}")

results['pairwise'] = {
    'z_m': [float(z) for z in pw_z],
    'beta_me': [float(b) for b in pw_b],
    'alpha': [float(a) for a in pw_a],
    'N_mid': [float(n) for n in pw_N],
    'alpha_global': float(alpha_global),
    'alpha_last5': float(alpha_last5),
}
save()
print("\nResults saved.")
