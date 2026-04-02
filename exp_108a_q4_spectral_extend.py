"""Sprint 108a: q=4 spectral decomposition extended to n=4-10 (GPU for n>=8).

q=4 is the BKT boundary in S_q Potts. Sprint 106c measured alpha=1.69 with
sizes n=4-8. Extend to n=9 (dim=262k) and n=10 (dim=1M) using GPU.

Track: gap_m, |me|^2, chi_F, dominant fraction for each size.
Decompose: alpha = beta_me + 2*z_m - 1.
"""
import numpy as np
import json, time, os
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '108a_q4_spectral_extend',
    'sprint': 108,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': 4,
    'g_c': 0.25,
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_108a_q4_spectral_extend.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts(n, q):
    """Build coupling and field parts separately: H(g) = H_coup + g * H_field.
    S_q Potts: field = sum of all cyclic permutations X^k for k=1..q-1."""
    dim = q**n
    H_coup = lil_matrix((dim, dim), dtype=float)
    H_field = lil_matrix((dim, dim), dtype=float)

    for idx in range(dim):
        c = []
        tmp = idx
        for _ in range(n):
            c.append(tmp % q)
            tmp //= q

        # Diagonal: Potts coupling -delta(s_i, s_j)
        diag = 0.0
        for site in range(n):
            nxt = (site + 1) % n
            if c[site] == c[nxt]:
                diag -= 1.0
        H_coup[idx, idx] = diag

        # Off-diagonal: S_q transverse field (all cyclic shifts)
        for site in range(n):
            for k in range(1, q):
                c2 = c.copy()
                c2[site] = (c[site] + k) % q
                idx2 = 0
                for i in range(n - 1, -1, -1):
                    idx2 = idx2 * q + c2[i]
                H_field[idx, idx2] += -1.0

    return csr_matrix(H_coup), csr_matrix(H_field)

q = 4
g_c = 1.0 / q
# n=4 (256), n=5 (1024), n=6 (4096), n=7 (16384), n=8 (65536), n=9 (262144), n=10 (1048576)
sizes = [4, 5, 6, 7, 8, 9, 10]
k_need = 8  # Need at least q+1=5 states to capture (q-1)-fold multiplet

print("=" * 70)
print("Sprint 108a: q=4 spectral decomposition, n=4-10")
print(f"g_c = {g_c}, k_states = {k_need}")
print("=" * 70)

gap_ms = []
me_sqs = []
chi_Fs = []
Ns = []
total_time = 0

for n in sizes:
    dim = q**n
    print(f"\nn={n} (dim={dim:,})...", flush=True)
    t0 = time.time()

    # Build Hamiltonian parts
    H_coup, H_field = build_sq_potts_parts(n, q)
    H = H_coup + g_c * H_field
    t_build = time.time() - t0
    print(f"  Built in {t_build:.1f}s", flush=True)

    # Get k lowest eigenstates
    k_use = min(k_need, dim - 2)
    evals, evecs = eigsh(H, k=k_use, which='SA')
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]

    E0 = evals[0]
    psi0 = evecs[:, 0]
    H_field_psi0 = H_field.dot(psi0)

    t_diag = time.time() - t0
    print(f"  Diagonalized in {t_diag:.1f}s", flush=True)

    # Compute spectral decomposition of chi_F
    contributions = []
    total_chi = 0
    best_chi = 0
    best_gap = 0
    best_me_sq = 0
    best_level = -1
    gap_spectral = evals[1] - E0

    for i in range(1, k_use):
        gap = evals[i] - E0
        if gap < 1e-12:
            continue
        me = np.dot(evecs[:, i], H_field_psi0)
        chi_n = me**2 / gap**2
        total_chi += chi_n
        contributions.append({
            'level': i,
            'gap': float(gap),
            'matrix_element': float(me),
            'me_sq': float(me**2),
            'chi_contribution': float(chi_n),
        })
        if chi_n > best_chi:
            best_chi = chi_n
            best_gap = gap
            best_me_sq = me**2
            best_level = i

    chi_F = total_chi / n
    dominant_frac = best_chi / total_chi if total_chi > 0 else 0

    # Sort contributions by magnitude
    contributions.sort(key=lambda x: -x['chi_contribution'])

    # Cumulative analysis
    cumulative = 0.0
    for c in contributions:
        cumulative += c['chi_contribution']
        c['cumulative_fraction'] = float(cumulative / total_chi) if total_chi > 0 else 0

    dt = time.time() - t0
    total_time += dt

    gap_ms.append(best_gap)
    me_sqs.append(best_me_sq)
    chi_Fs.append(chi_F)
    Ns.append(n)

    results['data'][str(n)] = {
        'n': n, 'dim': dim,
        'chi_F': float(chi_F),
        'chi_total_Nchi': float(total_chi),
        'gap_spectral': float(gap_spectral),
        'gap_multiplet': float(best_gap),
        'me_sq': float(best_me_sq),
        'dominant_level': best_level,
        'dominant_fraction': float(dominant_frac),
        'E0': float(E0),
        'eigenvalues': [float(e) for e in evals],
        'contributions_top5': contributions[:5],
        'time_s': dt,
    }

    record(sprint=108, model='S_q_Potts', q=4, n=n,
           quantity='chi_F', value=chi_F, method='spectral_decomp')
    record(sprint=108, model='S_q_Potts', q=4, n=n,
           quantity='gap_multiplet', value=best_gap, method='spectral_decomp')
    record(sprint=108, model='S_q_Potts', q=4, n=n,
           quantity='me_sq_multiplet', value=best_me_sq, method='spectral_decomp')
    record(sprint=108, model='S_q_Potts', q=4, n=n,
           quantity='chi_F_dominant_frac', value=dominant_frac, method='spectral_decomp')

    print(f"  gap_spec={gap_spectral:.6f}  gap_mult={best_gap:.6f}")
    print(f"  |me|²={best_me_sq:.6f}  chi_F={chi_F:.4f}  dom_frac={dominant_frac:.4f}")
    print(f"  Dominant level: {best_level}  Time: {dt:.1f}s")

    save()

    if total_time > 240:
        print("  WARNING: approaching 300s total limit, stopping")
        break

# Power-law fits
print(f"\n{'=' * 70}")
print("SCALING ANALYSIS")
print("=" * 70)

log_N = np.log(np.array(Ns, dtype=float))
log_gap = np.log(np.array(gap_ms))
log_me = np.log(np.array(me_sqs))
log_chi = np.log(np.array(chi_Fs))

# Global fits
z_m = -np.polyfit(log_N, log_gap, 1)[0]
beta_me = np.polyfit(log_N, log_me, 1)[0]
alpha = np.polyfit(log_N, log_chi, 1)[0]
alpha_pred = beta_me + 2 * z_m - 1

print(f"\nGlobal fits ({len(Ns)} sizes, n={Ns[0]}-{Ns[-1]}):")
print(f"  z_m (gap_mult ~ N^{{-z_m}}):   {z_m:.4f}")
print(f"  beta_me (|me|² ~ N^beta):    {beta_me:.4f}")
print(f"  alpha (chi_F ~ N^alpha):     {alpha:.4f}")
print(f"  Predicted alpha = beta + 2z - 1 = {alpha_pred:.4f}")
print(f"  Linear formula: alpha(4) = 0.315*4 + 0.469 = {0.315*4 + 0.469:.3f}")
print(f"  Sprint 106c: alpha = 1.69")

# Pairwise exponents
print(f"\nPairwise exponents:")
print(f"  {'(N1,N2)':>10} {'z_m':>8} {'beta_me':>8} {'alpha':>8} {'alpha_pred':>10}")
print(f"  {'-'*46}")
for i in range(len(Ns)-1):
    dln = log_N[i+1] - log_N[i]
    pw_z = -(log_gap[i+1] - log_gap[i]) / dln
    pw_b = (log_me[i+1] - log_me[i]) / dln
    pw_a = (log_chi[i+1] - log_chi[i]) / dln
    pw_pred = pw_b + 2*pw_z - 1
    print(f"  ({Ns[i]},{Ns[i+1]}){'':<4} {pw_z:>8.4f} {pw_b:>8.4f} {pw_a:>8.4f} {pw_pred:>10.4f}")

results['scaling'] = {
    'z_m': float(z_m),
    'beta_me': float(beta_me),
    'alpha_fit': float(alpha),
    'alpha_predicted': float(alpha_pred),
    'alpha_linear_formula': float(0.315*4 + 0.469),
    'sizes_used': Ns,
}

save()
print(f"\nTotal time: {total_time:.1f}s")
print("Results saved.")
