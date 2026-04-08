"""Sprint 121b: S_q Potts q=4 chi_F spectral decomposition at n=4-11.

Sprint 118 measured S_q q=4 with fewer sizes. This extends to n=4-11 (8 sizes)
to match Sprint 120b's hybrid data and compare pairwise alpha drift directly.

S_q q=4 is the marginal case: alpha should drift toward 2.0 from below.
Hybrid q=4 drifts downward (1.64->1.49). If S_q drifts upward, this confirms
different universality classes at q=4.

S_q g_c = 1/4 = 0.25 (Kramers-Wannier self-duality, exact).
"""
import numpy as np
import json, time, os
from gpu_utils import eigsh
from hamiltonian_utils import build_sq_potts_parts
from fss_utils import fit_power_law, pairwise_exponents, fit_spectral_exponents
from db_utils import record

results = {
    'experiment': '121b_sq_q4_extended',
    'sprint': 121,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'model': 'sq_potts',
    'q': 4,
    'g_c': 0.25,
    'data': [],
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results', 'sprint_121b_sq_q4_extended.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

q = 4
g_c = 0.25  # exact self-dual point
sizes = [4, 5, 6, 7, 8, 9, 10, 11]
k_need = 12

print("=" * 70)
print(f"S_q POTTS q={q}, g_c={g_c}, sizes={sizes}")
print("=" * 70)

chi_F_list, gap_m_list, me_sq_list, n_list = [], [], [], []

for n in sizes:
    dim = q**n
    if dim > 12_000_000:
        print(f"\n  n={n} (dim={dim:,}) -- SKIP (too large)")
        continue
    print(f"\n  n={n} (dim={dim:,})...", end="", flush=True)
    t0 = time.time()

    H_coup, H_field = build_sq_potts_parts(n, q)
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

    entry = {
        'n': n, 'dim': dim,
        'gap_multiplet': float(best_gap),
        'gap_spectral': float(spectral_gap) if spectral_gap else None,
        'dominant_level': best_level,
        'dominant_degeneracy': degen_count,
        'me_sq': best_me_sq,
        'chi_F': float(chi_F),
        'dominant_fraction': float(dominant_frac),
        'E0': float(E0),
        'time_s': round(dt, 1),
    }
    results['data'].append(entry)
    save()

    record(sprint=121, model='sq_potts', q=q, n=n,
           quantity='chi_F', value=float(chi_F), method='spectral_exact_gpu')
    record(sprint=121, model='sq_potts', q=q, n=n,
           quantity='gap_multiplet', value=float(best_gap), method='spectral_exact_gpu')
    record(sprint=121, model='sq_potts', q=q, n=n,
           quantity='me_sq', value=best_me_sq, method='spectral_exact_gpu')

    chi_F_list.append(chi_F)
    gap_m_list.append(best_gap)
    me_sq_list.append(best_me_sq)
    n_list.append(n)

# --- Scaling analysis ---
print(f"\n{'=' * 70}")
print("SCALING ANALYSIS: S_q q=4")
print("=" * 70)

sizes_arr = np.array(n_list, dtype=float)
chi_arr = np.array(chi_F_list)
gap_arr = np.array(gap_m_list)
me_arr = np.array(me_sq_list)

# Pairwise exponents
log_N = np.log(sizes_arr)
log_chi = np.log(chi_arr)
log_gap = np.log(gap_arr)
log_me = np.log(me_arr)

print(f"\n  {'Pair':>8s}  {'alpha':>8s}  {'z_m':>8s}  {'beta_me':>8s}  {'recon':>8s}")
pw_alphas, pw_zms = [], []
for i in range(len(sizes_arr)-1):
    dln = log_N[i+1] - log_N[i]
    pw_a = (log_chi[i+1] - log_chi[i]) / dln
    pw_z = -(log_gap[i+1] - log_gap[i]) / dln
    pw_b = (log_me[i+1] - log_me[i]) / dln
    recon = pw_b + 2*pw_z - 1
    pw_alphas.append(pw_a)
    pw_zms.append(pw_z)
    print(f"  ({int(sizes_arr[i]):2d},{int(sizes_arr[i+1]):2d})  {pw_a:8.4f}  {pw_z:8.4f}  {pw_b:8.4f}  {recon:8.4f}")

# Global fit
decomp = fit_spectral_exponents(sizes_arr, chi_arr, gap_arr, me_arr)
print(f"\n  Global fit: alpha={decomp['alpha']:.4f}±{decomp['alpha_err']:.4f}, "
      f"z_m={decomp['z_m']:.4f}±{decomp['z_m_err']:.4f}, "
      f"beta_me={decomp['beta_me']:.4f}±{decomp['beta_me_err']:.4f}")
print(f"  Reconstruction: {decomp['beta_me']+2*decomp['z_m']-1:.4f} vs alpha={decomp['alpha']:.4f}")

results['global_alpha'] = float(decomp['alpha'])
results['global_alpha_err'] = float(decomp['alpha_err'])
results['global_z_m'] = float(decomp['z_m'])
results['global_z_m_err'] = float(decomp['z_m_err'])
results['global_beta_me'] = float(decomp['beta_me'])
results['pairwise_alphas'] = [float(x) for x in pw_alphas]
results['pairwise_z_ms'] = [float(x) for x in pw_zms]

# Alpha drift comparison with hybrid
print(f"\n{'=' * 70}")
print("COMPARISON: alpha drift S_q vs Hybrid at q=4")
print("=" * 70)
# Hybrid pairwise alphas from Sprint 120b
hybrid_pw = [1.637, 1.577, 1.543, 1.522, 1.507, 1.497, 1.489]
hybrid_pairs = ['(4,5)', '(5,6)', '(6,7)', '(7,8)', '(8,9)', '(9,10)', '(10,11)']
print(f"  {'Pair':>8s}  {'S_q alpha':>10s}  {'Hybrid alpha':>13s}  {'Direction':>10s}")
for i, pair in enumerate(hybrid_pairs):
    sq_a = pw_alphas[i] if i < len(pw_alphas) else None
    hy_a = hybrid_pw[i]
    if sq_a is not None:
        sq_str = f"{sq_a:.4f}"
        direction = "UP" if i > 0 and pw_alphas[i] > pw_alphas[i-1] else "DOWN" if i > 0 else "-"
    else:
        sq_str = "N/A"
        direction = "-"
    print(f"  {pair:>8s}  {sq_str:>10s}  {hy_a:13.4f}  {direction:>10s}")

# Check if S_q alpha drifts upward (toward 2.0 as expected)
if len(pw_alphas) >= 3:
    first_half = np.mean(pw_alphas[:len(pw_alphas)//2])
    second_half = np.mean(pw_alphas[len(pw_alphas)//2:])
    drift = "UPWARD" if second_half > first_half else "DOWNWARD"
    print(f"\n  S_q alpha drift: {drift} ({first_half:.4f} -> {second_half:.4f})")
    print(f"  Hybrid alpha drift: DOWNWARD (1.637 -> 1.489)")
    results['sq_drift'] = drift
    results['sq_first_half_alpha'] = float(first_half)
    results['sq_second_half_alpha'] = float(second_half)

save()
print("\nResults saved.")
