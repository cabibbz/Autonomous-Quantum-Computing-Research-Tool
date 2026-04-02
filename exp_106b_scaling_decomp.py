"""Sprint 106b: Size scaling of chi_F spectral components.

From 106a: chi_F is dominated by ONE state (the (q-1)-fold degenerate multiplet).
chi_F ≈ |<multiplet|H_field|0>|^2 / (E_multiplet - E_0)^2 / N

Decompose: chi_F ~ N^alpha means:
  |me|^2 / gap_m^2 ~ N^(alpha+1)  [since chi_F = chi_total/N]

Track separately:
  gap_m(N) ~ N^{-z_m}  → gap_m^{-2} ~ N^{2z_m}
  |me|^2(N) ~ N^{beta_me}
  => alpha + 1 = beta_me + 2*z_m  => alpha = beta_me + 2*z_m - 1

This tells us whether super-scaling comes from gap closing, matrix element growth, or both.

Extended sizes using GPU for q=5.
"""
import numpy as np
import json, time, os
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '106b_scaling_decomp',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_106b_scaling_decomp.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts(n, q):
    """Build coupling and field parts separately."""
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

# For each q, get the dominant multiplet and track gap_m, |me|^2, chi_F vs N
# Need enough eigenstates to capture the (q-1)-fold multiplet
configs = {
    2: {'sizes': [6, 8, 10, 12, 14], 'k_need': 4},   # level 2, need k≥3
    3: {'sizes': [6, 8, 10], 'k_need': 6},             # level 3, need k≥4
    5: {'sizes': [6, 8, 9], 'k_need': 8},              # level 5, need k≥6 (GPU for n=9)
    7: {'sizes': [6, 7], 'k_need': 10},                # level 7, need k≥8
}

print("=" * 70)
print("Sprint 106b: Size scaling of chi_F spectral components")
print("=" * 70)
print("Decompose: chi_F ~ N^alpha via gap_m ~ N^{-z_m}, |me|^2 ~ N^{beta_me}")
print("=> alpha = beta_me + 2*z_m - 1\n")

for q in [2, 3, 5, 7]:
    cfg = configs[q]
    g_c = 1.0 / q
    k_need = cfg['k_need']

    gap_ms = []
    me_sqs = []
    chi_Fs = []
    Ns = []

    results['data'][f'q{q}'] = {'q': q, 'g_c': g_c, 'sizes': {}}

    for n in cfg['sizes']:
        dim = q**n
        if dim > 15_000_000:
            print(f"  q={q}, n={n}: dim={dim} too large, skipping")
            continue

        print(f"q={q}, n={n} (dim={dim})...", flush=True)
        t0 = time.time()

        H_coup, H_field = build_sq_potts_parts(n, q)
        H = H_coup + g_c * H_field

        # Get enough states to capture the multiplet
        k_use = min(k_need, dim - 2)
        evals, evecs = eigsh(H, k=k_use, which='SA')
        order = np.argsort(evals)
        evals = evals[order]
        evecs = evecs[:, order]

        E0 = evals[0]
        psi0 = evecs[:, 0]
        H_field_psi0 = H_field.dot(psi0)

        # Find the dominant state (largest chi contribution)
        best_chi = 0
        best_gap = 0
        best_me_sq = 0
        best_level = -1
        # Also track the spectral gap
        gap_spectral = evals[1] - E0

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
        dt = time.time() - t0

        gap_ms.append(best_gap)
        me_sqs.append(best_me_sq)
        chi_Fs.append(chi_F)
        Ns.append(n)

        results['data'][f'q{q}']['sizes'][str(n)] = {
            'dim': dim, 'chi_F': float(chi_F),
            'gap_spectral': float(gap_spectral),
            'gap_multiplet': float(best_gap),
            'me_sq': float(best_me_sq),
            'dominant_level': best_level,
            'dominant_fraction': float(best_chi / total_chi) if total_chi > 0 else 0,
            'time_s': dt,
        }

        record(sprint=106, model='S_q_Potts', q=q, n=n,
               quantity='gap_multiplet', value=best_gap, method='spectral_decomp')
        record(sprint=106, model='S_q_Potts', q=q, n=n,
               quantity='me_sq_multiplet', value=best_me_sq, method='spectral_decomp')

        print(f"  gap_spec={gap_spectral:.6f}  gap_mult={best_gap:.6f}  "
              f"|me|²={best_me_sq:.4f}  chi_F={chi_F:.4f}  "
              f"frac={best_chi/total_chi:.4f}  ({dt:.1f}s)")

        save()

        if dt > 100:
            print("  Time limit approaching, skipping remaining sizes")
            break

    # Power-law fits
    if len(Ns) >= 2:
        log_N = np.log(np.array(Ns, dtype=float))
        log_gap = np.log(np.array(gap_ms))
        log_me = np.log(np.array(me_sqs))
        log_chi = np.log(np.array(chi_Fs))

        # Fit: gap_m ~ N^{-z_m}
        z_m = -np.polyfit(log_N, log_gap, 1)[0]
        # Fit: |me|^2 ~ N^{beta_me}
        beta_me = np.polyfit(log_N, log_me, 1)[0]
        # Fit: chi_F ~ N^{alpha}
        alpha = np.polyfit(log_N, log_chi, 1)[0]
        # Predicted alpha
        alpha_pred = beta_me + 2 * z_m - 1

        # Pairwise exponents
        print(f"\n  q={q} SCALING FITS ({len(Ns)} sizes):")
        print(f"    z_m (gap_mult ~ N^{{-z_m}}):  {z_m:.4f}")
        print(f"    beta_me (|me|² ~ N^beta):    {beta_me:.4f}")
        print(f"    alpha (chi_F ~ N^alpha):     {alpha:.4f}")
        print(f"    Predicted alpha = beta + 2z - 1 = {alpha_pred:.4f}")
        print(f"    Known alpha (Sprint 103):     {0.315*q + 0.469:.3f}")

        # Pairwise
        print(f"    Pairwise z_m: ", end="")
        for i in range(len(Ns)-1):
            pw_z = -(log_gap[i+1] - log_gap[i]) / (log_N[i+1] - log_N[i])
            print(f"({Ns[i]},{Ns[i+1]})→{pw_z:.3f}  ", end="")
        print()
        print(f"    Pairwise beta: ", end="")
        for i in range(len(Ns)-1):
            pw_b = (log_me[i+1] - log_me[i]) / (log_N[i+1] - log_N[i])
            print(f"({Ns[i]},{Ns[i+1]})→{pw_b:.3f}  ", end="")
        print()
        print(f"    Pairwise alpha: ", end="")
        for i in range(len(Ns)-1):
            pw_a = (log_chi[i+1] - log_chi[i]) / (log_N[i+1] - log_N[i])
            print(f"({Ns[i]},{Ns[i+1]})→{pw_a:.3f}  ", end="")
        print()

        results['data'][f'q{q}']['z_m'] = float(z_m)
        results['data'][f'q{q}']['beta_me'] = float(beta_me)
        results['data'][f'q{q}']['alpha_fit'] = float(alpha)
        results['data'][f'q{q}']['alpha_predicted'] = float(alpha_pred)
        results['data'][f'q{q}']['alpha_known'] = float(0.315*q + 0.469)

    save()

# Final summary
print(f"\n{'=' * 70}")
print("SUMMARY: chi_F = |me|^2 / gap_m^2 / N decomposition")
print(f"{'q':>3} {'z_m':>8} {'beta_me':>8} {'alpha':>8} {'pred':>8} {'known':>8}")
print("-" * 50)
for q in [2, 3, 5, 7]:
    d = results['data'].get(f'q{q}', {})
    if 'z_m' in d:
        print(f"{q:>3} {d['z_m']:>8.4f} {d['beta_me']:>8.4f} "
              f"{d['alpha_fit']:>8.4f} {d['alpha_predicted']:>8.4f} {d['alpha_known']:>8.3f}")

save()
print("\nResults saved.")
