"""Sprint 106c: Mechanism analysis of chi_F super-scaling.

From 106b: alpha = beta_me + 2*z_m - 1, with both z_m and beta_me increasing with q.
Key questions:
1. What are z_m and beta_me for q=4 (the BKT crossover)?
2. How does z_m relate to the spectral gap exponent z_spec?
3. Can we extract z_m(q) and beta_me(q) as simple functions?
4. Is there a connection to the entanglement spectrum parameters?

Also: verify that gap_mult/gap_spec ratio is q-dependent (multiplet lives higher in spectrum
for walking, meaning it closes relatively slower).
"""
import numpy as np
import json, time, os
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '106c_mechanism',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_106c_mechanism.json')
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

# q=4 data + extended q=2 for comparison + q=6
configs = {
    4: {'sizes': [4, 5, 6, 7, 8], 'k_need': 8},
    6: {'sizes': [5, 6, 7], 'k_need': 10},
}

print("=" * 70)
print("Sprint 106c: Mechanism analysis — q=4 crossover and q=6 broken walking")
print("=" * 70)

for q in sorted(configs.keys()):
    cfg = configs[q]
    g_c = 1.0 / q
    k_need = cfg['k_need']

    gap_specs = []
    gap_mults = []
    me_sqs = []
    chi_Fs = []
    Ns = []

    results['data'][f'q{q}'] = {'q': q, 'g_c': g_c, 'sizes': {}}

    for n in cfg['sizes']:
        dim = q**n
        if dim > 5_000_000:
            print(f"  q={q}, n={n}: dim={dim}, attempting with GPU...", flush=True)
        else:
            print(f"q={q}, n={n} (dim={dim})...", flush=True)

        if dim > 15_000_000:
            print(f"  Skipping (too large)")
            continue

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
        dt = time.time() - t0

        gap_specs.append(gap_spectral)
        gap_mults.append(best_gap)
        me_sqs.append(best_me_sq)
        chi_Fs.append(chi_F)
        Ns.append(n)

        # Gap ratio
        gap_ratio = best_gap / gap_spectral

        results['data'][f'q{q}']['sizes'][str(n)] = {
            'dim': dim, 'chi_F': float(chi_F),
            'gap_spectral': float(gap_spectral),
            'gap_multiplet': float(best_gap),
            'gap_ratio': float(gap_ratio),
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
              f"ratio={gap_ratio:.3f}  |me|²={best_me_sq:.4f}  chi_F={chi_F:.4f}  ({dt:.1f}s)")

        save()

        if dt > 100:
            print("  Time limit, stopping for this q")
            break

    # Fits
    if len(Ns) >= 2:
        log_N = np.log(np.array(Ns, dtype=float))
        log_gap_m = np.log(np.array(gap_mults))
        log_gap_s = np.log(np.array(gap_specs))
        log_me = np.log(np.array(me_sqs))
        log_chi = np.log(np.array(chi_Fs))

        z_m = -np.polyfit(log_N, log_gap_m, 1)[0]
        z_s = -np.polyfit(log_N, log_gap_s, 1)[0]
        beta_me = np.polyfit(log_N, log_me, 1)[0]
        alpha = np.polyfit(log_N, log_chi, 1)[0]
        alpha_pred = beta_me + 2 * z_m - 1

        results['data'][f'q{q}']['z_m'] = float(z_m)
        results['data'][f'q{q}']['z_spec'] = float(z_s)
        results['data'][f'q{q}']['beta_me'] = float(beta_me)
        results['data'][f'q{q}']['alpha_fit'] = float(alpha)
        results['data'][f'q{q}']['alpha_predicted'] = float(alpha_pred)

        print(f"\n  q={q} SCALING:")
        print(f"    z_spec = {z_s:.4f}  (spectral gap)")
        print(f"    z_m    = {z_m:.4f}  (multiplet gap)")
        print(f"    beta_me= {beta_me:.4f}")
        print(f"    alpha  = {alpha:.4f}  (predicted: {alpha_pred:.4f})")

    save()

# Now compile all results from 106b + 106c
print(f"\n{'=' * 70}")
print("COMPLETE TABLE: chi_F mechanism across q=2-7")
print(f"{'q':>3} {'z_spec':>8} {'z_m':>8} {'z_m/z_s':>8} {'beta_me':>8} "
      f"{'alpha':>8} {'alpha_kn':>8} {'mechanism':>20}")
print("-" * 85)

# Include 106b results
all_q = {
    2: {'z_m': 0.9891, 'z_s': 1.0, 'beta_me': 0.0663, 'alpha': 1.0444},
    3: {'z_m': 1.0306, 'z_s': 1.0, 'beta_me': 0.3810, 'alpha': 1.4423},
    5: {'z_m': 1.1564, 'z_s': 1.0, 'beta_me': 0.7688, 'alpha': 2.0817},
    7: {'z_m': 1.3103, 'z_s': 1.0, 'beta_me': 1.0099, 'alpha': 2.6305},
}

for q in [4, 6]:
    d = results['data'].get(f'q{q}', {})
    if 'z_m' in d:
        all_q[q] = {
            'z_m': d['z_m'], 'z_s': d['z_spec'],
            'beta_me': d['beta_me'], 'alpha': d['alpha_fit'],
        }

for q in sorted(all_q.keys()):
    d = all_q[q]
    alpha_known = 0.315 * q + 0.469
    z_ratio = d['z_m'] / d.get('z_s', d['z_m']) if d.get('z_s', 0) != 0 else 0

    if d['beta_me'] < 0.1:
        mech = "gap only"
    elif d['z_m'] > 1.1:
        mech = "BOTH (walking)"
    else:
        mech = "mixed"

    print(f"{q:>3} {d.get('z_s', 0):>8.4f} {d['z_m']:>8.4f} {z_ratio:>8.3f} "
          f"{d['beta_me']:>8.4f} {d['alpha']:>8.4f} {alpha_known:>8.3f} {mech:>20}")

# Fit z_m(q) and beta_me(q)
qs = sorted(all_q.keys())
z_ms = [all_q[q]['z_m'] for q in qs]
betas = [all_q[q]['beta_me'] for q in qs]

# Linear fits
z_fit = np.polyfit(qs, z_ms, 1)
b_fit = np.polyfit(qs, betas, 1)

print(f"\nLinear fits:")
print(f"  z_m(q)    = {z_fit[0]:.4f}*q + {z_fit[1]:.4f}")
print(f"  beta_me(q)= {b_fit[0]:.4f}*q + {b_fit[1]:.4f}")
print(f"  alpha(q)  = 2*z_m + beta - 1 = {2*z_fit[0]+b_fit[0]:.4f}*q + {2*z_fit[1]+b_fit[1]-1:.4f}")
print(f"  Known:      alpha(q) = 0.315*q + 0.469")

results['z_m_fit'] = {'slope': float(z_fit[0]), 'intercept': float(z_fit[1])}
results['beta_me_fit'] = {'slope': float(b_fit[0]), 'intercept': float(b_fit[1])}

save()
print("\nResults saved.")
