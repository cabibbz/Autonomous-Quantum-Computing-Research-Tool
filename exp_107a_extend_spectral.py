"""Sprint 107a: Extend chi_F spectral decomposition to more sizes.

Hardening Sprint 106 mechanism: alpha = beta_me + 2*z_m - 1.
- q=5: add n=7 and n=10 (GPU, dim=9.8M) to existing n=6,8,9 → 5 sizes
- q=7: add n=8 (GPU, dim=5.7M) to existing n=6,7 → 3 sizes
- Also recompute existing sizes for consistency check

Track: gap_multiplet, |me|^2, chi_F for each, then fit z_m, beta_me, alpha.
"""
import numpy as np
import json, time, os
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '107a_extend_spectral',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_107a_extend_spectral.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts(n, q):
    """Build coupling and field parts of S_q Potts Hamiltonian separately."""
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

# Configurations: include existing sizes for consistency + new sizes
configs = {
    5: {'sizes': [6, 7, 8, 9, 10], 'k_need': 8},   # n=7,10 are NEW
    7: {'sizes': [6, 7, 8], 'k_need': 10},            # n=8 is NEW
}

print("=" * 70)
print("Sprint 107a: Extend chi_F spectral decomposition")
print("=" * 70)

for q in [5, 7]:
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
        print(f"\nq={q}, n={n} (dim={dim:,})...", flush=True)
        t0 = time.time()

        H_coup, H_field = build_sq_potts_parts(n, q)
        H = H_coup + g_c * H_field
        t_build = time.time() - t0
        print(f"  Build: {t_build:.1f}s", flush=True)

        # Get enough states to capture the (q-1)-fold multiplet at level q-1
        k_use = min(k_need, dim - 2)
        evals, evecs = eigsh(H, k=k_use, which='SA')
        order = np.argsort(evals)
        evals = evals[order]
        evecs = evecs[:, order]
        t_eig = time.time() - t0 - t_build
        print(f"  Eigsh: {t_eig:.1f}s (k={k_use})", flush=True)

        E0 = evals[0]
        psi0 = evecs[:, 0]
        H_field_psi0 = H_field.dot(psi0)
        gap_spectral = evals[1] - E0

        # Find dominant state and compute spectral decomposition
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
        dom_frac = best_chi / total_chi if total_chi > 0 else 0
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
            'dominant_fraction': float(dom_frac),
            'time_s': dt,
        }

        record(sprint=107, model='S_q_Potts', q=q, n=n,
               quantity='gap_multiplet', value=best_gap, method='spectral_decomp')
        record(sprint=107, model='S_q_Potts', q=q, n=n,
               quantity='me_sq_multiplet', value=best_me_sq, method='spectral_decomp')
        record(sprint=107, model='S_q_Potts', q=q, n=n,
               quantity='chi_F', value=chi_F, method='spectral_decomp')

        print(f"  gap_spec={gap_spectral:.6f}  gap_mult={best_gap:.6f}  "
              f"|me|²={best_me_sq:.4f}  chi_F={chi_F:.4f}  "
              f"frac={dom_frac:.4f}  level={best_level}  ({dt:.1f}s)")

        save()

    # Power-law fits
    if len(Ns) >= 3:
        log_N = np.log(np.array(Ns, dtype=float))
        log_gap = np.log(np.array(gap_ms))
        log_me = np.log(np.array(me_sqs))
        log_chi = np.log(np.array(chi_Fs))

        z_m = -np.polyfit(log_N, log_gap, 1)[0]
        beta_me = np.polyfit(log_N, log_me, 1)[0]
        alpha = np.polyfit(log_N, log_chi, 1)[0]
        alpha_pred = beta_me + 2 * z_m - 1
        alpha_known = 0.315 * q + 0.469

        print(f"\n  q={q} SCALING FITS ({len(Ns)} sizes, n={Ns}):")
        print(f"    z_m = {z_m:.4f}")
        print(f"    beta_me = {beta_me:.4f}")
        print(f"    alpha(fit) = {alpha:.4f}")
        print(f"    alpha(predicted = beta + 2z - 1) = {alpha_pred:.4f}")
        print(f"    alpha(known from Sprint 103) = {alpha_known:.3f}")

        # Pairwise exponents
        print(f"    Pairwise z_m: ", end="")
        for i in range(len(Ns)-1):
            pw_z = -(log_gap[i+1] - log_gap[i]) / (log_N[i+1] - log_N[i])
            print(f"({Ns[i]},{Ns[i+1]})→{pw_z:.4f}  ", end="")
        print()
        print(f"    Pairwise beta_me: ", end="")
        for i in range(len(Ns)-1):
            pw_b = (log_me[i+1] - log_me[i]) / (log_N[i+1] - log_N[i])
            print(f"({Ns[i]},{Ns[i+1]})→{pw_b:.4f}  ", end="")
        print()
        print(f"    Pairwise alpha: ", end="")
        for i in range(len(Ns)-1):
            pw_a = (log_chi[i+1] - log_chi[i]) / (log_N[i+1] - log_N[i])
            print(f"({Ns[i]},{Ns[i+1]})→{pw_a:.4f}  ", end="")
        print()

        results['data'][f'q{q}']['z_m'] = float(z_m)
        results['data'][f'q{q}']['beta_me'] = float(beta_me)
        results['data'][f'q{q}']['alpha_fit'] = float(alpha)
        results['data'][f'q{q}']['alpha_predicted'] = float(alpha_pred)
        results['data'][f'q{q}']['alpha_known'] = float(alpha_known)

        save()

# Final summary
print(f"\n{'=' * 70}")
print("SUMMARY: Extended spectral decomposition")
print(f"{'q':>3} {'#pts':>5} {'z_m':>8} {'beta_me':>8} {'alpha':>8} {'pred':>8} {'known':>8}")
print("-" * 55)
for q in [5, 7]:
    d = results['data'].get(f'q{q}', {})
    if 'z_m' in d:
        npts = len(d.get('sizes', {}))
        print(f"{q:>3} {npts:>5} {d['z_m']:>8.4f} {d['beta_me']:>8.4f} "
              f"{d['alpha_fit']:>8.4f} {d['alpha_predicted']:>8.4f} {d['alpha_known']:>8.3f}")

save()
print("\nResults saved to results/sprint_107a_extend_spectral.json")
