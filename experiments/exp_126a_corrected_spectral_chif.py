"""Sprint 126a: Corrected spectral chi_F with factor-2 fix + systematic exponent extraction.

Factor-2 correction: chi_F = (2/N) * sum_{n>0} |<n|V|0>|^2 / (E_n - E_0)^2
Sprint 125 found the spectral sum was missing the factor 2 from the second derivative
of fidelity. Exponents should be unaffected (constant prefactor), but we verify explicitly.

Tests:
- S_q Potts at q=3,4,5,7 (walking regime)
- Hybrid Potts-clock at q=3,4,5,7 (continuous regime)
- All feasible sizes with periodic BC
- Compare corrected spectral vs exact chi_F
- Extract alpha, z_m, beta_me with error bars via fss_utils
"""
import numpy as np
import json, time, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from gpu_utils import eigsh
from hamiltonian_utils import build_sq_potts_parts, build_hybrid_parts
from fss_utils import fit_power_law, pairwise_exponents, fit_spectral_exponents
from db_utils import record

results = {
    'experiment': '126a_corrected_spectral_chif',
    'sprint': 126,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'configs': {},
    'exponents': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results',
                           'sprint_126a_corrected_spectral_chif.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)


def chi_F_exact(H_coup, H_field, g_c, n, dg=1e-4):
    """Exact chi_F via central finite-difference of ground state overlap."""
    Hp = H_coup + (g_c + dg) * H_field
    Hm = H_coup + (g_c - dg) * H_field
    H0 = H_coup + g_c * H_field

    _, vecs0 = eigsh(H0, k=1, which='SA')
    _, vecsp = eigsh(Hp, k=1, which='SA')
    _, vecsm = eigsh(Hm, k=1, which='SA')

    psi0 = vecs0[:, 0]
    psip = vecsp[:, 0]
    psim = vecsm[:, 0]

    if np.dot(psi0, psip) < 0: psip = -psip
    if np.dot(psi0, psim) < 0: psim = -psim

    ov_p = np.dot(psi0, psip)**2
    ov_m = np.dot(psi0, psim)**2
    return (2.0 - ov_p - ov_m) / (dg**2 * n)


def chi_F_spectral_corrected(H_coup, H_field, g_c, n, k=20):
    """Corrected spectral chi_F with factor 2.

    chi_F = (2/N) * sum_{n>0} |<n|V|0>|^2 / (E_n - E_0)^2

    Returns dict with chi_F, dominant fraction, gap, dominant |me|^2, n_contributing.
    """
    H = H_coup + g_c * H_field
    dim = H.shape[0]
    k_use = min(k, dim - 2)

    evals, evecs = eigsh(H, k=k_use, which='SA')
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]

    E0 = evals[0]
    psi0 = evecs[:, 0]
    V_psi0 = H_field.dot(psi0)

    contributions = []
    gaps_all = []
    me_sq_all = []
    for i in range(1, k_use):
        gap = evals[i] - E0
        if gap < 1e-12:
            continue
        me = np.dot(evecs[:, i], V_psi0)
        me_sq = me**2
        chi_contrib = me_sq / gap**2
        if chi_contrib > 1e-20:
            contributions.append(chi_contrib)
            gaps_all.append(gap)
            me_sq_all.append(me_sq)

    if not contributions:
        return None

    total = sum(contributions)
    # CORRECTED: factor of 2
    chi_F_val = 2.0 * total / n

    # Dominant state info
    idx_dom = np.argmax(contributions)
    frac = contributions[idx_dom] / total if total > 0 else 0

    return {
        'chi_F': float(chi_F_val),
        'frac_dominant': float(frac),
        'gap_dominant': float(gaps_all[idx_dom]),
        'me_sq_dominant': float(me_sq_all[idx_dom]),
        'n_contributing': len(contributions),
        'total_lehmann_sum': float(total),
    }


# ============================================================
# Test configurations
# ============================================================
configs = [
    # (model, q, g_c, sizes)
    # Keep sizes where 4 eigsh calls per size fit in ~60s budget
    ('sq', 3, 1.0/3, [4, 6, 8, 10, 12]),
    ('sq', 4, 0.25, [4, 6, 8, 10]),
    ('sq', 5, 0.20, [4, 6, 8]),
    ('sq', 7, 1.0/7, [4, 6, 7]),
    ('hybrid', 3, 0.333, [4, 6, 8, 10, 12]),
    ('hybrid', 4, 0.393, [4, 6, 8, 10]),
    ('hybrid', 5, 0.438, [4, 6, 8]),
    ('hybrid', 7, 0.535, [4, 6, 7]),
]

print("=" * 85)
print("Sprint 126a: CORRECTED spectral chi_F (factor-2 fix) + exponent extraction")
print("=" * 85)

for model, q, g_c, sizes in configs:
    label = f"{model}_q{q}"
    print(f"\n{'='*70}")
    print(f"  {model.upper()} q={q}, g_c={g_c:.4f}")
    print(f"{'='*70}")

    size_data = []
    chi_F_arr, gap_arr, me_sq_arr, sizes_used = [], [], [], []

    for n in sizes:
        dim = q**n
        if dim > 12_000_000:
            print(f"  n={n}: dim={dim:,} > 12M, skipping")
            continue

        t0 = time.time()
        if model == 'sq':
            H_coup, H_field = build_sq_potts_parts(n, q)
        else:
            H_coup, H_field = build_hybrid_parts(n, q)

        # Exact chi_F
        exact = chi_F_exact(H_coup, H_field, g_c, n)

        # Corrected spectral chi_F
        spec = chi_F_spectral_corrected(H_coup, H_field, g_c, n, k=20)
        dt = time.time() - t0

        if spec is None:
            print(f"  n={n}: spectral failed")
            continue

        ratio = spec['chi_F'] / exact if exact > 0 else 0
        print(f"  n={n:2d} dim={dim:>10,}  exact={exact:10.4f}  spec={spec['chi_F']:10.4f}  "
              f"ratio={ratio:.4f}  frac={spec['frac_dominant']:.4f}  "
              f"states={spec['n_contributing']}  ({dt:.1f}s)")

        entry = {
            'n': n, 'dim': dim, 'chi_F_exact': float(exact),
            'chi_F_spectral': float(spec['chi_F']),
            'ratio': float(ratio),
            'frac_dominant': float(spec['frac_dominant']),
            'gap_dominant': float(spec['gap_dominant']),
            'me_sq_dominant': float(spec['me_sq_dominant']),
            'n_contributing': spec['n_contributing'],
            'time_s': round(dt, 1),
        }
        size_data.append(entry)

        # For exponent extraction: use corrected spectral values
        sizes_used.append(n)
        chi_F_arr.append(spec['chi_F'])
        gap_arr.append(spec['gap_dominant'])
        me_sq_arr.append(spec['me_sq_dominant'])

        # Record to DB
        record(sprint=126, model=model, q=q, n=n,
               quantity='chi_F_corrected', value=float(spec['chi_F']),
               method='spectral_corrected')
        record(sprint=126, model=model, q=q, n=n,
               quantity='chi_F_exact', value=float(exact),
               method='finite_diff')
        record(sprint=126, model=model, q=q, n=n,
               quantity='frac_dominant', value=float(spec['frac_dominant']),
               method='spectral_corrected')

    results['configs'][label] = size_data
    save()

    # Extract exponents if we have 3+ sizes
    if len(sizes_used) >= 3:
        sizes_np = np.array(sizes_used, dtype=float)
        chi_np = np.array(chi_F_arr)
        gap_np = np.array(gap_arr)
        me_np = np.array(me_sq_arr)

        decomp = fit_spectral_exponents(sizes_np, chi_np, gap_np, me_np)

        print(f"\n  --- Exponents (fss_utils, {len(sizes_used)} sizes) ---")
        print(f"  alpha    = {decomp['alpha']:.4f} ± {decomp['alpha_err']:.4f}  (R²={decomp['r_squared']['chi_F']:.6f})")
        print(f"  z_m      = {decomp['z_m']:.4f} ± {decomp['z_m_err']:.4f}  (R²={decomp['r_squared']['gap']:.6f})")
        print(f"  beta_me  = {decomp['beta_me']:.4f} ± {decomp['beta_me_err']:.4f}  (R²={decomp['r_squared']['me_sq']:.6f})")
        print(f"  recon    = {decomp['reconstruction']:.4f} vs alpha={decomp['alpha']:.4f} (err={decomp['recon_error']:.4f})")

        # Pairwise alpha drift
        pairs = pairwise_exponents(sizes_np, chi_np)
        print(f"  pairwise alpha: ", end="")
        for p in pairs:
            print(f"({p['n1']},{p['n2']})={p['alpha']:.3f}  ", end="")
        print()

        results['exponents'][label] = {
            'alpha': decomp['alpha'], 'alpha_err': decomp['alpha_err'],
            'z_m': decomp['z_m'], 'z_m_err': decomp['z_m_err'],
            'beta_me': decomp['beta_me'], 'beta_me_err': decomp['beta_me_err'],
            'reconstruction': decomp['reconstruction'],
            'recon_error': decomp['recon_error'],
            'r_squared': decomp['r_squared'],
            'pairwise': [{'n1': p['n1'], 'n2': p['n2'], 'alpha': p['alpha']} for p in pairs],
            'sizes_used': sizes_used,
        }

        # Record exponents to DB
        record(sprint=126, model=model, q=q, n=max(sizes_used),
               quantity='alpha_corrected', value=decomp['alpha'],
               error=decomp['alpha_err'], method='spectral_corrected_fss')
        record(sprint=126, model=model, q=q, n=max(sizes_used),
               quantity='z_m_corrected', value=decomp['z_m'],
               error=decomp['z_m_err'], method='spectral_corrected_fss')
        record(sprint=126, model=model, q=q, n=max(sizes_used),
               quantity='beta_me_corrected', value=decomp['beta_me'],
               error=decomp['beta_me_err'], method='spectral_corrected_fss')

        save()

# ============================================================
# Summary table
# ============================================================
print(f"\n{'='*85}")
print("SUMMARY: Corrected exponents with error bars")
print(f"{'='*85}")
print(f"{'Model':>8} {'q':>3} {'alpha':>12} {'z_m':>12} {'beta_me':>12} {'recon':>8} {'R²_chi':>8} {'sizes':>10}")
for label, exp in results['exponents'].items():
    model_q = label.split('_')
    model = model_q[0]
    q_val = model_q[1]
    alpha_str = f"{exp['alpha']:.4f}±{exp['alpha_err']:.4f}"
    zm_str = f"{exp['z_m']:.4f}±{exp['z_m_err']:.4f}"
    bm_str = f"{exp['beta_me']:.4f}±{exp['beta_me_err']:.4f}"
    print(f"{model:>8} {q_val:>3} {alpha_str:>12} {zm_str:>12} {bm_str:>12} "
          f"{exp['recon_error']:>7.4f} {exp['r_squared']['chi_F']:>8.5f} {str(exp['sizes_used']):>10}")

# Verification: spectral/exact ratio
print(f"\n{'='*85}")
print("VERIFICATION: Corrected spectral / exact chi_F ratios (should be ~1.00)")
print(f"{'='*85}")
print(f"{'Config':>12} {'n':>4} {'exact':>10} {'spectral':>10} {'ratio':>8} {'frac':>6}")
for label, data in results['configs'].items():
    for entry in data:
        print(f"{label:>12} {entry['n']:>4} {entry['chi_F_exact']:>10.4f} "
              f"{entry['chi_F_spectral']:>10.4f} {entry['ratio']:>8.4f} {entry['frac_dominant']:>6.3f}")

save()
print(f"\nResults saved.")
