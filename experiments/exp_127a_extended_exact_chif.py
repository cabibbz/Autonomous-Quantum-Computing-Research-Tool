"""Sprint 127a: Extend exact chi_F to GPU sizes for tighter exponent fits.

Adds new sizes beyond Sprint 126's reach:
- q=3: n=14 (dim=4.8M, GPU)
- q=4: n=9 (262k), n=11 (4.2M, GPU)
- q=5: n=9 (2.0M, GPU), n=10 (9.8M, GPU)
- q=7: n=8 (5.8M, GPU)

Combines with all prior exact chi_F data for complete refit.
Uses exact chi_F (finite-difference) only — the gold standard per Sprint 126.
"""
import numpy as np
import json, time, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from gpu_utils import eigsh
from hamiltonian_utils import build_sq_potts_parts, build_hybrid_parts
from fss_utils import fit_power_law, pairwise_exponents
from db_utils import record, query

results = {
    'experiment': '127a_extended_exact_chif',
    'sprint': 127,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'new_points': {},
    'combined_fits': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results',
                           'sprint_127a_extended_exact_chif.json')
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


# Prior exact chi_F data from Sprint 126 (DB query)
prior_data = {
    # (model, q): {n: chi_F}
    ('sq', 3): {4: 4.086420, 6: 7.717202, 8: 11.868345, 10: 16.456744, 12: 21.427516},
    ('sq', 4): {4: 12.839321, 6: 27.138030, 8: 45.527791, 10: 67.784201},
    ('sq', 5): {4: 30.943214, 6: 72.624695, 8: 132.424150},
    ('sq', 7): {4: 117.243073, 6: 333.985930, 7: 501.054584},
    ('hybrid', 3): {4: 4.094578, 6: 7.732527, 8: 11.891721, 10: 16.488781, 12: 21.468582},
    ('hybrid', 4): {4: 4.312760, 6: 8.354260, 8: 13.064844, 10: 18.341038},
    ('hybrid', 5): {4: 4.003316, 6: 7.390136, 8: 10.988427},
    ('hybrid', 7): {4: 2.775864, 6: 4.330642, 7: 4.943292},
}

# New sizes to compute this sprint
# Each entry: (model, q, g_c, [new_sizes])
new_configs = [
    ('sq', 3, 1.0/3, [14]),           # 3^14 = 4.8M
    ('sq', 4, 0.25, [9, 11]),         # 4^9=262k, 4^11=4.2M
    ('sq', 5, 0.20, [9, 10]),         # 5^9=2.0M, 5^10=9.8M
    ('sq', 7, 1.0/7, [8]),            # 7^8=5.8M
    ('hybrid', 3, 0.333, [14]),
    ('hybrid', 4, 0.393, [9, 11]),
    ('hybrid', 5, 0.438, [9, 10]),
    ('hybrid', 7, 0.535, [8]),
]

print("=" * 85)
print("Sprint 127a: Extended exact chi_F to GPU sizes")
print("=" * 85)

# Phase 1: Compute new data points
for model, q, g_c, new_sizes in new_configs:
    label = f"{model}_q{q}"
    print(f"\n{'='*70}")
    print(f"  {model.upper()} q={q}, g_c={g_c:.4f} — NEW sizes: {new_sizes}")
    print(f"{'='*70}")

    for n in new_sizes:
        dim = q**n
        print(f"  n={n:2d} dim={dim:>12,} ...", end="", flush=True)
        t0 = time.time()

        if model == 'sq':
            H_coup, H_field = build_sq_potts_parts(n, q)
        else:
            H_coup, H_field = build_hybrid_parts(n, q)

        chi_val = chi_F_exact(H_coup, H_field, g_c, n)
        dt = time.time() - t0

        print(f"  chi_F = {chi_val:12.6f}  ({dt:.1f}s)")

        # Store
        key = (model, q)
        if key not in prior_data:
            prior_data[key] = {}
        prior_data[key][n] = chi_val

        if label not in results['new_points']:
            results['new_points'][label] = []
        results['new_points'][label].append({
            'n': n, 'dim': dim, 'chi_F_exact': float(chi_val), 'time_s': round(dt, 1)
        })

        # Record to DB
        record(sprint=127, model=model, q=q, n=n,
               quantity='chi_F_exact', value=float(chi_val),
               method='finite_diff')
        save()

# Phase 2: Combined power-law fits with ALL data
print(f"\n{'='*85}")
print("Phase 2: Combined power-law fits (all sizes, exact chi_F)")
print(f"{'='*85}")

g_c_map = {'sq': {3: 1.0/3, 4: 0.25, 5: 0.20, 7: 1.0/7},
            'hybrid': {3: 0.333, 4: 0.393, 5: 0.438, 7: 0.535}}

for (model, q), data in sorted(prior_data.items()):
    label = f"{model}_q{q}"
    sizes = np.array(sorted(data.keys()), dtype=float)
    chi_vals = np.array([data[int(n)] for n in sizes])

    print(f"\n  {model.upper()} q={q}: sizes={[int(s) for s in sizes]}")
    for i, n in enumerate(sizes):
        print(f"    n={int(n):2d}: chi_F={chi_vals[i]:.6f}")

    if len(sizes) >= 3:
        fit = fit_power_law(sizes, chi_vals)
        print(f"  -> alpha = {fit['alpha']:.4f} +/- {fit['alpha_err']:.4f}  "
              f"(R2={fit['r_squared']:.8f}, A={fit['A']:.4f})")

        pairs = pairwise_exponents(sizes, chi_vals)
        print(f"    pairwise: ", end="")
        for p in pairs:
            print(f"({p['n1']},{p['n2']})={p['alpha']:.4f}  ", end="")
        print()

        # Store fit
        results['combined_fits'][label] = {
            'alpha': fit['alpha'], 'alpha_err': fit['alpha_err'],
            'A': fit['A'], 'r_squared': fit['r_squared'],
            'sizes': [int(s) for s in sizes],
            'chi_F': [float(v) for v in chi_vals],
            'pairwise': [{'n1': p['n1'], 'n2': p['n2'], 'alpha': p['alpha']} for p in pairs],
        }

        # Record to DB
        record(sprint=127, model=model, q=q, n=int(max(sizes)),
               quantity='alpha_exact', value=fit['alpha'],
               error=fit['alpha_err'], method='exact_chif_power_law',
               notes=f"sizes={[int(s) for s in sizes]}")
    else:
        print(f"  -> only {len(sizes)} sizes, need 3+ for fit")

save()

# Phase 3: alpha(q) refit
print(f"\n{'='*85}")
print("Phase 3: alpha(q) curve refit with extended exact chi_F data")
print(f"{'='*85}")

for model_name in ['sq', 'hybrid']:
    qs = []
    alphas = []
    errs = []
    for (m, q), data in sorted(prior_data.items()):
        if m != model_name:
            continue
        sizes = np.array(sorted(data.keys()), dtype=float)
        chi_vals = np.array([data[int(n)] for n in sizes])
        if len(sizes) >= 3:
            fit = fit_power_law(sizes, chi_vals)
            qs.append(q)
            alphas.append(fit['alpha'])
            errs.append(fit['alpha_err'])

    qs = np.array(qs, dtype=float)
    alphas = np.array(alphas)
    errs = np.array(errs)

    print(f"\n  {model_name.upper()} alpha(q):")
    for i, q in enumerate(qs):
        print(f"    q={int(q)}: alpha={alphas[i]:.4f} +/- {errs[i]:.4f}")

    if len(qs) >= 3:
        # Fit: alpha = a * ln(q) + b
        log_q = np.log(qs)
        weights = 1.0 / errs
        # Weighted linear regression in log(q)
        W = np.diag(weights**2)
        X = np.column_stack([log_q, np.ones_like(log_q)])
        beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ alphas)
        residuals = alphas - X @ beta
        chi2 = np.sum((residuals * weights)**2)
        print(f"  -> Log fit: alpha(q) = {beta[0]:.4f} * ln(q) + ({beta[1]:.4f})")
        print(f"    chi2/dof = {chi2/(len(qs)-2):.4f}")
        predictions = X @ beta
        for i, q in enumerate(qs):
            print(f"    q={int(q)}: measured={alphas[i]:.4f}, predicted={predictions[i]:.4f}, "
                  f"resid={residuals[i]:.4f}")

        results['combined_fits'][f'{model_name}_alpha_q_logfit'] = {
            'a': float(beta[0]), 'b': float(beta[1]),
            'formula': f'alpha(q) = {beta[0]:.4f} * ln(q) + ({beta[1]:.4f})',
            'qs': [int(q) for q in qs],
            'alphas': [float(a) for a in alphas],
            'errors': [float(e) for e in errs],
            'chi2_dof': float(chi2/(len(qs)-2)),
        }

        # Also try: alpha = a * ln(q) + b * ln(ln(q)) + c (loglog model from Sprint 116)
        if all(q > 1 for q in qs):
            log_log_q = np.log(np.log(qs))
            X2 = np.column_stack([log_q, log_log_q, np.ones_like(log_q)])
            if len(qs) > 3:
                beta2 = np.linalg.solve(X2.T @ W @ X2, X2.T @ W @ alphas)
                pred2 = X2 @ beta2
                resid2 = alphas - pred2
                chi2_2 = np.sum((resid2 * weights)**2)
                print(f"  -> Log+loglog fit: alpha(q) = {beta2[0]:.4f}*ln(q) + {beta2[1]:.4f}*ln(ln(q)) + ({beta2[2]:.4f})")
                print(f"    chi2/dof = {chi2_2/(len(qs)-3):.4f}")
                results['combined_fits'][f'{model_name}_alpha_q_loglogfit'] = {
                    'a': float(beta2[0]), 'b': float(beta2[1]), 'c': float(beta2[2]),
                    'chi2_dof': float(chi2_2/(len(qs)-3)) if len(qs) > 3 else None,
                }

save()

# Summary table
print(f"\n{'='*85}")
print("SUMMARY: Extended exact chi_F exponents (Sprint 127)")
print(f"{'='*85}")
print(f"{'Model':>8} {'q':>3} {'alpha':>14} {'#sizes':>7} {'largest_n':>10} {'R²':>10}")
for label, fit in sorted(results['combined_fits'].items()):
    if 'alpha' not in fit:
        continue
    parts = label.split('_')
    if len(parts) == 2 and parts[1].startswith('q'):
        model = parts[0]
        q_val = parts[1]
        n_sizes = len(fit.get('sizes', []))
        max_n = max(fit.get('sizes', [0]))
        alpha_str = f"{fit['alpha']:.4f}+/-{fit['alpha_err']:.4f}"
        print(f"{model:>8} {q_val:>3} {alpha_str:>14} {n_sizes:>7} {max_n:>10} {fit['r_squared']:>10.8f}")

# Comparison with Sprint 126
print(f"\n{'='*85}")
print("COMPARISON: Sprint 126 vs Sprint 127 alpha values")
print(f"{'='*85}")
s126_alphas = {
    ('sq', 3): (1.481, 0.014), ('sq', 4): (1.800, 0.010),
    ('sq', 5): (2.094, 0.005), ('sq', 7): (2.606, 0.019),
    ('hybrid', 3): (1.481, 0.014), ('hybrid', 4): (1.557, 0.019),
    ('hybrid', 5): (1.436, 0.043), ('hybrid', 7): (1.028, 0.067),
}
for label, fit in sorted(results['combined_fits'].items()):
    if 'alpha' not in fit:
        continue
    parts = label.split('_')
    if len(parts) != 2 or not parts[1].startswith('q'):
        continue
    model = parts[0]
    q = int(parts[1][1:])
    key = (model, q)
    if key in s126_alphas:
        old_a, old_e = s126_alphas[key]
        new_a, new_e = fit['alpha'], fit['alpha_err']
        shift = new_a - old_a
        print(f"  {model:>8} q={q}: {old_a:.4f}+/-{old_e:.4f} -> {new_a:.4f}+/-{new_e:.4f}  "
              f"(shift={shift:+.4f}, err reduction={1-new_e/old_e:.0%})")

save()
print(f"\nResults saved to results/sprint_127a_extended_exact_chif.json")
