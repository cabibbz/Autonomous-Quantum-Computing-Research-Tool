"""Sprint 125a: Is single-multiplet dominance a k-truncation artifact?

Compare EXACT chi_F (finite-difference) vs spectral chi_F at k=12,50,100,200.
Test at q=3,4,5 for sizes where both methods are feasible.

The spectral decomposition: chi_F = (1/N) * sum_{n>0} |<n|V|0>|^2 / (E_n-E_0)^2
Only k states are computed. If frac=1.000 is real, increasing k shouldn't change anything.
If it's a truncation artifact, we'll see the captured fraction drop below 1.000.

Ground truth: chi_F_exact = (1/N) * (2*|<psi0|psi0'>|^2 - 2) / dg^2
via finite-difference of ground state overlap at g_c +/- dg.
"""
import numpy as np
import json, time, os
from scipy.sparse import csr_matrix
from gpu_utils import eigsh
from hamiltonian_utils import build_sq_potts_parts, build_hybrid_parts
from db_utils import record

results = {
    'experiment': '125a_truncation_test',
    'sprint': 125,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results',
                           'sprint_125a_truncation_test.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)


def chi_F_finite_diff(H_coup, H_field, g_c, n, dg=1e-4):
    """Exact chi_F via finite-difference of ground state overlap."""
    H0 = H_coup + g_c * H_field
    Hp = H_coup + (g_c + dg) * H_field
    Hm = H_coup + (g_c - dg) * H_field

    _, vecs0 = eigsh(H0, k=1, which='SA')
    _, vecsp = eigsh(Hp, k=1, which='SA')
    _, vecsm = eigsh(Hm, k=1, which='SA')

    psi0 = vecs0[:, 0]
    psip = vecsp[:, 0]
    psim = vecsm[:, 0]

    # Fix sign ambiguity
    if np.dot(psi0, psip) < 0: psip = -psip
    if np.dot(psi0, psim) < 0: psim = -psim

    overlap_p = np.dot(psi0, psip)**2
    overlap_m = np.dot(psi0, psim)**2

    # chi_F = -d^2 ln|<psi0|psi(g)>|^2 / dg^2 / N
    # ≈ (2 - overlap_p - overlap_m) / dg^2 / N  (central difference)
    chi_F = (2.0 - overlap_p - overlap_m) / (dg**2 * n)
    return chi_F


def chi_F_spectral(H_coup, H_field, g_c, n, k_vals):
    """Spectral chi_F at multiple k values. Returns dict {k: {chi_F, frac, dominant_me2, dominant_gap}}."""
    H = H_coup + g_c * H_field
    dim = H.shape[0]
    k_max = max(k_vals)
    k_use = min(k_max, dim - 2)

    evals, evecs = eigsh(H, k=k_use, which='SA')
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]

    E0 = evals[0]
    psi0 = evecs[:, 0]
    V_psi0 = H_field.dot(psi0)

    # Compute all contributions
    contributions = []
    for i in range(1, k_use):
        gap = evals[i] - E0
        if gap < 1e-12:
            continue
        me = np.dot(evecs[:, i], V_psi0)
        chi_n = me**2 / gap**2
        contributions.append(chi_n)

    contributions.sort(reverse=True)

    out = {}
    for k in k_vals:
        k_eff = min(k - 1, len(contributions))  # k includes ground state
        total = sum(contributions[:k_eff])
        chi_F_k = total / n
        best = contributions[0] if contributions else 0
        frac = best / total if total > 0 else 0
        out[k] = {
            'chi_F': float(chi_F_k),
            'total_sum': float(total),
            'dominant_contrib': float(best),
            'frac': float(frac),
            'n_nonzero': sum(1 for c in contributions[:k_eff] if c > 1e-20),
        }
    return out


# Test configurations
configs = [
    # (model, q, g_c, sizes, source)
    ('sq', 3, 1.0/3, [4, 6, 8, 10], 'S_q q=3, exact alpha=1.40'),
    ('sq', 4, 0.25, [4, 6, 8], 'S_q q=4, alpha~1.77'),
    ('sq', 5, 0.20, [4, 6, 8], 'S_q q=5, walking alpha~2.09'),
    ('hybrid', 5, 0.438, [4, 6, 8], 'Hybrid q=5, alpha~1.41'),
]

k_vals = [12, 50, 100, 200]

print("=" * 80)
print("Sprint 125a: TRUNCATION TEST — Is frac=1.000 real or artifact?")
print("Compare exact chi_F (finite-diff) vs spectral at k=12,50,100,200")
print("=" * 80)

for model, q, g_c, sizes, desc in configs:
    print(f"\n{'='*70}")
    print(f"  {desc}")
    print(f"  Model={model}, q={q}, g_c={g_c}")
    print(f"{'='*70}")

    for n in sizes:
        dim = q**n
        # Limit k to what's feasible
        k_test = [k for k in k_vals if k < dim - 1]
        if not k_test:
            k_test = [min(12, dim - 2)]

        print(f"\n  n={n} (dim={dim:,}), testing k={k_test}...")
        t0 = time.time()

        # Build Hamiltonian
        if model == 'sq':
            H_coup, H_field = build_sq_potts_parts(n, q)
        else:
            H_coup, H_field = build_hybrid_parts(n, q)

        # Exact chi_F
        chi_exact = chi_F_finite_diff(H_coup, H_field, g_c, n)

        # Spectral chi_F at multiple k
        spectral = chi_F_spectral(H_coup, H_field, g_c, n, k_test)

        dt = time.time() - t0

        print(f"    EXACT chi_F = {chi_exact:.6f}")
        print(f"    {'k':>5s}  {'chi_F_spec':>12s}  {'captured%':>10s}  {'frac_dom':>10s}  {'n_nonzero':>10s}")
        for k in k_test:
            s = spectral[k]
            captured = s['chi_F'] / chi_exact * 100 if chi_exact > 0 else 0
            print(f"    {k:5d}  {s['chi_F']:12.6f}  {captured:9.2f}%  {s['frac']:10.6f}  {s['n_nonzero']:10d}")

        # Store
        key = f"{model}_q{q}_n{n}"
        results['data'][key] = {
            'model': model, 'q': q, 'n': n, 'dim': dim, 'g_c': g_c,
            'chi_F_exact': float(chi_exact),
            'spectral': {str(k): spectral[k] for k in k_test},
            'time_s': round(dt, 1),
        }
        save()

        # Record key result
        best_k = max(k_test)
        captured_best = spectral[best_k]['chi_F'] / chi_exact * 100 if chi_exact > 0 else 0
        record(sprint=125, model=model, q=q, n=n,
               quantity='chi_F_exact', value=float(chi_exact),
               method='finite_diff')
        record(sprint=125, model=model, q=q, n=n,
               quantity='captured_pct_k200', value=float(captured_best),
               method=f'spectral_k{best_k}')

        print(f"    ({dt:.1f}s)")

# Summary
print(f"\n{'='*80}")
print("SUMMARY: Dominant state capture of TRUE chi_F")
print(f"{'='*80}")
print(f"{'Config':>25s}  {'n':>3s}  {'Exact':>10s}  {'k=12':>8s}  {'k=best':>8s}  {'captured':>9s}")
for model, q, g_c, sizes, desc in configs:
    for n in sizes:
        key = f"{model}_q{q}_n{n}"
        d = results['data'].get(key)
        if not d:
            continue
        exact = d['chi_F_exact']
        s12 = d['spectral'].get('12', {})
        best_k_str = max(d['spectral'].keys(), key=int)
        s_best = d['spectral'][best_k_str]
        cap12 = s12.get('chi_F', 0) / exact * 100 if exact > 0 else 0
        cap_best = s_best['chi_F'] / exact * 100 if exact > 0 else 0
        label = f"{model} q={q}"
        print(f"{label:>25s}  {n:3d}  {exact:10.4f}  {cap12:7.1f}%  {cap_best:7.1f}%  (k={best_k_str})")

save()
print(f"\nResults saved.")
