"""Sprint 128a: S_q q=6 exact chi_F across sizes n=4-9.

g_c = 1/6 (exact, Kramers-Wannier self-duality).
Sizes: n=4 (1296), n=5 (7776), n=6 (46656), n=7 (280k), n=8 (1.7M GPU), n=9 (10M GPU).
Time estimate: n=8 ~5s GPU, n=9 ~42s GPU (based on q=5 n=10 timing).
"""
import numpy as np
import json, time, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from gpu_utils import eigsh
from hamiltonian_utils import build_sq_potts_parts
from fss_utils import fit_power_law, pairwise_exponents
from db_utils import record

q = 6
g_c = 1.0 / q
results = {
    'experiment': '128a_sq_q6_chif',
    'sprint': 128,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'model': 'sq',
    'q': q,
    'g_c': g_c,
    'data': [],
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results',
                           'sprint_128a_sq_q6_chif.json')
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


sizes = [4, 5, 6, 7, 8, 9]
print("=" * 80)
print(f"Sprint 128a: S_q q={q} exact chi_F, g_c={g_c:.6f}")
print("=" * 80)

# Timing test on smallest size first
n_test = sizes[0]
dim_test = q**n_test
print(f"\nTiming test: n={n_test}, dim={dim_test:,}")
t0 = time.time()
H_coup, H_field = build_sq_potts_parts(n_test, q)
dt_build = time.time() - t0
print(f"  Build time: {dt_build:.2f}s")
t0 = time.time()
chi_test = chi_F_exact(H_coup, H_field, g_c, n_test)
dt_solve = time.time() - t0
print(f"  Solve time: {dt_solve:.2f}s, chi_F={chi_test:.6f}")

chi_data = {}

for n in sizes:
    dim = q**n
    print(f"\nn={n:2d}, dim={dim:>12,} ...", flush=True)

    t0 = time.time()
    H_coup, H_field = build_sq_potts_parts(n, q)
    dt_build = time.time() - t0
    print(f"  Build: {dt_build:.1f}s", end="", flush=True)

    t0 = time.time()
    chi_val = chi_F_exact(H_coup, H_field, g_c, n)
    dt_solve = time.time() - t0
    print(f"  Solve: {dt_solve:.1f}s  chi_F = {chi_val:.6f}")

    chi_data[n] = chi_val
    results['data'].append({
        'n': n, 'dim': dim,
        'chi_F_exact': float(chi_val),
        'build_time_s': round(dt_build, 1),
        'solve_time_s': round(dt_solve, 1),
    })

    record(sprint=128, model='sq', q=q, n=n,
           quantity='chi_F_exact', value=float(chi_val),
           method='finite_diff')
    save()

    # Abort if taking too long (safety for n=9)
    if dt_solve > 250:
        print(f"  WARNING: solve took {dt_solve:.0f}s, skipping larger sizes")
        break

# Power-law fit
print(f"\n{'='*80}")
print("Power-law fit: chi_F = A * N^alpha")
print(f"{'='*80}")

ns = np.array(sorted(chi_data.keys()), dtype=float)
chis = np.array([chi_data[int(n)] for n in ns])

for i, n in enumerate(ns):
    print(f"  n={int(n):2d}: chi_F = {chis[i]:.6f}")

if len(ns) >= 3:
    fit = fit_power_law(ns, chis)
    print(f"\nalpha = {fit['alpha']:.4f} +/- {fit['alpha_err']:.4f}")
    print(f"A = {fit['A']:.4f}, R² = {fit['r_squared']:.8f}")

    pairs = pairwise_exponents(ns, chis)
    print("\nPairwise exponents:")
    for p in pairs:
        print(f"  ({p['n1']},{p['n2']}): alpha = {p['alpha']:.4f}")

    results['fit'] = {
        'alpha': fit['alpha'], 'alpha_err': fit['alpha_err'],
        'A': fit['A'], 'r_squared': fit['r_squared'],
        'pairwise': [{'n1': p['n1'], 'n2': p['n2'], 'alpha': p['alpha']} for p in pairs],
    }

    record(sprint=128, model='sq', q=q, n=int(max(ns)),
           quantity='alpha_exact', value=fit['alpha'],
           error=fit['alpha_err'], method='exact_chif_power_law',
           notes=f"sizes={[int(s) for s in ns]}")

save()
print(f"\nResults saved.")
