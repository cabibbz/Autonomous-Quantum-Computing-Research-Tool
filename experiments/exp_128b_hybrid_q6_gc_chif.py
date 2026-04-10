"""Sprint 128b: Hybrid q=6 — find g_c via gap×N crossing, then compute chi_F.

Known hybrid g_c: q=5→0.438, q=7→0.535. Interpolation: q=6 ≈ 0.486.
Phase 1: Scan g ∈ [0.45, 0.52] at n=6,8 to find gap×N crossing.
Phase 2: Compute exact chi_F at g_c for n=4-8.
"""
import numpy as np
import json, time, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from gpu_utils import eigsh
from hamiltonian_utils import build_hybrid_parts
from fss_utils import fit_power_law, pairwise_exponents
from db_utils import record

q = 6
results = {
    'experiment': '128b_hybrid_q6_gc_chif',
    'sprint': 128,
    'model': 'hybrid',
    'q': q,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'gc_scan': {},
    'gc_estimate': None,
    'chi_F_data': [],
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results',
                           'sprint_128b_hybrid_q6_gc_chif.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)


def chi_F_exact(H_coup, H_field, g_c, n, dg=1e-4):
    """Exact chi_F via central finite-difference."""
    Hp = H_coup + (g_c + dg) * H_field
    Hm = H_coup + (g_c - dg) * H_field
    H0 = H_coup + g_c * H_field
    _, vecs0 = eigsh(H0, k=1, which='SA')
    _, vecsp = eigsh(Hp, k=1, which='SA')
    _, vecsm = eigsh(Hm, k=1, which='SA')
    psi0 = vecs0[:, 0]; psip = vecsp[:, 0]; psim = vecsm[:, 0]
    if np.dot(psi0, psip) < 0: psip = -psip
    if np.dot(psi0, psim) < 0: psim = -psim
    ov_p = np.dot(psi0, psip)**2
    ov_m = np.dot(psi0, psim)**2
    return (2.0 - ov_p - ov_m) / (dg**2 * n)


# ============================================================
# Phase 1: Find g_c via gap×N crossing
# ============================================================
print("=" * 80)
print(f"Sprint 128b: Hybrid q={q} — g_c scan + chi_F")
print("=" * 80)

g_values = np.linspace(0.44, 0.54, 25)
scan_sizes = [6, 8]  # 6^6=46656, 6^8=1.68M

print("\nPhase 1: gap×N crossing scan")
gap_data = {}

for n in scan_sizes:
    dim = q**n
    print(f"\n  n={n}, dim={dim:,}")
    gap_data[n] = []

    for g in g_values:
        H_coup, H_field = build_hybrid_parts(n, q)
        H = H_coup + g * H_field
        evals, _ = eigsh(H, k=2, which='SA')
        gap = evals[1] - evals[0]
        gap_N = gap * n
        gap_data[n].append((float(g), float(gap), float(gap_N)))
        print(f"    g={g:.4f}: gap={gap:.6f}, gap×N={gap_N:.4f}")

    results['gc_scan'][str(n)] = gap_data[n]
    save()

# Find crossing: where gap×N curves for n=6 and n=8 cross
g_vals = np.array([d[0] for d in gap_data[scan_sizes[0]]])
gapN_small = np.array([d[2] for d in gap_data[scan_sizes[0]]])
gapN_large = np.array([d[2] for d in gap_data[scan_sizes[1]]])

diff = gapN_small - gapN_large
# Find zero crossing
sign_changes = np.where(np.diff(np.sign(diff)))[0]
if len(sign_changes) > 0:
    idx = sign_changes[0]
    # Linear interpolation
    g1, g2 = g_vals[idx], g_vals[idx + 1]
    d1, d2 = diff[idx], diff[idx + 1]
    g_c = g1 - d1 * (g2 - g1) / (d2 - d1)
    print(f"\n  gap*N crossing at g_c ~ {g_c:.4f} (between g={g1:.4f} and g={g2:.4f})")
else:
    # Fallback: interpolation from known values
    g_c = 0.486
    print(f"\n  No crossing found in range. Using interpolation g_c ≈ {g_c:.4f}")

results['gc_estimate'] = float(g_c)
record(sprint=128, model='hybrid', q=q, n=max(scan_sizes),
       quantity='g_c', value=float(g_c), method='gap_N_crossing')
save()

# ============================================================
# Phase 2: Exact chi_F at g_c
# ============================================================
print(f"\n{'='*80}")
print(f"Phase 2: Exact chi_F at g_c={g_c:.4f}")
print(f"{'='*80}")

chi_sizes = [4, 5, 6, 7, 8]
chi_data = {}

for n in chi_sizes:
    dim = q**n
    print(f"\n  n={n:2d}, dim={dim:>10,} ...", end="", flush=True)
    t0 = time.time()
    H_coup, H_field = build_hybrid_parts(n, q)
    chi_val = chi_F_exact(H_coup, H_field, g_c, n)
    dt = time.time() - t0
    print(f"  chi_F = {chi_val:.6f}  ({dt:.1f}s)")

    chi_data[n] = chi_val
    results['chi_F_data'].append({
        'n': n, 'dim': dim, 'chi_F_exact': float(chi_val), 'time_s': round(dt, 1)
    })
    record(sprint=128, model='hybrid', q=q, n=n,
           quantity='chi_F_exact', value=float(chi_val),
           method='finite_diff', notes=f'g_c={g_c:.4f}')
    save()

# Fit
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

    record(sprint=128, model='hybrid', q=q, n=int(max(ns)),
           quantity='alpha_exact', value=fit['alpha'],
           error=fit['alpha_err'], method='exact_chif_power_law',
           notes=f"sizes={[int(s) for s in ns]}, g_c={g_c:.4f}")

save()
print(f"\nResults saved.")
