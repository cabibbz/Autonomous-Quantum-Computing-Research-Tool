"""Sprint 113c: Periodic vs Open BC chi_F comparison.

113a/b found open BC alpha << periodic BC alpha (q=4: 1.51 vs 1.77, q=5: 1.60 vs 2.09).
Is this a BC effect or a method artifact?

Test: compute overlap-based chi_F with PERIODIC BC exact diag at same sizes.
If periodic overlap matches periodic spectral, the method is correct and
the difference is purely a boundary condition effect.

Also: quantify boundary correction by fitting chi_F = A*N^alpha + B.
"""
import numpy as np
import json, time, os
from scipy.sparse import coo_matrix, csr_matrix
from scipy.optimize import curve_fit
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '113c_bc_comparison',
    'sprint': 113,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'periodic': {},
    'open': {},
    'comparison': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_113c_bc_comparison.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)


def build_sq_potts(n, q_val, g, periodic=True):
    """Build S_q Potts Hamiltonian (periodic or open BC)."""
    dim = q_val**n
    all_idx = np.arange(dim, dtype=np.int64)
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q_val
        tmp //= q_val
    powers = q_val ** np.arange(n, dtype=np.int64)
    n_bonds = n if periodic else (n - 1)
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n_bonds):
        nxt = (site + 1) % n
        diag_vals -= (digits[:, site] == digits[:, nxt]).astype(np.float64)
    rows_list = [all_idx]; cols_list = [all_idx]; vals_list = [diag_vals]
    for site in range(n):
        pw = powers[site]
        old_digit = digits[:, site]
        for k in range(1, q_val):
            new_digit = (old_digit + k) % q_val
            delta = (new_digit.astype(np.int64) - old_digit.astype(np.int64)) * pw
            new_idx = all_idx + delta
            rows_list.append(new_idx)
            cols_list.append(all_idx)
            vals_list.append(np.full(dim, -g, dtype=np.float64))
    return csr_matrix(coo_matrix((np.concatenate(vals_list),
                                   (np.concatenate(rows_list), np.concatenate(cols_list))),
                                  shape=(dim, dim)))


def overlap_chif(n, q, g_c, dg=1e-3, periodic=True):
    """Compute chi_F via ground state overlap."""
    H0 = build_sq_potts(n, q, g_c, periodic=periodic)
    H1 = build_sq_potts(n, q, g_c + dg, periodic=periodic)
    _, evecs0 = eigsh(H0, k=1, which='SA')
    _, evecs1 = eigsh(H1, k=1, which='SA')
    psi0 = evecs0[:, 0]
    psi1 = evecs1[:, 0]
    overlap = abs(np.dot(psi0, psi1))**2
    chi_F_total = (1.0 - overlap) / dg**2
    return chi_F_total / n, overlap


def spectral_chif(n, q, g_c, k_states=10):
    """Compute chi_F via spectral decomposition (periodic BC) — for cross-check."""
    dim = q**n
    # Build coupling and field separately
    all_idx = np.arange(dim, dtype=np.int64)
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q_val
        tmp //= q_val
    powers = q_val ** np.arange(n, dtype=np.int64)
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n):
        nxt = (site + 1) % n
        diag_vals -= (digits[:, site] == digits[:, nxt]).astype(np.float64)
    H_coup = csr_matrix((diag_vals, (all_idx, all_idx)), shape=(dim, dim))
    rows_list, cols_list, vals_list = [], [], []
    for site in range(n):
        pw = powers[site]
        old_digit = digits[:, site]
        for k in range(1, q_val):
            new_digit = (old_digit + k) % q_val
            delta = (new_digit.astype(np.int64) - old_digit.astype(np.int64)) * pw
            new_idx = all_idx + delta
            rows_list.append(new_idx)
            cols_list.append(all_idx)
            vals_list.append(np.full(dim, -1.0, dtype=np.float64))
    H_field = coo_matrix((np.concatenate(vals_list),
                           (np.concatenate(rows_list), np.concatenate(cols_list))),
                          shape=(dim, dim)).tocsr()
    H = H_coup + g_c * H_field
    k_use = min(k_states, dim - 2)
    evals, evecs = eigsh(H, k=k_use, which='SA')
    order = np.argsort(evals)
    evals = evals[order]; evecs = evecs[:, order]
    E0 = evals[0]; psi0 = evecs[:, 0]
    H_field_psi0 = H_field.dot(psi0)
    total_chi = 0.0
    for i in range(1, k_use):
        gap = evals[i] - E0
        if gap < 1e-12: continue
        me = np.dot(evecs[:, i], H_field_psi0)
        total_chi += me**2 / gap**2
    return total_chi / n


# ============================================================
dg = 1e-3

print("=" * 70)
print("Sprint 113c: Periodic vs Open BC chi_F comparison")
print("=" * 70)

for q_val in [2, 4, 5]:
    g_c = 1.0 / q_val
    print(f"\n{'='*60}")
    print(f"q={q_val}, g_c={g_c:.4f}")
    print(f"{'='*60}")

    # Determine feasible sizes
    if q_val == 2:
        sizes = [6, 8, 10, 12, 14]
    elif q_val == 4:
        sizes = [6, 7, 8, 9]
    elif q_val == 5:
        sizes = [6, 7, 8]

    periodic_data = []
    open_data = []

    for n in sizes:
        dim = q_val**n
        if dim > 500000:
            print(f"  n={n}: dim={dim:,}, skipping (too large)")
            continue

        # Periodic BC overlap
        t0 = time.time()
        chi_per, ov_per = overlap_chif(n, q_val, g_c, dg=dg, periodic=True)
        dt_per = time.time() - t0

        # Open BC overlap
        t0 = time.time()
        chi_open, ov_open = overlap_chif(n, q_val, g_c, dg=dg, periodic=False)
        dt_open = time.time() - t0

        ratio = chi_per / chi_open if chi_open > 0 else 0
        print(f"  n={n:2d}: periodic={chi_per:.6f}, open={chi_open:.6f}, "
              f"ratio={ratio:.4f}, t={dt_per+dt_open:.1f}s")

        periodic_data.append({'n': n, 'chi_F': float(chi_per)})
        open_data.append({'n': n, 'chi_F': float(chi_open)})

    results['periodic'][str(q_val)] = periodic_data
    results['open'][str(q_val)] = open_data
    save()

    # Pairwise alpha
    for label, data in [("periodic", periodic_data), ("open", open_data)]:
        if len(data) < 2:
            continue
        Ns = np.array([d['n'] for d in data], dtype=float)
        chis = np.array([d['chi_F'] for d in data])
        log_N = np.log(Ns)
        log_chi = np.log(chis)

        pw = []
        for i in range(len(Ns) - 1):
            a = (log_chi[i+1] - log_chi[i]) / (log_N[i+1] - log_N[i])
            pw.append(a)

        p = np.polyfit(log_N, log_chi, 1)
        alpha_g = p[0]
        print(f"    {label:10s}: pairwise = {['%.3f' % a for a in pw]}, global = {alpha_g:.4f}")

        results['comparison'][f'q{q_val}_{label}'] = {
            'pairwise_alpha': [float(a) for a in pw],
            'alpha_global': float(alpha_g),
        }

    # Periodic/Open ratio scaling
    if len(periodic_data) >= 2 and len(open_data) >= 2:
        Ns_match = [d['n'] for d in periodic_data]
        ratios = [periodic_data[i]['chi_F'] / open_data[i]['chi_F']
                  for i in range(len(periodic_data))]
        print(f"    ratio periodic/open: {['%.3f' % r for r in ratios]}")
        results['comparison'][f'q{q_val}_ratio'] = {
            'sizes': Ns_match, 'ratios': [float(r) for r in ratios],
        }

save()

# === Summary ===
print(f"\n{'='*70}")
print("SUMMARY: Alpha comparison (periodic vs open BC)")
print("=" * 70)
print(f"  {'q':>3} {'periodic':>12} {'open':>12} {'diff':>8} {'prior_periodic':>16}")
print(f"  {'-'*3} {'-'*12} {'-'*12} {'-'*8} {'-'*16}")

prior = {2: 1.0, 4: 1.77, 5: 2.09}
for q_val in [2, 4, 5]:
    pk = f'q{q_val}_periodic'
    ok = f'q{q_val}_open'
    if pk in results['comparison'] and ok in results['comparison']:
        a_p = results['comparison'][pk]['alpha_global']
        a_o = results['comparison'][ok]['alpha_global']
        diff = a_p - a_o
        p = prior.get(q_val, 0)
        print(f"  {q_val:3d} {a_p:12.4f} {a_o:12.4f} {diff:8.4f} {p:16.2f}")

print(f"\nResults saved.")
