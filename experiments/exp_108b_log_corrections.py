"""Sprint 108b: Test for BKT log corrections in q=4 spectral exponents.

q=4 Potts is BKT-like and should have multiplicative log corrections.
If gap_m ~ N^{-z_m} * (ln N)^p, then pairwise z_m drifts as ~p/ln(N).

Compare pairwise drift in z_m and beta_me across q=3 (continuous), q=4 (BKT), q=5 (walking).
Use existing 108a data for q=4, compute fresh for q=3 (n=4-10) and q=5 (n=4-10).
"""
import numpy as np
import json, time, os
from scipy.sparse import lil_matrix, csr_matrix
from scipy.optimize import curve_fit
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '108b_log_corrections',
    'sprint': 108,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_108b_log_corrections.json')
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

def get_spectral_data(q, sizes, k_need):
    """Get gap_m and |me|^2 for each size."""
    g_c = 1.0 / q
    gap_ms = []
    me_sqs = []
    chi_Fs = []
    Ns = []

    for n in sizes:
        dim = q**n
        if dim > 12_000_000:
            print(f"  q={q}, n={n}: dim={dim:,} too large, skipping")
            continue
        print(f"  q={q}, n={n} (dim={dim:,})...", end="", flush=True)
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

        best_chi = 0
        best_gap = 0
        best_me_sq = 0
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

        chi_F = total_chi / n
        dt = time.time() - t0
        print(f" gap_m={best_gap:.6f}, |me|²={best_me_sq:.4f}, chi_F={chi_F:.4f} ({dt:.1f}s)")

        gap_ms.append(best_gap)
        me_sqs.append(best_me_sq)
        chi_Fs.append(chi_F)
        Ns.append(n)

    return np.array(Ns), np.array(gap_ms), np.array(me_sqs), np.array(chi_Fs)

print("=" * 70)
print("Sprint 108b: BKT log correction test for z_m and beta_me")
print("=" * 70)

# q=3: continuous, no log corrections expected
print("\n--- q=3 (continuous CFT, no log corrections) ---")
N3, gap3, me3, chi3 = get_spectral_data(3, [4, 5, 6, 7, 8, 9, 10], k_need=6)

# q=4: BKT, log corrections expected
print("\n--- q=4 (BKT boundary, log corrections expected) ---")
# Load from 108a instead of recomputing
try:
    with open('results/sprint_108a_q4_spectral_extend.json') as f:
        d108a = json.load(f)
    N4 = []
    gap4 = []
    me4 = []
    chi4 = []
    for key in sorted(d108a['data'].keys(), key=int):
        dd = d108a['data'][key]
        N4.append(dd['n'])
        gap4.append(dd['gap_multiplet'])
        me4.append(dd['me_sq'])
        chi4.append(dd['chi_F'])
    N4 = np.array(N4)
    gap4 = np.array(gap4)
    me4 = np.array(me4)
    chi4 = np.array(chi4)
    print(f"  Loaded {len(N4)} sizes from 108a")
except FileNotFoundError:
    print("  108a data not found, computing fresh...")
    N4, gap4, me4, chi4 = get_spectral_data(4, [4, 5, 6, 7, 8, 9, 10], k_need=8)

# q=5: walking, no log corrections expected
print("\n--- q=5 (walking regime) ---")
N5, gap5, me5, chi5 = get_spectral_data(5, [4, 5, 6, 7, 8, 9], k_need=8)

save()

# Analysis: pairwise exponent drift
print(f"\n{'=' * 70}")
print("PAIRWISE EXPONENT DRIFT ANALYSIS")
print("=" * 70)

for label, Ns, gaps, mes, chis in [
    ("q=3", N3, gap3, me3, chi3),
    ("q=4", N4, gap4, me4, chi4),
    ("q=5", N5, gap5, me5, chi5),
]:
    if len(Ns) < 3:
        continue
    log_N = np.log(Ns.astype(float))
    log_gap = np.log(gaps)
    log_me = np.log(mes)
    log_chi = np.log(chis)

    # Pairwise exponents
    pw_z = []
    pw_b = []
    pw_a = []
    pw_N_mid = []
    for i in range(len(Ns)-1):
        dln = log_N[i+1] - log_N[i]
        pw_z.append(-(log_gap[i+1] - log_gap[i]) / dln)
        pw_b.append((log_me[i+1] - log_me[i]) / dln)
        pw_a.append((log_chi[i+1] - log_chi[i]) / dln)
        pw_N_mid.append(0.5 * (Ns[i] + Ns[i+1]))

    pw_z = np.array(pw_z)
    pw_b = np.array(pw_b)
    pw_a = np.array(pw_a)
    pw_N_mid = np.array(pw_N_mid)

    # Drift = change in pairwise exponent per unit 1/ln(N)
    # If z_m has log correction: z_m_eff(N) = z_m_inf + p/ln(N)
    print(f"\n{label}:")
    print(f"  Pairwise z_m:     {['%.4f' % z for z in pw_z]}")
    print(f"  Pairwise beta_me: {['%.4f' % b for b in pw_b]}")
    print(f"  Pairwise alpha:   {['%.4f' % a for a in pw_a]}")

    # z_m drift: total change from first to last pair
    z_drift = pw_z[-1] - pw_z[0]
    b_drift = pw_b[-1] - pw_b[0]
    a_drift = pw_a[-1] - pw_a[0]
    print(f"  z_m drift (last-first): {z_drift:+.4f}")
    print(f"  beta_me drift:          {b_drift:+.4f}")
    print(f"  alpha drift:            {a_drift:+.4f}")

    # Fit z_m_pw vs 1/ln(N_mid) to test log correction
    inv_ln_N = 1.0 / np.log(pw_N_mid)
    if len(pw_z) >= 3:
        p_z = np.polyfit(inv_ln_N, pw_z, 1)
        z_inf = p_z[1]  # extrapolated z_m at N→∞
        z_log_coeff = p_z[0]  # coefficient of 1/ln(N)
        print(f"  z_m = {z_inf:.4f} + {z_log_coeff:.3f}/ln(N)  (R²=", end="")
        ss_res = np.sum((pw_z - np.polyval(p_z, inv_ln_N))**2)
        ss_tot = np.sum((pw_z - np.mean(pw_z))**2)
        R2 = 1 - ss_res/ss_tot if ss_tot > 1e-15 else 0
        print(f"{R2:.4f})")

        p_b = np.polyfit(inv_ln_N, pw_b, 1)
        b_inf = p_b[1]
        b_log_coeff = p_b[0]
        print(f"  beta = {b_inf:.4f} + {b_log_coeff:.3f}/ln(N)  (R²=", end="")
        ss_res = np.sum((pw_b - np.polyval(p_b, inv_ln_N))**2)
        ss_tot = np.sum((pw_b - np.mean(pw_b))**2)
        R2 = 1 - ss_res/ss_tot if ss_tot > 1e-15 else 0
        print(f"{R2:.4f})")

        p_a = np.polyfit(inv_ln_N, pw_a, 1)
        a_inf = p_a[1]
        a_log_coeff = p_a[0]
        print(f"  alpha = {a_inf:.4f} + {a_log_coeff:.3f}/ln(N)  (R²=", end="")
        ss_res = np.sum((pw_a - np.polyval(p_a, inv_ln_N))**2)
        ss_tot = np.sum((pw_a - np.mean(pw_a))**2)
        R2 = 1 - ss_res/ss_tot if ss_tot > 1e-15 else 0
        print(f"{R2:.4f})")

        results['data'][label] = {
            'Ns': [int(n) for n in Ns],
            'gap_ms': [float(g) for g in gaps],
            'me_sqs': [float(m) for m in mes],
            'chi_Fs': [float(c) for c in chis],
            'pairwise_z': [float(z) for z in pw_z],
            'pairwise_beta': [float(b) for b in pw_b],
            'pairwise_alpha': [float(a) for a in pw_a],
            'z_inf': float(z_inf), 'z_log_coeff': float(z_log_coeff),
            'b_inf': float(b_inf), 'b_log_coeff': float(b_log_coeff),
            'a_inf': float(a_inf), 'a_log_coeff': float(a_log_coeff),
        }

save()

# Summary comparison table
print(f"\n{'=' * 70}")
print("SUMMARY: Log correction strength comparison")
print(f"{'':>5} {'z_m_inf':>8} {'z_log':>8} {'b_inf':>8} {'b_log':>8} {'a_inf':>8} {'a_log':>8}")
print("-" * 55)
for label in ['q=3', 'q=4', 'q=5']:
    d = results['data'].get(label, {})
    if d:
        print(f"{label:>5} {d['z_inf']:>8.4f} {d['z_log_coeff']:>8.3f} "
              f"{d['b_inf']:>8.4f} {d['b_log_coeff']:>8.3f} "
              f"{d['a_inf']:>8.4f} {d['a_log_coeff']:>8.3f}")

save()
print("\nResults saved.")
