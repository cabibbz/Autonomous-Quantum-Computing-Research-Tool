#!/usr/bin/env python3
"""Sprint 102b: FSS analysis of fidelity susceptibility peaks.

chi_F_max ~ N^(2/nu - 1) for continuous transition
FWHM ~ N^(-1/nu) for continuous, N^(-1) for first-order

Also: q=2 n=14 (dim=16k, fast) to extend scaling range.
And: finer g-scan near peak for better peak fitting.
"""
import numpy as np
import json, time
from scipy.sparse import lil_matrix, csr_matrix
from scipy.optimize import curve_fit
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '102b_fidelity_fss',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_102b_fidelity_fss.json", "w") as f:
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

def fidelity_scan(H_coup, H_field, n, g_vals, dg=1e-4):
    """Compute chi_F at each g value."""
    chi_vals = []
    for g in g_vals:
        H1 = H_coup + g * H_field
        H2 = H_coup + (g + dg) * H_field
        _, v1 = eigsh(H1, k=1, which='SA')
        _, v2 = eigsh(H2, k=1, which='SA')
        overlap = abs(np.vdot(v1[:, 0], v2[:, 0]))
        overlap = min(overlap, 1.0 - 1e-15)
        chi_F = 2.0 * (1.0 - overlap) / (dg**2 * n)
        chi_vals.append(float(chi_F))
    return chi_vals

# ===== Part 1: Additional data point q=2 n=14 =====
print("Part 1: q=2 n=14 scan (dim=16384)")
t0 = time.time()
n, q = 14, 2
g_c = 0.5
g_vals = np.linspace(0.35, 0.65, 25)
H_coup, H_field = build_sq_potts_parts(n, q)
chi_vals = fidelity_scan(H_coup, H_field, n, g_vals)
i_peak = np.argmax(chi_vals)
g_peak = g_vals[i_peak]
chi_peak = chi_vals[i_peak]
above = np.array(chi_vals) > chi_peak / 2
left, right = np.where(above)[0][0], np.where(above)[0][-1]
fwhm = g_vals[right] - g_vals[left]
dt = time.time() - t0
print(f"  g_peak={g_peak:.4f}, chi_peak={chi_peak:.2f}, FWHM={fwhm:.4f}, time={dt:.1f}s")
record(sprint=102, model='sq_potts', q=2, n=14, quantity='chi_F_max', value=chi_peak, method='fidelity_suscept')
record(sprint=102, model='sq_potts', q=2, n=14, quantity='chi_F_FWHM', value=fwhm, method='fidelity_suscept')

results['data']['q2_n14'] = {
    'g_peak': float(g_peak), 'chi_peak': float(chi_peak), 'fwhm': float(fwhm), 'time': dt
}
save()

# ===== Part 2: FSS analysis =====
print("\nPart 2: Finite-size scaling analysis")
print("=" * 60)

# Collected data: (q, sizes, chi_peaks, fwhms)
fss_data = {
    2: {
        'sizes': [6, 8, 10, 12, 14],
        'chi_peaks': [0.73, 0.95, 1.19, 1.43, chi_peak],
        'fwhms': [0.2400, 0.2250, 0.2100, 0.1800, fwhm],
        'nu_exact': 1.0,  # Ising
    },
    3: {
        'sizes': [6, 8, 10],
        'chi_peaks': [4.09, 6.12, 8.27],
        'fwhms': [0.1333, 0.0889, 0.0778],
        'nu_exact': 5/6,  # 3-state Potts
    },
    5: {
        'sizes': [6, 8],
        'chi_peaks': [36.29, 66.18],
        'fwhms': [0.0375, 0.0225],
        'nu_exact': 0.83,  # Walking estimate
    },
    7: {
        'sizes': [6],
        'chi_peaks': [166.88],
        'fwhms': [0.0171],
        'nu_exact': None,
    },
}

print(f"\n{'q':>3} {'nu_exact':>8} {'alpha_chi':>10} {'nu_chi':>8} {'alpha_w':>10} {'nu_w':>8}")
print("-" * 55)

for q, d in fss_data.items():
    sizes = np.array(d['sizes'], dtype=float)
    chi = np.array(d['chi_peaks'])
    fwhm = np.array(d['fwhms'])

    if len(sizes) >= 2:
        # chi_max ~ N^alpha where alpha = 2/nu - 1
        # log-log fit
        log_n = np.log(sizes)
        log_chi = np.log(chi)
        if len(sizes) >= 3:
            p_chi = np.polyfit(log_n, log_chi, 1)
            alpha_chi = p_chi[0]
        else:
            # Two points: slope
            alpha_chi = (log_chi[-1] - log_chi[0]) / (log_n[-1] - log_n[0])
        nu_chi = 2.0 / (alpha_chi + 1)

        # FWHM ~ N^(-1/nu) → log(FWHM) = -1/nu * log(N) + const
        log_fwhm = np.log(fwhm)
        if len(sizes) >= 3:
            p_w = np.polyfit(log_n, log_fwhm, 1)
            alpha_w = p_w[0]
        else:
            alpha_w = (log_fwhm[-1] - log_fwhm[0]) / (log_n[-1] - log_n[0])
        nu_w = -1.0 / alpha_w

        nu_str = f"{d['nu_exact']:.3f}" if d['nu_exact'] else "?"
        print(f"{q:>3} {nu_str:>8} {alpha_chi:>10.4f} {nu_chi:>8.3f} {alpha_w:>10.4f} {nu_w:>8.3f}")

        results['data'][f'fss_q{q}'] = {
            'alpha_chi': float(alpha_chi),
            'nu_from_chi': float(nu_chi),
            'alpha_fwhm': float(alpha_w),
            'nu_from_fwhm': float(nu_w),
            'nu_exact': d['nu_exact'],
        }

        record(sprint=102, model='sq_potts', q=q, n=max(d['sizes']),
               quantity='nu_fidelity', value=nu_chi, method='chi_F_peak_scaling')
        record(sprint=102, model='sq_potts', q=q, n=max(d['sizes']),
               quantity='nu_fwhm', value=nu_w, method='chi_F_FWHM_scaling')
    else:
        print(f"{q:>3} {'?':>8} {'—':>10} {'—':>8} {'—':>10} {'—':>8} (only 1 size)")

# ===== Part 3: Pairwise nu extraction =====
print("\nPairwise nu from consecutive sizes:")
print(f"{'q':>3} {'n1':>4} {'n2':>4} {'nu(chi)':>8} {'nu(fwhm)':>10}")
print("-" * 40)
for q, d in fss_data.items():
    sizes = d['sizes']
    chi = d['chi_peaks']
    fwhm = d['fwhms']
    for i in range(len(sizes) - 1):
        n1, n2 = sizes[i], sizes[i+1]
        # alpha = ln(chi2/chi1) / ln(n2/n1)
        a_chi = np.log(chi[i+1] / chi[i]) / np.log(n2 / n1)
        nu_pw = 2.0 / (a_chi + 1)
        a_w = np.log(fwhm[i+1] / fwhm[i]) / np.log(n2 / n1)
        nu_pw_w = -1.0 / a_w
        print(f"{q:>3} {n1:>4} {n2:>4} {nu_pw:>8.3f} {nu_pw_w:>10.3f}")

# ===== Part 4: chi_F(g_c) ratio between q values =====
print("\nchi_F at g_c, n=6 (all q):")
chi_n6 = {2: 0.73, 3: 4.09, 5: 36.29, 7: 166.88}
for q, chi in chi_n6.items():
    print(f"  q={q}: chi_F={chi:.2f}, ratio to q=2: {chi/0.73:.1f}x")

results['data']['chi_ratio_n6'] = {str(q): float(chi/0.73) for q, chi in chi_n6.items()}

save()
print("\nResults saved to results/sprint_102b_fidelity_fss.json")
