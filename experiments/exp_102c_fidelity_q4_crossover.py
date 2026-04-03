#!/usr/bin/env python3
"""Sprint 102c: Fidelity susceptibility at q=4 (walking boundary) and q=5 extended.

q=4 is the real-to-complex CFT boundary. Does chi_F show the crossover?
Also: q=5 finer scan at n=6,8 for better peak characterization.
"""
import numpy as np
import json, time
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '102c_fidelity_q4_crossover',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_102c_fidelity_q4.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_parts(n, q):
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

dg = 1e-4

# ===== q=4 scan =====
print("=" * 60)
print("q=4: S_4 Potts fidelity susceptibility")
print("=" * 60)

g_c_4 = 0.25  # 1/q
q4_data = {}

for n in [6, 8, 10]:
    dim = 4**n
    print(f"\nq=4, n={n} (dim={dim})...", flush=True)
    t0 = time.time()

    H_coup, H_field = build_sq_potts_parts(n, 4)
    print(f"  Built in {time.time()-t0:.1f}s", flush=True)

    g_vals = np.linspace(g_c_4 - 0.08, g_c_4 + 0.08, 21)
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

    dt = time.time() - t0

    i_peak = np.argmax(chi_vals)
    g_peak = g_vals[i_peak]
    chi_peak = chi_vals[i_peak]
    above = np.array(chi_vals) > chi_peak / 2
    left, right = np.where(above)[0][0], np.where(above)[0][-1]
    fwhm = g_vals[right] - g_vals[left]

    q4_data[n] = {'chi_peak': chi_peak, 'g_peak': g_peak, 'fwhm': fwhm}
    print(f"  g_peak={g_peak:.4f}, chi_peak={chi_peak:.2f}, FWHM={fwhm:.4f}, time={dt:.1f}s")

    record(sprint=102, model='sq_potts', q=4, n=n, quantity='chi_F_max', value=chi_peak, method='fidelity_suscept')
    record(sprint=102, model='sq_potts', q=4, n=n, quantity='chi_F_FWHM', value=fwhm, method='fidelity_suscept')

    results['data'][f'q4_n{n}'] = {
        'g_peak': float(g_peak), 'chi_peak': float(chi_peak),
        'fwhm': float(fwhm), 'time': dt,
        'g_vals': g_vals.tolist(), 'chi_F': chi_vals,
    }
    save()

    if dt > 50:
        print(f"  Skipping larger sizes (took {dt:.0f}s)")
        break

# FSS for q=4
sizes_4 = sorted(q4_data.keys())
if len(sizes_4) >= 2:
    chi_arr = np.array([q4_data[n]['chi_peak'] for n in sizes_4])
    log_n = np.log(np.array(sizes_4, dtype=float))
    log_chi = np.log(chi_arr)
    if len(sizes_4) >= 3:
        p = np.polyfit(log_n, log_chi, 1)
        alpha = p[0]
    else:
        alpha = (log_chi[-1] - log_chi[0]) / (log_n[-1] - log_n[0])
    nu_4 = 2.0 / (alpha + 1)
    print(f"\nq=4 FSS: alpha={alpha:.3f}, nu={nu_4:.3f} (exact: 2/3=0.667)")
    record(sprint=102, model='sq_potts', q=4, n=max(sizes_4),
           quantity='nu_fidelity', value=nu_4, method='chi_F_peak_scaling')
    results['data']['q4_fss'] = {'alpha': float(alpha), 'nu': float(nu_4)}

# ===== Summary: nu(q) from fidelity =====
print("\n" + "=" * 60)
print("COMPLETE nu(q) FROM FIDELITY SUSCEPTIBILITY")
print(f"{'q':>3} {'nu_exact':>10} {'nu_fidelity':>12} {'ratio':>8}")
print("-" * 40)

nu_exact = {2: 1.0, 3: 5/6, 4: 2/3, 5: 0.83, 7: None}
nu_fid = {2: 1.009, 3: 0.841}
if len(sizes_4) >= 2:
    nu_fid[4] = nu_4
# q=5 from 102b
nu_fid[5] = 0.648

for q in [2, 3, 4, 5]:
    if q in nu_fid:
        ne = nu_exact.get(q)
        ne_str = f"{ne:.3f}" if ne else "?"
        ratio = nu_fid[q] / ne if ne else float('nan')
        print(f"{q:>3} {ne_str:>10} {nu_fid[q]:>12.3f} {ratio:>8.3f}")

# ===== chi_F(g_c, n=6) vs q =====
print("\nchi_F(g_c, n=6) exponential growth:")
chi_n6 = {2: 0.73, 3: 4.09, 4: q4_data[6]['chi_peak'], 5: 36.29, 7: 166.88}
for q, chi in sorted(chi_n6.items()):
    print(f"  q={q}: chi_F={chi:.2f}")

# Fit: chi_F(q) ~ exp(a*q + b)
qs = np.array(sorted(chi_n6.keys()), dtype=float)
chis = np.array([chi_n6[q] for q in sorted(chi_n6.keys())])
log_chis = np.log(chis)
p = np.polyfit(qs, log_chis, 1)
print(f"  Exponential fit: chi_F ~ exp({p[0]:.2f}*q + {p[1]:.2f})")
print(f"  Growth rate per unit q: {np.exp(p[0]):.2f}x")

results['data']['chi_vs_q'] = {
    'q_vals': qs.tolist(), 'chi_n6': chis.tolist(),
    'exp_rate': float(p[0]), 'growth_per_q': float(np.exp(p[0])),
}

save()
print("\nResults saved.")
