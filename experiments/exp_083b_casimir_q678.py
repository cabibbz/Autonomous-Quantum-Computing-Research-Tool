#!/usr/bin/env python3
"""Sprint 083b: Casimir energy extraction for q=4,6,7,8 on periodic S_q Potts chains.

Extends 083a to cover q=4 (marginal Ashkin-Teller) and the walking-broken regime q=6,7,8.
Key question: Does the Casimir formula show larger deviations for q>5 (walking breakdown)?

Sizes:
  q=4: n=4,6,8 (max dim=65k)
  q=6: n=4,5,6,7 (max dim=280k)
  q=7: n=4,5,6,7 (max dim=824k)
  q=8: n=4,5,6 (max dim=262k)
"""
import numpy as np
import json, time
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh

results = {
    'experiment': '083b_casimir_q678',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_083b_casimir_q678.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_periodic(n, q, g):
    """Build S_q Potts Hamiltonian on periodic chain."""
    dim = q**n
    H = lil_matrix((dim, dim), dtype=complex)

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
        H[idx, idx] += diag

        for site in range(n):
            for k in range(1, q):
                c2 = c.copy()
                c2[site] = (c[site] + k) % q
                idx2 = 0
                for i in range(n - 1, -1, -1):
                    idx2 = idx2 * q + c2[i]
                H[idx, idx2] += -g

    return csr_matrix(H)

# Reference data from Sprint 082 and complex CFT
# c values: q=4 exact 1.0, q=6-8 use c_eff from KNOWLEDGE (best DMRG)
# For Casimir extraction, we use c to convert vc -> v
ref_data = {
    4: {'c': 1.000, 'x_sigma': 0.135, 'v_082': 0.81, 'gap_N': 0.686},
    6: {'c': 1.253, 'x_sigma': 0.135, 'v_082': 0.71, 'gap_N': 0.601},  # Re(c) for reference
    7: {'c': 1.351, 'x_sigma': 0.132, 'v_082': 0.68, 'gap_N': 0.560},
    8: {'c': 1.438, 'x_sigma': 0.130, 'v_082': 0.66, 'gap_N': 0.536},
}

# Also try with c_eff (measured) instead of Re(c):
c_eff_data = {
    4: 1.000,  # exact
    6: 1.147,  # Sprint 080a n=8
    7: 1.111,  # Sprint 079a n=8
    8: 1.062,  # Sprint 080b n=7
}

size_lists = {
    4: [4, 6, 8],
    6: [4, 5, 6, 7],
    7: [4, 5, 6, 7],
    8: [4, 5, 6],
}

print("Sprint 083b: Casimir energy extraction for S_q Potts (q=4,6,7,8)")
print("=" * 70, flush=True)

for q in [4, 6, 7, 8]:
    gc = 1.0 / q
    sizes = size_lists[q]
    ref = ref_data[q]
    print(f"\n{'='*70}")
    print(f"q={q}, g_c=1/{q}={gc:.4f}, Re(c)={ref['c']}, c_eff={c_eff_data[q]}")
    print(f"{'='*70}", flush=True)

    q_data = {'q': q, 'gc': gc, 'sizes': [], 'E0_per_N': [], 'gaps': [], 'gap_N': []}

    for n in sizes:
        dim = q**n
        print(f"  n={n}, dim={dim:,} ... ", end='', flush=True)

        t0 = time.time()
        H = build_sq_potts_periodic(n, q, gc)
        t_build = time.time() - t0

        t0 = time.time()
        evals, _ = eigsh(H, k=4, which='SA')
        t_eig = time.time() - t0

        evals = np.sort(evals)
        E0 = float(evals[0])
        gap = float(evals[1] - evals[0])
        E0_per_N = E0 / n
        gapN = gap * n

        q_data['sizes'].append(n)
        q_data['E0_per_N'].append(E0_per_N)
        q_data['gaps'].append(gap)
        q_data['gap_N'].append(gapN)

        print(f"E0/N={E0_per_N:.8f}, gap={gap:.6f}, gap*N={gapN:.4f}, "
              f"build={t_build:.1f}s, eig={t_eig:.1f}s", flush=True)

    # Fit Casimir formula
    N_arr = np.array(q_data['sizes'], dtype=float)
    y = np.array(q_data['E0_per_N'])
    x = 1.0 / N_arr**2

    A = np.vstack([x, np.ones_like(x)]).T
    (slope, intercept), _, _, _ = np.linalg.lstsq(A, y, rcond=None)

    eps_inf = intercept
    vc_casimir = -slope * 6 / np.pi

    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    # v using Re(c) and c_eff
    v_rec = vc_casimir / ref['c']
    v_ceff = vc_casimir / c_eff_data[q]

    # v from gap
    gap_N_best = q_data['gap_N'][-1]
    v_gap = gap_N_best / (2 * np.pi * ref['x_sigma'])

    # Pairwise
    pairwise = []
    for i in range(len(N_arr) - 1):
        N1, N2 = N_arr[i], N_arr[i+1]
        dE = y[i+1] - y[i]
        dx = 1/N2**2 - 1/N1**2
        vc_pair = -dE / dx * 6 / np.pi
        pairwise.append({
            'pair': f'({int(N1)},{int(N2)})',
            'vc': float(vc_pair),
            'v_rec': float(vc_pair / ref['c']),
            'v_ceff': float(vc_pair / c_eff_data[q]),
        })

    print(f"\n  Casimir fit: E0/N = {eps_inf:.8f} - {-slope:.6f}/N²")
    print(f"  R² = {R2:.8f}")
    print(f"  v*c (Casimir) = {vc_casimir:.5f}")
    print(f"  v (Re(c))     = {v_rec:.4f}")
    print(f"  v (c_eff)     = {v_ceff:.4f}")
    print(f"  v (gap/x_σ)   = {v_gap:.4f}  (gap×N={gap_N_best:.4f})")
    print(f"  v (Sprint 082)= {ref['v_082']:.4f}")

    print(f"\n  Pairwise v*c:")
    for p in pairwise:
        print(f"    {p['pair']}: v*c={p['vc']:.5f}, v(Re(c))={p['v_rec']:.4f}, v(c_eff)={p['v_ceff']:.4f}")

    print(f"\n  Residuals:")
    for i, n in enumerate(q_data['sizes']):
        resid = y[i] - y_pred[i]
        print(f"    n={n}: resid={resid:+.2e}")

    q_data['fit'] = {
        'eps_inf': float(eps_inf),
        'slope': float(slope),
        'vc_casimir': float(vc_casimir),
        'v_rec': float(v_rec),
        'v_ceff': float(v_ceff),
        'v_gap': float(v_gap),
        'v_082': ref['v_082'],
        'R2': float(R2),
        'c_rec': ref['c'],
        'c_eff': c_eff_data[q],
    }
    q_data['pairwise'] = pairwise

    results['data'][str(q)] = q_data
    save()

# Grand summary
print(f"\n{'='*70}")
print("GRAND SUMMARY: Casimir velocity for q=4,6,7,8")
print(f"{'='*70}")
print(f"{'q':>3} {'vc':>8} {'v(Re c)':>8} {'v(c_eff)':>9} {'v(gap)':>7} {'v(082)':>7} {'R²':>10}")
for q in [4, 6, 7, 8]:
    d = results['data'][str(q)]
    f = d['fit']
    print(f"{q:3d} {f['vc_casimir']:8.4f} {f['v_rec']:8.4f} {f['v_ceff']:9.4f} "
          f"{f['v_gap']:7.4f} {f['v_082']:7.4f} {f['R2']:10.6f}")

save()
print(f"\nSaved to results/sprint_083b_casimir_q678.json")

from db_utils import record
for q in [4, 6, 7, 8]:
    d = results['data'][str(q)]
    f = d['fit']
    record(sprint=83, model='sq_potts', q=q, n=max(d['sizes']),
           quantity='v_casimir', value=f['v_ceff'],
           method='casimir_energy',
           notes=f'vc={f["vc_casimir"]:.5f}, c_eff={f["c_eff"]}, R2={f["R2"]:.6f}')
    record(sprint=83, model='sq_potts', q=q, n=max(d['sizes']),
           quantity='eps_inf', value=f['eps_inf'],
           method='casimir_fit',
           notes=f'from {len(d["sizes"])} sizes, R2={f["R2"]:.6f}')
print("Recorded to DB.")
