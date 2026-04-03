#!/usr/bin/env python3
"""Sprint 083a: Casimir energy extraction for q=2,3,5 on periodic S_q Potts chains.

CFT Casimir formula for periodic chain of N sites:
  E_0(N)/N = eps_inf - pi*v*c/(6*N^2) + O(1/N^4)

From multiple sizes, fit E_0/N vs 1/N^2 to get eps_inf (intercept) and v*c (slope).
Then v = v*c / c_known.

Also extract v from gap: v_gap = gap*N / (2*pi*x_sigma).

Sizes chosen to stay under ~60s per eigsh call:
  q=2: n=6,8,10,12,14 (max dim=16k)
  q=3: n=4,6,8,10 (max dim=59k)
  q=5: n=4,6,8 (max dim=390k)
"""
import numpy as np
import json, time
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh

results = {
    'experiment': '083a_casimir_q235',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_083a_casimir_q235.json", "w") as f:
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

        # Diagonal: coupling
        diag = 0.0
        for site in range(n):
            nxt = (site + 1) % n
            if c[site] == c[nxt]:
                diag -= 1.0
        H[idx, idx] += diag

        # Off-diagonal: transverse field (S_q: all cyclic shifts)
        for site in range(n):
            for k in range(1, q):
                c2 = c.copy()
                c2[site] = (c[site] + k) % q
                idx2 = 0
                for i in range(n - 1, -1, -1):
                    idx2 = idx2 * q + c2[i]
                H[idx, idx2] += -g

    return csr_matrix(H)

# Known values from Sprint 082 and exact CFT
ref_data = {
    2: {'c_exact': 0.500, 'x_sigma': 0.123, 'v_082': 1.02, 'gap_N': 0.786},
    3: {'c_exact': 0.800, 'x_sigma': 0.132, 'v_082': 0.89, 'gap_N': 0.733},
    5: {'c_exact': 1.138, 'x_sigma': 0.136, 'v_082': 0.75, 'gap_N': 0.639},
}

# Size lists per q
size_lists = {
    2: [6, 8, 10, 12, 14],
    3: [4, 6, 8, 10],
    5: [4, 6, 8],
}

print("Sprint 083a: Casimir energy extraction for S_q Potts (q=2,3,5)")
print("=" * 70, flush=True)

for q in [2, 3, 5]:
    gc = 1.0 / q
    sizes = size_lists[q]
    ref = ref_data[q]
    print(f"\n{'='*70}")
    print(f"q={q}, g_c=1/{q}={gc:.4f}, c_exact={ref['c_exact']}")
    print(f"{'='*70}", flush=True)

    q_data = {'q': q, 'gc': gc, 'sizes': [], 'E0_per_N': [], 'gaps': [], 'gap_N': []}

    for n in sizes:
        dim = q**n
        print(f"  n={n}, dim={dim:,} ... ", end='', flush=True)

        t0 = time.time()
        H = build_sq_potts_periodic(n, q, gc)
        t_build = time.time() - t0

        t0 = time.time()
        evals, evecs = eigsh(H, k=4, which='SA')
        t_eig = time.time() - t0

        order = np.argsort(evals)
        evals = evals[order]

        E0 = float(evals[0])
        E1 = float(evals[1])
        gap = E1 - E0
        E0_per_N = E0 / n
        gapN = gap * n

        q_data['sizes'].append(n)
        q_data['E0_per_N'].append(E0_per_N)
        q_data['gaps'].append(gap)
        q_data['gap_N'].append(gapN)

        print(f"E0/N={E0_per_N:.8f}, gap={gap:.6f}, gap*N={gapN:.4f}, "
              f"build={t_build:.1f}s, eig={t_eig:.1f}s", flush=True)

    # Fit Casimir formula: E0/N = eps_inf - pi*v*c/(6*N^2)
    # Linear regression: y = a + b*x where y=E0/N, x=1/N^2, b=-pi*v*c/6
    N_arr = np.array(q_data['sizes'], dtype=float)
    y = np.array(q_data['E0_per_N'])
    x = 1.0 / N_arr**2

    # Linear fit
    A = np.vstack([x, np.ones_like(x)]).T
    (slope, intercept), residuals, _, _ = np.linalg.lstsq(A, y, rcond=None)

    eps_inf = intercept
    vc_casimir = -slope * 6 / np.pi  # v*c = -slope * 6/pi

    # Predicted values and R^2
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    # Extract v
    c_known = ref['c_exact']
    v_casimir = vc_casimir / c_known

    # Also compute pairwise v*c from consecutive size pairs
    pairwise_vc = []
    for i in range(len(N_arr) - 1):
        N1, N2 = N_arr[i], N_arr[i+1]
        dE = y[i+1] - y[i]  # E0(N2)/N2 - E0(N1)/N1
        dx = 1/N2**2 - 1/N1**2
        vc_pair = -dE / dx * 6 / np.pi
        v_pair = vc_pair / c_known
        pairwise_vc.append({
            'pair': f'({int(N1)},{int(N2)})',
            'vc': float(vc_pair),
            'v': float(v_pair),
        })

    # v from gap (Sprint 082 method)
    # Use largest available size
    gap_N_best = q_data['gap_N'][-1]
    x_sigma = ref['x_sigma']
    v_gap = gap_N_best / (2 * np.pi * x_sigma)

    print(f"\n  Casimir fit: E0/N = {eps_inf:.8f} - {-slope:.6f}/N²")
    print(f"  R² = {R2:.8f}")
    print(f"  v*c (Casimir) = {vc_casimir:.5f}")
    print(f"  v (Casimir)   = {v_casimir:.4f}  (using c={c_known})")
    print(f"  v (gap/x_σ)   = {v_gap:.4f}  (gap×N={gap_N_best:.4f}, x_σ={x_sigma})")
    print(f"  v (Sprint 082)= {ref['v_082']:.4f}")
    print(f"  Agreement: Casimir vs gap = {100*abs(v_casimir-v_gap)/v_gap:+.1f}%")
    print(f"  Agreement: Casimir vs 082 = {100*(v_casimir-ref['v_082'])/ref['v_082']:+.1f}%")

    print(f"\n  Pairwise v*c from consecutive sizes:")
    for p in pairwise_vc:
        print(f"    {p['pair']}: v*c={p['vc']:.5f}, v={p['v']:.4f}")

    # Per-point residuals from Casimir fit
    print(f"\n  Per-point residuals from Casimir fit:")
    for i, n in enumerate(q_data['sizes']):
        resid = y[i] - y_pred[i]
        print(f"    n={n}: E0/N={y[i]:.8f}, pred={y_pred[i]:.8f}, resid={resid:+.2e}")

    q_data['fit'] = {
        'eps_inf': float(eps_inf),
        'slope': float(slope),
        'vc_casimir': float(vc_casimir),
        'v_casimir': float(v_casimir),
        'v_gap': float(v_gap),
        'v_082': ref['v_082'],
        'R2': float(R2),
        'c_used': c_known,
    }
    q_data['pairwise_vc'] = pairwise_vc

    results['data'][str(q)] = q_data
    save()

# Grand summary
print(f"\n{'='*70}")
print("GRAND SUMMARY: Casimir velocity vs correlator velocity")
print(f"{'='*70}")
print(f"{'q':>3} {'v(Casimir)':>11} {'v(gap/xσ)':>11} {'v(082)':>8} {'Cas/082':>8} {'R²':>10}")
for q in [2, 3, 5]:
    d = results['data'][str(q)]
    f = d['fit']
    ratio = f['v_casimir'] / f['v_082']
    print(f"{q:3d} {f['v_casimir']:11.5f} {f['v_gap']:11.5f} {f['v_082']:8.4f} {ratio:8.4f} {f['R2']:10.8f}")

save()
print(f"\nSaved to results/sprint_083a_casimir_q235.json")

from db_utils import record
for q in [2, 3, 5]:
    d = results['data'][str(q)]
    f = d['fit']
    record(sprint=83, model='sq_potts', q=q, n=max(d['sizes']),
           quantity='v_casimir', value=f['v_casimir'],
           method='casimir_energy',
           notes=f'vc={f["vc_casimir"]:.5f}, c={f["c_used"]}, R2={f["R2"]:.6f}')
    record(sprint=83, model='sq_potts', q=q, n=max(d['sizes']),
           quantity='eps_inf', value=f['eps_inf'],
           method='casimir_fit',
           notes=f'from {len(d["sizes"])} sizes, R2={f["R2"]:.6f}')
print("Recorded to DB.")
