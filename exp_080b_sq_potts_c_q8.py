#!/usr/bin/env python3
"""Sprint 080b: c_eff(q=8) from entropy profile at g_c = 1/8 (S_q Potts self-duality).

Exact diag at n=5,6,7 (dim=32k, 262k, 2.1M). GPU for n=7.
Also try q=9 at n=6 (dim=531k) to get one more data point.

Key question: Where does q=8 fall? Between q=7 (0.82) and q=10 (0.60)?
"""
import numpy as np
import json, time
from scipy.sparse import kron, eye, csr_matrix
from scipy.linalg import svd
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '080b_sq_potts_c_q8_q9',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

def save():
    with open("results/sprint_080b_sq_potts_c_q8.json", "w") as f:
        json.dump(results, f, indent=2, default=str)


def build_sq_potts_H(n, q, g):
    """Build S_q Potts Hamiltonian."""
    dim = q ** n
    I_q = eye(q, format='csr')
    Sq_field = csr_matrix(np.ones((q, q)) - np.eye(q))
    H = csr_matrix((dim, dim), dtype=float)
    for site in range(n - 1):
        for a in range(q):
            Pa = csr_matrix(np.diag([1.0 if s == a else 0.0 for s in range(q)]))
            op = csr_matrix(eye(1))
            for j in range(n):
                if j == site or j == site + 1:
                    op = kron(op, Pa, format='csr')
                else:
                    op = kron(op, I_q, format='csr')
            H -= op
    for site in range(n):
        op = csr_matrix(eye(1))
        for j in range(n):
            if j == site:
                op = kron(op, Sq_field, format='csr')
            else:
                op = kron(op, I_q, format='csr')
        H -= g * op
    return H


def entropy_profile(psi, n, q):
    S_list = []
    psi_vec = psi.reshape(-1)
    for cut in range(1, n):
        d_left = q ** cut
        d_right = q ** (n - cut)
        psi_mat = psi_vec.reshape(d_left, d_right)
        sv = svd(psi_mat, compute_uv=False)
        sv = sv[sv > 1e-15]
        sv2 = sv ** 2
        sv2 /= sv2.sum()
        S = -np.sum(sv2 * np.log(sv2))
        S_list.append(float(S))
    return S_list


def fit_cc(L, S_profile):
    x_vals = np.arange(1, L)
    chord = np.log((2 * L / np.pi) * np.sin(np.pi * x_vals / L))
    S_arr = np.array(S_profile)
    n_pts = len(x_vals)
    lo = n_pts // 4
    hi = 3 * n_pts // 4
    if hi - lo < 2:
        lo = 0; hi = n_pts
    A = np.vstack([chord[lo:hi], np.ones(hi - lo)]).T
    (slope, s0), _, _, _ = np.linalg.lstsq(A, S_arr[lo:hi], rcond=None)
    c = 6 * slope
    S_pred = slope * chord[lo:hi] + s0
    ss_res = np.sum((S_arr[lo:hi] - S_pred) ** 2)
    ss_tot = np.sum((S_arr[lo:hi] - np.mean(S_arr[lo:hi])) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    return float(c), float(R2)


def complex_cft_c(q):
    if q <= 4:
        p = np.pi / np.arccos(np.sqrt(q) / 2)
        return 1 - 6 / (p * (p - 1)), 0.0
    else:
        alpha = np.arccosh(np.sqrt(q) / 2)
        Re_c = 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)
        Im_c = -6 * alpha**3 / (np.pi * (np.pi**2 + alpha**2))
        return float(Re_c), float(Im_c)


# ---- MAIN ----
print("Sprint 080b: S_q Potts c_eff at q=8 and q=9")
print("=" * 60, flush=True)

for q in [8, 9]:
    gc = 1.0 / q
    Re_c, Im_c = complex_cft_c(q)
    print(f"\n--- q={q}, g_c=1/{q}={gc:.6f}, Re(c)={Re_c:.4f} ---", flush=True)

    results[f'q{q}'] = {'gc': gc, 'Re_c': Re_c, 'Im_c': Im_c, 'sizes': {}}

    # Determine sizes
    if q == 8:
        sizes = [5, 6, 7]
    else:  # q=9
        sizes = [5, 6]

    for n in sizes:
        dim = q ** n
        if dim > 10_000_000:
            print(f"  n={n} dim={dim:,} — skipping (>10M)", flush=True)
            continue

        print(f"  n={n}, dim={dim:,}...", end=" ", flush=True)
        t0 = time.time()
        H = build_sq_potts_H(n, q, gc)
        t_build = time.time() - t0

        t0 = time.time()
        evals, evecs = eigsh(H, k=4, which='SA')
        t_eig = time.time() - t0
        psi = evecs[:, 0]

        gaps = evals[1:] - evals[0]
        gap_N = gaps[0] * n

        S = entropy_profile(psi, n, q)
        c, R2 = fit_cc(n, S)
        S_mid = S[n // 2 - 1]
        ratio = c / Re_c
        dt = t_build + t_eig

        print(f"c={c:.4f}, R²={R2:.6f}, c/Re(c)={ratio:.4f}, gap*N={gap_N:.4f} [{dt:.1f}s]", flush=True)
        print(f"    gaps: {[f'{g:.6f}' for g in gaps]}", flush=True)
        print(f"    S: {[f'{s:.4f}' for s in S]}", flush=True)

        results[f'q{q}']['sizes'][str(n)] = {
            'dim': dim, 'c': c, 'R2': R2, 'S_mid': S_mid, 'S': S,
            'E0': float(evals[0]), 'gaps': [float(g) for g in gaps],
            'gap_N': float(gap_N), 'ratio_c_Re_c': ratio,
            'time': dt,
        }
        save()

        record(sprint=80, model='sq_potts', q=q, n=n,
               quantity='c_eff', value=c, method='exact_diag_CC',
               notes=f'gc=1/{q}, R2={R2:.4f}')
        record(sprint=80, model='sq_potts', q=q, n=n,
               quantity='gap_N', value=gap_N, method='exact_diag',
               notes=f'gc=1/{q}')

# Summary
print(f"\n{'='*60}")
print("SUMMARY: c_eff/Re(c) at n=best")
print(f"{'='*60}")
for q in [8, 9]:
    d = results[f'q{q}']
    sizes = d['sizes']
    best_n = max(int(k) for k in sizes.keys())
    r = sizes[str(best_n)]
    print(f"  q={q}: c_eff={r['c']:.4f}, Re(c)={d['Re_c']:.4f}, ratio={r['ratio_c_Re_c']:.4f} (n={best_n})")

save()
print("\nDone. Saved to results/sprint_080b_sq_potts_c_q8.json")
