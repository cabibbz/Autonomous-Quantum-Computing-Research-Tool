#!/usr/bin/env python3
"""Sprint 080a: c_eff(q=6) from entropy profile at g_c = 1/6 (S_q Potts self-duality).

Exact diag at n=6,8 (dim=46k, 1.68M). GPU for n=8.
Also test n=10 (dim=60.5M) timing — likely infeasible but worth a single-point test.

Key question: Where does q=6 fall on the c_eff/Re(c) curve between q=5 (1.00) and q=7 (0.82)?
"""
import numpy as np
import json, time
from scipy.sparse import kron, eye, csr_matrix
from scipy.linalg import svd
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '080a_sq_potts_c_q6',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': 6, 'gc': 1/6,
}

def save():
    with open("results/sprint_080a_sq_potts_c_q6.json", "w") as f:
        json.dump(results, f, indent=2, default=str)


def build_sq_potts_H(n, q, g):
    """Build S_q Potts Hamiltonian: H = -sum delta(s_i,s_j) - g * sum (sum_k X^k)."""
    dim = q ** n
    I_q = eye(q, format='csr')
    # S_q field: sum_{k=1}^{q-1} X^k = J - I (all-ones minus identity)
    Sq_field = csr_matrix(np.ones((q, q)) - np.eye(q))
    H = csr_matrix((dim, dim), dtype=float)
    # Coupling: -delta(s_i, s_j) for nearest neighbors
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
    # Field: -g * sum_i Sq_field_i
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
    """Entanglement entropy for all bipartitions 1..n-1."""
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
    """Calabrese-Cardy fit: S = (c/6) ln[(2L/pi) sin(pi x/L)] + s0."""
    x_vals = np.arange(1, L)
    chord = np.log((2 * L / np.pi) * np.sin(np.pi * x_vals / L))
    S_arr = np.array(S_profile)
    # Use central half for fit
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
    """Complex CFT central charge from Coulomb gas."""
    if q <= 4:
        p = np.pi / np.arccos(np.sqrt(q) / 2)
        return 1 - 6 / (p * (p - 1)), 0.0
    else:
        alpha = np.arccosh(np.sqrt(q) / 2)
        Re_c = 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)
        Im_c = -6 * alpha**3 / (np.pi * (np.pi**2 + alpha**2))
        return float(Re_c), float(Im_c)


# ---- MAIN ----
q = 6
gc = 1.0 / q
Re_c, Im_c = complex_cft_c(q)
print(f"Sprint 080a: S_q Potts c_eff at q={q}, g_c=1/{q}={gc:.6f}")
print(f"Complex CFT: Re(c) = {Re_c:.4f}, Im(c) = {Im_c:.4f}")
print("=" * 60, flush=True)

results['Re_c'] = Re_c
results['Im_c'] = Im_c
results['sizes'] = {}

for n in [6, 8]:
    dim = q ** n
    print(f"\nn={n}, dim={dim:,}", flush=True)
    t0 = time.time()
    H = build_sq_potts_H(n, q, gc)
    t_build = time.time() - t0
    print(f"  H built [{t_build:.1f}s]", flush=True)

    t0 = time.time()
    evals, evecs = eigsh(H, k=4, which='SA')
    t_eig = time.time() - t0
    psi = evecs[:, 0]
    print(f"  eigsh k=4 [{t_eig:.1f}s]", flush=True)
    print(f"  E0 = {evals[0]:.8f}", flush=True)

    # Gaps
    gaps = evals[1:] - evals[0]
    gap_N = gaps[0] * n
    print(f"  gap = {gaps[0]:.6f}, gap*N = {gap_N:.4f}", flush=True)
    print(f"  gaps: {gaps}", flush=True)

    # Degeneracy check
    degen = [1]
    for i in range(1, len(gaps)):
        if abs(gaps[i] - gaps[i-1]) / max(abs(gaps[i]), 1e-10) < 0.01:
            degen[-1] += 1
        else:
            degen.append(1)
    print(f"  degeneracy pattern: {degen}", flush=True)

    # Entropy profile
    t0 = time.time()
    S = entropy_profile(psi, n, q)
    t_ent = time.time() - t0
    print(f"  entropy [{t_ent:.1f}s]: {[f'{s:.4f}' for s in S]}", flush=True)

    c, R2 = fit_cc(n, S)
    S_mid = S[n // 2 - 1]
    print(f"  c_eff = {c:.4f}, R² = {R2:.6f}, S_mid = {S_mid:.4f}", flush=True)

    ratio = c / Re_c
    print(f"  c_eff/Re(c) = {ratio:.4f}", flush=True)

    results['sizes'][str(n)] = {
        'dim': dim, 'c': c, 'R2': R2, 'S_mid': S_mid, 'S': S,
        'E0': float(evals[0]), 'gaps': [float(g) for g in gaps],
        'gap_N': float(gap_N), 'degen': degen,
        'time_build': t_build, 'time_eig': t_eig, 'time_ent': t_ent,
        'ratio_c_Re_c': ratio,
    }
    save()

    # Record to DB
    record(sprint=80, model='sq_potts', q=q, n=n,
           quantity='c_eff', value=c, method='exact_diag_CC',
           notes=f'gc=1/{q}, R2={R2:.4f}')
    record(sprint=80, model='sq_potts', q=q, n=n,
           quantity='gap_N', value=gap_N, method='exact_diag',
           notes=f'gc=1/{q}')

print(f"\n{'='*60}")
print(f"RESULT: q=6 c_eff/Re(c) = {results['sizes'].get('8', results['sizes'].get('6', {})).get('ratio_c_Re_c', 'N/A'):.4f}")
print(f"{'='*60}")
save()
print("Done. Saved to results/sprint_080a_sq_potts_c_q6.json")
