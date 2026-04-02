#!/usr/bin/env python3
"""Sprint 079b: Exact diag c(q=7,10) for S_q Potts — no chi truncation.

DMRG at q=7 is chi-limited (chi=64 vs 2401 exact center-bond states).
Use GPU exact diag at accessible sizes to get TRUE entropy profile.

q=7: n=6 (dim=118k), n=8 (dim=5.8M) — both GPU-feasible
q=10: n=5 (dim=100k), n=6 (dim=1M) — GPU-feasible

Compare to DMRG results to quantify truncation error.
"""
import numpy as np
import json, time
from scipy.sparse import kron, eye, csr_matrix
from scipy.linalg import svd
from gpu_utils import eigsh

results = {
    'experiment': '079b_sq_exact_c_q7q10',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

def save():
    with open("results/sprint_079b_sq_exact_c_q7q10.json", "w") as f:
        json.dump(results, f, indent=2, default=str)


def build_sq_potts_H(n, q, g):
    """Build S_q Potts Hamiltonian: H = -J*Σ δ(s_i,s_j) - g*Σ Sq_field_i.

    Sq_field = ones(q,q) - eye(q) = Σ_{k=1}^{q-1} X^k (full S_q transverse field).
    """
    dim = q ** n
    I_q = eye(q, format='csr')

    # S_q field: ones(q,q) - I
    Sq_field = csr_matrix(np.ones((q, q)) - np.eye(q))

    H = csr_matrix((dim, dim), dtype=float)

    # Coupling: -J * Σ_<ij> δ(s_i, s_j) = -J * Σ_<ij> Σ_a P_a⊗P_a
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

    # Field: -g * Σ_i Sq_field_i
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
    """Compute entanglement entropy at each bond cut via SVD."""
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


def fit_cc_central(L, S_profile):
    """Fit Calabrese-Cardy formula, central half only."""
    x_vals = np.arange(1, L)
    chord = np.log((2 * L / np.pi) * np.sin(np.pi * x_vals / L))
    S_arr = np.array(S_profile)

    n_pts = len(x_vals)
    lo = n_pts // 4
    hi = 3 * n_pts // 4
    if hi - lo < 2:
        lo = 0
        hi = n_pts  # fall back to full range for small n

    chord_c = chord[lo:hi]
    S_c = S_arr[lo:hi]

    A = np.vstack([chord_c, np.ones_like(chord_c)]).T
    (slope, s0), _, _, _ = np.linalg.lstsq(A, S_c, rcond=None)
    c = 6 * slope

    S_pred = slope * chord_c + s0
    ss_res = np.sum((S_c - S_pred) ** 2)
    ss_tot = np.sum((S_c - np.mean(S_c)) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return float(c), float(s0), float(R2)


def complex_cft_c(q):
    """Correct complex CFT central charge.

    Standard parametrization: √Q = 2*cos(π/p), c = 1 - 6/[p*(p-1)].
    For Q > 4: p is complex. α = arccosh(√Q/2), p = iπ/α.
    Re(c) = 1 + 6α²/(π² + α²).
    """
    if q <= 4:
        p = np.pi / np.arccos(np.sqrt(q) / 2)
        c = 1 - 6 / (p * (p - 1))
        return float(c), 0.0
    else:
        alpha = np.arccosh(np.sqrt(q) / 2)
        Re_c = 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)
        Im_c = -6 * alpha**3 / (np.pi * (np.pi**2 + alpha**2))
        return float(Re_c), float(Im_c)


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 079b: Exact diag c for S_q Potts at q=7,10", flush=True)
print("=" * 60, flush=True)

# Complex CFT predictions (correct formula)
for q_test in [2, 3, 5, 7, 10]:
    Re_c, Im_c = complex_cft_c(q_test)
    print(f"  q={q_test:2d}: Re(c) = {Re_c:.4f}, Im(c) = {Im_c:.4f}", flush=True)
results['complex_cft'] = {}
for q_test in [5, 7, 10]:
    Re_c, Im_c = complex_cft_c(q_test)
    results['complex_cft'][f'q{q_test}'] = {'Re_c': Re_c, 'Im_c': Im_c}
save()

# q=7 exact diag
print(f"\n{'='*60}", flush=True)
print("q=7 S_q Potts, g_c = 1/7", flush=True)
print(f"{'='*60}", flush=True)

q7_data = {}
for n in [6, 8]:
    dim = 7 ** n
    gc = 1.0 / 7
    print(f"\n  n={n}, dim={dim:,}...", flush=True)
    t0 = time.time()
    H = build_sq_potts_H(n, 7, gc)
    t_build = time.time() - t0
    print(f"  H built in {t_build:.1f}s", flush=True)

    t0 = time.time()
    evals, evecs = eigsh(H, k=1, which='SA')
    t_eig = time.time() - t0
    E0 = evals[0]
    psi = evecs[:, 0]
    print(f"  eigsh in {t_eig:.1f}s, E0 = {E0:.8f}", flush=True)

    t0 = time.time()
    S = entropy_profile(psi, n, 7)
    t_ent = time.time() - t0
    c, s0, R2 = fit_cc_central(n, S)
    S_mid = S[n // 2 - 1]
    print(f"  Entropy in {t_ent:.1f}s", flush=True)
    print(f"  c = {c:.4f}, R² = {R2:.6f}, S_mid = {S_mid:.4f}", flush=True)
    print(f"  S profile: {[f'{s:.4f}' for s in S]}", flush=True)

    q7_data[f'n{n}'] = {
        'n': n, 'dim': dim, 'E0': float(E0),
        'c': c, 's0': s0, 'R2': R2, 'S_mid': S_mid,
        'S_profile': S,
        'time_build': t_build, 'time_eig': t_eig, 'time_ent': t_ent,
    }
    results['q7'] = q7_data
    save()

    total_t = t_build + t_eig + t_ent
    if total_t > 250:
        print(f"  Time limit reached.", flush=True)
        break

# q=10 exact diag
print(f"\n{'='*60}", flush=True)
print("q=10 S_q Potts, g_c = 1/10", flush=True)
print(f"{'='*60}", flush=True)

q10_data = {}
for n in [5, 6]:
    dim = 10 ** n
    gc = 0.1
    print(f"\n  n={n}, dim={dim:,}...", flush=True)
    t0 = time.time()
    H = build_sq_potts_H(n, 10, gc)
    t_build = time.time() - t0
    print(f"  H built in {t_build:.1f}s", flush=True)

    t0 = time.time()
    evals, evecs = eigsh(H, k=1, which='SA')
    t_eig = time.time() - t0
    E0 = evals[0]
    psi = evecs[:, 0]
    print(f"  eigsh in {t_eig:.1f}s, E0 = {E0:.8f}", flush=True)

    t0 = time.time()
    S = entropy_profile(psi, n, 10)
    t_ent = time.time() - t0
    c, s0, R2 = fit_cc_central(n, S)
    S_mid = S[n // 2 - 1]
    print(f"  Entropy in {t_ent:.1f}s", flush=True)
    print(f"  c = {c:.4f}, R² = {R2:.6f}, S_mid = {S_mid:.4f}", flush=True)
    print(f"  S profile: {[f'{s:.4f}' for s in S]}", flush=True)

    q10_data[f'n{n}'] = {
        'n': n, 'dim': dim, 'E0': float(E0),
        'c': c, 's0': s0, 'R2': R2, 'S_mid': S_mid,
        'S_profile': S,
        'time_build': t_build, 'time_eig': t_eig, 'time_ent': t_ent,
    }
    results['q10'] = q10_data
    save()

    total_t = t_build + t_eig + t_ent
    if total_t > 250:
        print(f"  Time limit reached.", flush=True)
        break

# Summary
print(f"\n{'='*60}", flush=True)
print("SUMMARY: Exact diag c_eff vs complex CFT Re(c)", flush=True)
print(f"{'='*60}", flush=True)

Re_c_7, _ = complex_cft_c(7)
Re_c_10, _ = complex_cft_c(10)
Re_c_5, _ = complex_cft_c(5)
print(f"\nComplex CFT predictions: q=5 Re(c)={Re_c_5:.4f}, q=7 Re(c)={Re_c_7:.4f}, q=10 Re(c)={Re_c_10:.4f}", flush=True)

print(f"\nq=7 (CFT Re(c)={Re_c_7:.4f}):", flush=True)
for k, v in q7_data.items():
    print(f"  {k}: c = {v['c']:.4f} (R² = {v['R2']:.6f}), ratio c/Re(c) = {v['c']/Re_c_7:.3f}", flush=True)

print(f"\nq=10 (CFT Re(c)={Re_c_10:.4f}):", flush=True)
for k, v in q10_data.items():
    print(f"  {k}: c = {v['c']:.4f} (R² = {v['R2']:.6f}), ratio c/Re(c) = {v['c']/Re_c_10:.3f}", flush=True)

# DMRG comparison for q=7
print(f"\nDMRG vs exact comparison (q=7):", flush=True)
print(f"  DMRG n=8 chi=56: c = 1.1108", flush=True)
if 'n8' in q7_data:
    print(f"  Exact n=8:        c = {q7_data['n8']['c']:.4f}", flush=True)
    print(f"  DMRG truncation error: {abs(1.1108 - q7_data['n8']['c']):.4f} ({abs(1.1108 - q7_data['n8']['c'])/q7_data['n8']['c']*100:.1f}%)", flush=True)

save()
print("\nSaved to results/sprint_079b_sq_exact_c_q7q10.json", flush=True)

from db_utils import record
for q_val, data in [(7, q7_data), (10, q10_data)]:
    for k, v in data.items():
        record(sprint=79, model='sq_potts', q=q_val, n=v['n'],
               quantity='c_eff', value=v['c'],
               method='exact_diag_CC', notes=f'gc=1/{q_val}, R2={v["R2"]:.4f}, exact (no chi truncation)')
print("Recorded to DB.", flush=True)
