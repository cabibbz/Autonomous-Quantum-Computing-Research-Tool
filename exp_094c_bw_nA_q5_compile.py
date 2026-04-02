#!/usr/bin/env python3
"""Sprint 094c: BW R² vs nA for q=5 + compile all q, fit gamma(q).

q=5: nA=3,4,5 (n=6,8,10, max dim=9765625 — GPU edge)
Then compile q=2-5 data from 094a/094b/094c and fit power-law gamma(q).
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.optimize import curve_fit
from gpu_utils import eigsh
from db_utils import record
import json, time

t0 = time.time()

results = {
    'experiment': '094c_bw_nA_q5_compile',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
    'compilation': {},
}

def save():
    with open("results/sprint_094c_bw_nA_q5_compile.json", "w") as f:
        json.dump(results, f, indent=2, default=str)


def build_sq_potts_periodic(n, q, g):
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


def get_rho_A(psi, n, q, nA):
    dimA = q**nA
    dimB = q**(n - nA)
    psi_mat = psi.reshape(dimA, dimB)
    rho_A = psi_mat @ psi_mat.conj().T
    return np.array(rho_A)


def entanglement_hamiltonian(rho_A):
    evals, evecs = np.linalg.eigh(rho_A)
    evals = np.clip(evals, 1e-30, None)
    H_E = -evecs @ np.diag(np.log(evals)) @ evecs.conj().T
    return np.real(H_E)


def potts_delta_bond(q, nA, si, sj):
    dimA = q**nA
    op = np.zeros((dimA, dimA))
    for idx in range(dimA):
        c = []
        tmp = idx
        for _ in range(nA):
            c.append(tmp % q)
            tmp //= q
        if c[si] == c[sj]:
            op[idx, idx] = 1.0
    return op


def clock_field_site(q, nA, site):
    dimA = q**nA
    op = np.zeros((dimA, dimA))
    for idx in range(dimA):
        c = []
        tmp = idx
        for _ in range(nA):
            c.append(tmp % q)
            tmp //= q
        for k in range(1, q):
            c2 = c.copy()
            c2[site] = (c[site] + k) % q
            idx2 = 0
            for i in range(nA - 1, -1, -1):
                idx2 = idx2 * q + c2[i]
            op[idx, idx2] += 1.0
    return op


def bw_envelope_bond(N, nA, bond_x):
    return (N / np.pi) * np.sin(np.pi * (bond_x + 1) / N) * \
           np.sin(np.pi * (nA - bond_x - 1) / N) / np.sin(np.pi * nA / N)


def bw_envelope_site(N, nA, site_x):
    return (N / np.pi) * np.sin(np.pi * (site_x + 0.5) / N) * \
           np.sin(np.pi * (nA - site_x - 0.5) / N) / np.sin(np.pi * nA / N)


def build_H_BW(q, n, nA, g):
    dimA = q**nA
    H_BW = np.zeros((dimA, dimA))
    for x in range(nA - 1):
        beta = bw_envelope_bond(n, nA, x)
        delta = potts_delta_bond(q, nA, x, x + 1)
        H_BW -= beta * delta
    for x in range(nA):
        beta = bw_envelope_site(n, nA, x)
        field = clock_field_site(q, nA, x)
        H_BW -= g * beta * field
    return H_BW


def frobenius_R2(H_target, H_model):
    dim = H_target.shape[0]
    t_mean = np.trace(H_target) / dim
    m_mean = np.trace(H_model) / dim
    H_t = H_target - t_mean * np.eye(dim)
    H_m = H_model - m_mean * np.eye(dim)
    alpha = np.sum(H_t * H_m) / np.sum(H_m * H_m)
    residual = H_t - alpha * H_m
    ss_res = np.sum(residual**2)
    ss_tot = np.sum(H_t**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return R2, alpha


# ============ q=5 ============
q = 5
g = 1.0 / q
nA_list = [3, 4]  # nA=5 means n=10, dim=5^10=9765625 — try separately

for nA in nA_list:
    n = 2 * nA
    dim = q**n
    dimA = q**nA
    key = f"q5_nA{nA}"
    print(f"\n=== q={q}, nA={nA}, n={n}, dim={dim}, dimA={dimA} ===")
    t1 = time.time()

    H = build_sq_potts_periodic(n, q, g)
    print(f"  H built, {time.time()-t1:.1f}s")
    evals, evecs = eigsh(H, k=1, which='SA')
    psi = evecs[:, 0]
    E0 = float(evals[0])
    print(f"  E0 = {E0:.6f}, diag time = {time.time()-t1:.1f}s")

    rho_A = get_rho_A(psi, n, q, nA)
    H_E = entanglement_hamiltonian(rho_A)
    H_BW = build_H_BW(q, n, nA, g)
    R2, alpha = frobenius_R2(H_E, H_BW)
    one_minus_R2 = 1.0 - R2
    print(f"  BW R² = {R2:.8f}, 1-R² = {one_minus_R2:.2e}, alpha = {alpha:.4f}")

    entry = {
        'q': q, 'nA': nA, 'n': n, 'dim': dim, 'dimA': dimA,
        'E0': E0, 'R2': float(R2), 'one_minus_R2': float(one_minus_R2),
        'alpha': float(alpha), 'time': time.time() - t1,
    }
    results['data'][key] = entry
    record(sprint=94, model='sq_potts', q=q, n=n,
           quantity='bw_R2_nA_scaling', value=R2, method='094c_nA_scan')
    save()

# Try nA=5 (dim=9.8M, GPU)
nA = 5
n = 2 * nA
dim = q**n
dimA = q**nA
key = f"q5_nA{nA}"
print(f"\n=== q={q}, nA={nA}, n={n}, dim={dim}, dimA={dimA} (GPU) ===")
t1 = time.time()
try:
    H = build_sq_potts_periodic(n, q, g)
    print(f"  H built, {time.time()-t1:.1f}s")
    evals, evecs = eigsh(H, k=1, which='SA')
    psi = evecs[:, 0]
    E0 = float(evals[0])
    print(f"  E0 = {E0:.6f}, diag time = {time.time()-t1:.1f}s")

    rho_A = get_rho_A(psi, n, q, nA)
    H_E = entanglement_hamiltonian(rho_A)

    # For nA=5, dimA=3125 — build_H_BW creates 3125x3125 dense matrices, feasible
    H_BW = build_H_BW(q, n, nA, g)
    R2, alpha = frobenius_R2(H_E, H_BW)
    one_minus_R2 = 1.0 - R2
    print(f"  BW R² = {R2:.8f}, 1-R² = {one_minus_R2:.2e}, alpha = {alpha:.4f}")

    entry = {
        'q': q, 'nA': nA, 'n': n, 'dim': dim, 'dimA': dimA,
        'E0': E0, 'R2': float(R2), 'one_minus_R2': float(one_minus_R2),
        'alpha': float(alpha), 'time': time.time() - t1,
    }
    results['data'][key] = entry
    record(sprint=94, model='sq_potts', q=q, n=n,
           quantity='bw_R2_nA_scaling', value=R2, method='094c_nA_scan')
    save()
except Exception as e:
    print(f"  q=5 nA=5 FAILED: {e}")
    results['data'][key] = {'error': str(e), 'time': time.time() - t1}
    save()

# ============ COMPILATION ============
print("\n\n========== COMPILATION ==========")

# Load results from all three experiments
all_data = {}
for fname in ['results/sprint_094a_bw_nA_q2.json',
              'results/sprint_094b_bw_nA_q34.json',
              'results/sprint_094c_bw_nA_q5_compile.json']:
    with open(fname) as f:
        d = json.load(f)
    for k, v in d['data'].items():
        if 'R2' in v:
            all_data[k] = v

# Organize by q
by_q = {}
for k, v in all_data.items():
    q_val = v.get('q', 2)  # 094a didn't store q explicitly
    if 'nA' not in k and k.startswith('nA'):
        q_val = 2
    elif k.startswith('q3'):
        q_val = 3
    elif k.startswith('q4'):
        q_val = 4
    elif k.startswith('q5'):
        q_val = 5
    if q_val not in by_q:
        by_q[q_val] = []
    by_q[q_val].append((v['nA'], v['one_minus_R2'], v['alpha']))

print("\nAll BW data (1-R² vs nA):")
print(f"{'q':>3} {'nA':>4} {'1-R²':>12} {'alpha':>8}")
print("-" * 35)
for q_val in sorted(by_q.keys()):
    pts = sorted(by_q[q_val])
    for nA, omr, alpha in pts:
        print(f"{q_val:>3} {nA:>4} {omr:>12.2e} {alpha:>8.3f}")

# Fit power law for each q: 1-R² = A * nA^gamma
print("\n\nPower-law fits: 1-R² = A * nA^gamma")
print(f"{'q':>3} {'gamma':>8} {'A':>12} {'fit_R2':>8} {'pts':>4}")
print("-" * 45)

compilation = {}
for q_val in sorted(by_q.keys()):
    pts = sorted(by_q[q_val])
    nA_arr = np.array([p[0] for p in pts])
    omr_arr = np.array([p[1] for p in pts])

    if len(pts) < 3:
        # Log-log linear fit
        log_nA = np.log(nA_arr)
        log_omr = np.log(omr_arr)
        coeffs = np.polyfit(log_nA, log_omr, 1)
        gamma = coeffs[0]
        A = np.exp(coeffs[1])
        pred = A * nA_arr**gamma
        ss_res = np.sum((omr_arr - pred)**2)
        ss_tot = np.sum((omr_arr - np.mean(omr_arr))**2)
        fit_R2 = 1 - ss_res / ss_tot if ss_tot > 1e-30 else 1.0
        method = 'log-log'
    else:
        # Curve fit
        log_nA = np.log(nA_arr)
        log_omr = np.log(omr_arr)
        coeffs = np.polyfit(log_nA, log_omr, 1)
        gamma = coeffs[0]
        A = np.exp(coeffs[1])
        pred = A * nA_arr**gamma
        ss_res = np.sum((omr_arr - pred)**2)
        ss_tot = np.sum((omr_arr - np.mean(omr_arr))**2)
        fit_R2 = 1 - ss_res / ss_tot if ss_tot > 1e-30 else 1.0
        method = 'log-log'

    print(f"{q_val:>3} {gamma:>8.2f} {A:>12.2e} {fit_R2:>8.4f} {len(pts):>4}")
    compilation[f'q{q_val}'] = {
        'gamma': float(gamma), 'A': float(A), 'fit_R2': float(fit_R2),
        'n_points': len(pts), 'method': method,
        'nA_values': [int(x) for x in nA_arr],
        'one_minus_R2': [float(x) for x in omr_arr],
    }

    record(sprint=94, model='sq_potts', q=q_val, n=0,
           quantity='bw_gamma_nA', value=gamma, method='094c_powerlaw_fit')

results['compilation'] = compilation

# Gamma vs q
print("\n\nGamma(q) summary:")
q_vals = sorted(compilation.keys(), key=lambda x: int(x[1:]))
for qk in q_vals:
    c = compilation[qk]
    print(f"  {qk}: gamma = {c['gamma']:.2f}, A = {c['A']:.2e}")

# Check if gamma increases with q
gammas = [compilation[f'q{q_val}']['gamma'] for q_val in sorted(by_q.keys())]
q_list = sorted(by_q.keys())
print(f"\n  q = {q_list}")
print(f"  gamma = {[f'{g:.2f}' for g in gammas]}")

# Does gamma increase monotonically?
monotonic = all(gammas[i] <= gammas[i+1] for i in range(len(gammas)-1))
print(f"  Monotonically increasing: {monotonic}")

# Fit gamma(q)
if len(q_list) >= 3:
    q_arr = np.array(q_list, dtype=float)
    g_arr = np.array(gammas)
    coeffs_gq = np.polyfit(q_arr, g_arr, 1)
    print(f"  Linear fit: gamma(q) = {coeffs_gq[0]:.2f}*q + {coeffs_gq[1]:.2f}")
    results['gamma_vs_q'] = {
        'q_values': [int(x) for x in q_list],
        'gamma_values': [float(x) for x in gammas],
        'linear_slope': float(coeffs_gq[0]),
        'linear_intercept': float(coeffs_gq[1]),
    }

save()
print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done!")
