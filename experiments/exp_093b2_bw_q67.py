#!/usr/bin/env python3
"""Sprint 093b extension: Global BW fit for q=6,7 at nA=3 (n=6).
Also re-run q=2-5 at nA=3 for consistent comparison.
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record
import json, time

t0 = time.time()

results = {
    'experiment': '093b2_bw_q67',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_093b2_bw_q67.json", "w") as f:
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


nA = 3

for q in [2, 3, 4, 5, 6, 7]:
    n = 2 * nA
    g = 1.0 / q
    dimA = q**nA
    key = f"q{q}"
    print(f"\n=== q={q}, n={n}, nA={nA}, g_c={g:.4f}, dim={q**n} ===")
    t1 = time.time()

    H = build_sq_potts_periodic(n, q, g)
    evals, evecs = eigsh(H, k=1, which='SA')
    psi = evecs[:, 0]
    E0 = float(evals[0])

    rho_A = get_rho_A(psi, n, q, nA)
    H_E = entanglement_hamiltonian(rho_A)

    H_BW = build_H_BW(q, n, nA, g)
    R2_1p, alpha_1p = frobenius_R2(H_E, H_BW)

    # Also extract bond/field coefficients for profile
    bond_coeffs = []
    for x in range(nA - 1):
        delta = potts_delta_bond(q, nA, x, x + 1)
        c = np.sum(H_E * delta) / np.sum(delta * delta)
        bond_coeffs.append(float(c))

    field_coeffs = []
    for x in range(nA):
        field = clock_field_site(q, nA, x)
        c = np.sum(H_E * field) / np.sum(field * field)
        field_coeffs.append(float(c))

    print(f"  E0 = {E0:.6f}")
    print(f"  1-param BW: alpha={alpha_1p:.6f}, R²={R2_1p:.8f}")
    print(f"  1-R² = {1-R2_1p:.2e}")
    print(f"  Bond coeffs: {[f'{c:.4f}' for c in bond_coeffs]}")
    print(f"  Field coeffs: {[f'{c:.4f}' for c in field_coeffs]}")

    entry = {
        'q': q, 'n': n, 'nA': nA, 'g_c': g, 'E0': E0,
        'R2_1param': float(R2_1p),
        'alpha_1param': float(alpha_1p),
        'one_minus_R2': float(1 - R2_1p),
        'bond_coeffs': bond_coeffs,
        'field_coeffs': field_coeffs,
        'time': time.time() - t1,
    }
    results['data'][key] = entry

    record(sprint=93, model='sq_potts', q=q, n=n,
           quantity='bw_R2_nA3', value=R2_1p, method='global_bw_fit_nA3')

    save()
    print(f"  Time: {time.time()-t1:.1f}s")

print(f"\nTotal time: {time.time()-t0:.1f}s")
save()
print("Done!")
