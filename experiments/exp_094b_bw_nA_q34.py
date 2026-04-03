#!/usr/bin/env python3
"""Sprint 094b: BW R² vs subsystem size nA for q=3 and q=4.

q=3: nA=3,4,5,6 (n=6,8,10,12, max dim=531441)
q=4: nA=3,4,5 (n=6,8,10, max dim=1048576 — GPU)
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.optimize import curve_fit
from gpu_utils import eigsh
from db_utils import record
import json, time

t0 = time.time()

results = {
    'experiment': '094b_bw_nA_q34',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_094b_bw_nA_q34.json", "w") as f:
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


# ============ q=3 ============
q = 3
g = 1.0 / q
nA_list_q3 = [3, 4, 5, 6]

for nA in nA_list_q3:
    n = 2 * nA
    dim = q**n
    dimA = q**nA
    key = f"q3_nA{nA}"
    print(f"\n=== q={q}, nA={nA}, n={n}, dim={dim}, dimA={dimA} ===")
    t1 = time.time()

    if dim > 1200000:
        print(f"  SKIP: dim={dim} too large")
        continue

    H = build_sq_potts_periodic(n, q, g)
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
           quantity='bw_R2_nA_scaling', value=R2, method='094b_nA_scan')
    save()

# ============ q=4 ============
q = 4
g = 1.0 / q
nA_list_q4 = [3, 4, 5]

for nA in nA_list_q4:
    n = 2 * nA
    dim = q**n
    dimA = q**nA
    key = f"q4_nA{nA}"
    print(f"\n=== q={q}, nA={nA}, n={n}, dim={dim}, dimA={dimA} ===")
    t1 = time.time()

    if dim > 10000000:
        print(f"  SKIP: dim={dim} too large even for GPU")
        continue

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
           quantity='bw_R2_nA_scaling', value=R2, method='094b_nA_scan')
    save()

print(f"\nTotal time: {time.time()-t0:.1f}s")
save()
print("Done!")
