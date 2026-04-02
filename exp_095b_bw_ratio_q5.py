#!/usr/bin/env python3
"""Sprint 095b: BW R² at fixed nA with varying n for q=5 (periodic).

q=5 has walking (complex CFT). Sprint 094 showed nA*=4 (threshold shifted
down by 1 vs q=2). Test: does increasing n at fixed nA improve BW R²?

q=5 periodic: nA=3 with n=7,8,9,10; nA=4 with n=9,10.
Dim limit: 5^10 ≈ 10M (GPU feasible).
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record
import json, time

t0 = time.time()
q = 5
g = 1.0 / q

results = {
    'experiment': '095b_bw_ratio_q5',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'g_c': g,
    'data': {},
}

def save():
    with open("results/sprint_095b_bw_ratio_q5.json", "w") as f:
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


# q=5: dim grows as 5^n. Feasible: n≤10 (GPU, ~10M).
# nA must be < n/2.
scan = {
    3: [7, 8, 9, 10],   # nA=3: n=7..10, dim up to 9.8M
    4: [9, 10],          # nA=4: n=9..10, dim up to 9.8M
}

print(f"Sprint 095b: BW R² vs n at fixed nA for q={q} (periodic)")
print("Testing whether walking amplification depends on nA/n ratio")
print("=" * 70, flush=True)

for nA in sorted(scan.keys()):
    n_list = scan[nA]
    print(f"\n{'='*60}")
    print(f"nA = {nA} (dimA = {q**nA})")
    print(f"{'='*60}", flush=True)

    results['data'][f'nA{nA}'] = {}

    for n in n_list:
        dim = q**n
        ratio = nA / n
        key = f"n{n}"
        print(f"\n  n={n}, dim={dim}, nA/n={ratio:.3f}", flush=True)
        t1 = time.time()

        H = build_sq_potts_periodic(n, q, g)
        dt_build = time.time() - t1
        print(f"    H build: {dt_build:.1f}s", flush=True)

        if dt_build > 200:
            print(f"    SKIP: H build too slow")
            continue

        evals, evecs = eigsh(H, k=1, which='SA')
        psi = evecs[:, 0]
        E0 = float(evals[0])
        dt_diag = time.time() - t1
        print(f"    E0 = {E0:.6f}, diag: {dt_diag:.1f}s", flush=True)

        rho_A = get_rho_A(psi, n, q, nA)
        H_E = entanglement_hamiltonian(rho_A)
        H_BW = build_H_BW(q, n, nA, g)
        R2, alpha = frobenius_R2(H_E, H_BW)
        one_minus_R2 = 1.0 - R2

        dt_total = time.time() - t1
        print(f"    BW R²={R2:.8f}, 1-R²={one_minus_R2:.2e}, alpha={alpha:.4f}, time={dt_total:.1f}s")

        entry = {
            'nA': nA, 'n': n, 'dim': dim, 'dimA': q**nA,
            'ratio_nA_n': float(ratio),
            'E0': E0, 'R2': float(R2), 'one_minus_R2': float(one_minus_R2),
            'alpha': float(alpha), 'time': dt_total,
        }
        results['data'][f'nA{nA}'][key] = entry
        record(sprint=95, model='sq_potts', q=q, n=n,
               quantity='bw_R2_ratio', value=R2,
               method=f'095b_nA{nA}_n{n}')
        save()

# Summary
print(f"\n\n{'='*70}")
print("SUMMARY: 1-R²(nA, n) for q=5")
for nA in sorted(scan.keys()):
    entries = results['data'].get(f'nA{nA}', {})
    print(f"\nnA={nA}:")
    for key in sorted(entries.keys()):
        if key == 'fit':
            continue
        e = entries[key]
        print(f"  n={e['n']:>2}, nA/n={e['ratio_nA_n']:.3f}: 1-R²={e['one_minus_R2']:.2e}, alpha={e['alpha']:.4f}")

# Compare q=2 vs q=5 at same nA
print(f"\n\nCOMPARISON: q=5 walking amplification at different nA/n ratios")
print("Sprint 094 (n=2*nA): q=2 nA=3: 2.9e-4, q=5 nA=3: 1.0e-3 → ratio 3.5x")
print("Sprint 094 (n=2*nA): q=2 nA=4: 5.3e-4, q=5 nA=4: 3.4e-2 → ratio 64x")

save()
print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done!")
