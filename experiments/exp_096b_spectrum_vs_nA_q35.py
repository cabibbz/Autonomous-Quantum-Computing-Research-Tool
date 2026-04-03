#!/usr/bin/env python3
"""Sprint 096b: Entanglement spectrum vs nA for q=3,5 (periodic).

Same analysis as 096a but for q=3 (n=10) and q=5 (n=8).
Check if smooth-spectrum-sharp-BW pattern is universal.
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record
import json, time

t0 = time.time()

results = {
    'experiment': '096b_spectrum_vs_nA_q35',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_096b_spectrum_vs_nA_q35.json", "w") as f:
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


def bw_envelope_bond(N, nA, x):
    return (N / np.pi) * np.sin(np.pi * (x + 1) / N) * \
           np.sin(np.pi * (nA - x - 1) / N) / np.sin(np.pi * nA / N)


def bw_envelope_site(N, nA, x):
    return (N / np.pi) * np.sin(np.pi * (x + 0.5) / N) * \
           np.sin(np.pi * (nA - x - 0.5) / N) / np.sin(np.pi * nA / N)


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


# q=3 n=10: nA=3,4,5
# q=5 n=8: nA=3,4
configs = [
    (3, 10, [3, 4, 5]),
    (5, 8, [3, 4]),
]

for q, n, nA_list in configs:
    g = 1.0 / q
    dim = q**n
    print(f"\n{'#'*60}")
    print(f"q={q}, n={n}, dim={dim}, g_c={g}")
    print(f"{'#'*60}", flush=True)

    t1 = time.time()
    H = build_sq_potts_periodic(n, q, g)
    print(f"H built: {time.time()-t1:.1f}s", flush=True)

    evals_H, evecs_H = eigsh(H, k=1, which='SA')
    psi = evecs_H[:, 0]
    E0 = float(evals_H[0])
    print(f"E0 = {E0:.8f}, diag: {time.time()-t1:.1f}s", flush=True)

    q_key = f'q{q}'
    results['data'][q_key] = {'n': n, 'q': q, 'g_c': float(g), 'E0': E0, 'nA_data': {}}

    for nA in nA_list:
        ratio = nA / n
        dimA = q**nA
        print(f"\n  nA={nA}, dimA={dimA}, nA/n={ratio:.3f}", flush=True)

        rho_A = get_rho_A(psi, n, q, nA)
        schmidt = np.sort(np.linalg.eigvalsh(rho_A))[::-1]
        schmidt = np.clip(schmidt, 0, None)
        schmidt_norm = schmidt / schmidt.sum()

        S_vN = float(-np.sum(schmidt_norm[schmidt_norm > 1e-30] *
                              np.log(schmidt_norm[schmidt_norm > 1e-30])))

        lev0_weight = float(schmidt_norm[0])
        lev1_weight = float(np.sum(schmidt_norm[1:q]))
        tail_weight = float(np.sum(schmidt_norm[q:]))

        S_lev0 = float(-schmidt_norm[0] * np.log(schmidt_norm[0])) if schmidt_norm[0] > 1e-30 else 0
        S_lev1 = float(-np.sum(schmidt_norm[1:q] * np.log(np.clip(schmidt_norm[1:q], 1e-30, None))))
        S_tail = S_vN - S_lev0 - S_lev1

        ent_gap = float(-np.log(schmidt_norm[1] / schmidt_norm[0])) if schmidt_norm[1] > 1e-30 else float('inf')
        n_significant = int(np.sum(schmidt_norm > 1e-6))
        participation = float(1.0 / np.sum(schmidt_norm**2))

        H_E = entanglement_hamiltonian(rho_A)
        H_BW = build_H_BW(q, n, nA, g)
        R2, alpha = frobenius_R2(H_E, H_BW)
        one_minus_R2 = 1.0 - R2

        print(f"    S_vN={S_vN:.6f}, tail_w={tail_weight:.2e}, %S(tail)={S_tail/S_vN:.4f}")
        print(f"    ent_gap={ent_gap:.4f}, N_sig={n_significant}, participation={participation:.4f}")
        print(f"    BW R²={R2:.8f}, 1-R²={one_minus_R2:.2e}, alpha={alpha:.4f}")

        entry = {
            'nA': nA, 'dimA': dimA, 'ratio_nA_n': float(ratio),
            'S_vN': float(S_vN),
            'lev0_weight': lev0_weight,
            'lev1_weight': lev1_weight,
            'tail_weight': tail_weight,
            'pct_S_lev0': float(S_lev0/S_vN),
            'pct_S_lev1': float(S_lev1/S_vN),
            'pct_S_tail': float(S_tail/S_vN),
            'ent_gap': ent_gap,
            'n_significant': n_significant,
            'participation': participation,
            'R2': float(R2),
            'one_minus_R2': float(one_minus_R2),
            'alpha': float(alpha),
            'top_10_schmidt': [float(x) for x in schmidt_norm[:min(10, len(schmidt_norm))]],
        }
        results['data'][q_key]['nA_data'][f'nA{nA}'] = entry

        record(sprint=96, model='sq_potts', q=q, n=n,
               quantity='bw_R2_spectrum', value=R2,
               method=f'096b_nA{nA}', notes=f'tail={tail_weight:.2e}')
        save()

# Cross-q summary
print(f"\n\n{'='*70}")
print("CROSS-Q SUMMARY: Spectrum + BW")
print(f"{'q':>3} {'nA':>4} {'nA/n':>6} {'tail_w':>10} {'%S(tail)':>10} {'ent_gap':>10} {'1-R²':>10}")
print("-" * 62)

# Include q=2 data from 096a
for nA in [3, 4, 5, 6, 7]:
    pass  # 096a data is in a different file

for q_key in sorted(results['data'].keys()):
    qd = results['data'][q_key]
    q_val = qd['q']
    for nA_key in sorted(qd['nA_data'].keys()):
        d = qd['nA_data'][nA_key]
        print(f"{q_val:>3} {d['nA']:>4} {d['ratio_nA_n']:>6.3f} {d['tail_weight']:>10.2e} "
              f"{d['pct_S_tail']:>10.4f} {d['ent_gap']:>10.4f} {d['one_minus_R2']:>10.2e}")

save()
print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done!")
