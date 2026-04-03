#!/usr/bin/env python3
"""Sprint 095a: BW R² at fixed nA with varying n for q=2 (periodic).

Sprint 094 always used n=2*nA (equal bipartition). BW theory predicts
accuracy improves when nA/n → 0. Test: at fixed nA, does BW R² improve
as n grows (nA/n decreases)?

q=2 periodic, nA=3,4,5,6 with n ranging from 2*nA up to feasible limit.
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record
import json, time

t0 = time.time()
q = 2
g = 1.0 / q

results = {
    'experiment': '095a_bw_ratio_q2',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'g_c': g,
    'data': {},
}

def save():
    with open("results/sprint_095a_bw_ratio_q2.json", "w") as f:
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
    """Get reduced density matrix for first nA sites (periodic chain)."""
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
    """BW envelope for periodic chain of length N, subsystem size nA."""
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


# Scan: nA fixed, vary n. Keep nA < n/2 (proper subsystem).
# n must be even for clean half-chain, but nA<n/2 is the real constraint.
scan = {
    3: [7, 8, 9, 10, 12, 14, 16],     # nA=3: n=7..16, dim up to 65536
    4: [9, 10, 12, 14, 16, 18],         # nA=4: n=9..18, dim up to 262144
    5: [11, 12, 14, 16, 18],            # nA=5: n=11..18, dim up to 262144
    6: [13, 14, 16, 18],                # nA=6: n=13..18, dim up to 262144
}

print(f"Sprint 095a: BW R² vs n at fixed nA for q={q} (periodic)")
print("Testing whether BW threshold depends on nA/n ratio")
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

        # Build H, get ground state
        H = build_sq_potts_periodic(n, q, g)
        dt_build = time.time() - t1
        print(f"    H build: {dt_build:.1f}s", flush=True)

        evals, evecs = eigsh(H, k=1, which='SA')
        psi = evecs[:, 0]
        E0 = float(evals[0])
        dt_diag = time.time() - t1
        print(f"    E0 = {E0:.6f}, diag: {dt_diag:.1f}s", flush=True)

        # Get rho_A and compute BW fidelity
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
               method=f'095a_nA{nA}_n{n}')
        save()

# Summary table
print(f"\n\n{'='*70}")
print("SUMMARY: 1-R²(nA, n)")
print(f"{'nA':>4} | ", end="")
all_n = sorted(set(n for ns in scan.values() for n in ns))
for n in all_n:
    print(f"{'n='+str(n):>10}", end="")
print()
print("-" * (6 + 10 * len(all_n)))

for nA in sorted(scan.keys()):
    print(f"{nA:>4} | ", end="")
    for n in all_n:
        key = f"n{n}"
        d = results['data'].get(f'nA{nA}', {}).get(key)
        if d:
            print(f"{d['one_minus_R2']:>10.2e}", end="")
        else:
            print(f"{'—':>10}", end="")
    print()

# For each nA, fit 1-R² vs (nA/n)
print(f"\n\nFIT: 1-R² vs nA/n at each nA")
for nA in sorted(scan.keys()):
    entries = results['data'].get(f'nA{nA}', {})
    if len(entries) < 3:
        continue
    ratios = np.array([e['ratio_nA_n'] for e in entries.values()])
    omr = np.array([e['one_minus_R2'] for e in entries.values()])
    # Log-log fit
    mask = omr > 0
    if np.sum(mask) >= 3:
        log_r = np.log(ratios[mask])
        log_o = np.log(omr[mask])
        coeffs = np.polyfit(log_r, log_o, 1)
        slope = coeffs[0]
        # R² of fit
        pred = np.polyval(coeffs, log_r)
        ss_res = np.sum((log_o - pred)**2)
        ss_tot = np.sum((log_o - np.mean(log_o))**2)
        fit_R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        print(f"  nA={nA}: 1-R² ~ (nA/n)^{slope:.2f}, fit R²={fit_R2:.3f}")
        results['data'][f'nA{nA}']['fit'] = {
            'slope': float(slope), 'fit_R2': float(fit_R2),
        }

save()
print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done!")
