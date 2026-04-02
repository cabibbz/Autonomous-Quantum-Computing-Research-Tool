#!/usr/bin/env python3
"""Sprint 096a: Entanglement spectrum vs nA at fixed n for q=2 (periodic).

Connect entanglement spectrum structure to BW fidelity threshold.
n=14 periodic chain at g_c=1/2. Subsystems nA=3,4,5,6,7.
Track: Schmidt spectrum, tail weight, entanglement gap, BW R².
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record
import json, time

t0 = time.time()
q = 2
g = 1.0 / q
n = 14

results = {
    'experiment': '096a_spectrum_vs_nA_q2',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'n': n, 'g_c': g,
    'data': {},
}

def save():
    with open("results/sprint_096a_spectrum_vs_nA_q2.json", "w") as f:
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


# Build H and get ground state
print(f"Sprint 096a: Entanglement spectrum vs nA for q={q}, n={n} (periodic)")
print(f"dim = {q**n}", flush=True)

H = build_sq_potts_periodic(n, q, g)
print(f"H built: {time.time()-t0:.1f}s", flush=True)

evals_H, evecs_H = eigsh(H, k=1, which='SA')
psi = evecs_H[:, 0]
E0 = float(evals_H[0])
print(f"E0 = {E0:.8f}, diag: {time.time()-t0:.1f}s", flush=True)

results['E0'] = E0

# For each nA, compute entanglement spectrum and BW R²
nA_list = [3, 4, 5, 6, 7]

for nA in nA_list:
    ratio = nA / n
    dimA = q**nA
    print(f"\n{'='*60}")
    print(f"nA = {nA}, dimA = {dimA}, nA/n = {ratio:.3f}", flush=True)

    t1 = time.time()

    # Reduced density matrix
    rho_A = get_rho_A(psi, n, q, nA)

    # Schmidt spectrum (eigenvalues of rho_A, sorted descending)
    schmidt = np.sort(np.linalg.eigvalsh(rho_A))[::-1]
    schmidt = np.clip(schmidt, 0, None)  # numerical cleanup
    schmidt_norm = schmidt / schmidt.sum()  # normalize

    # Entanglement entropy
    S_vN = float(-np.sum(schmidt_norm[schmidt_norm > 1e-30] *
                          np.log(schmidt_norm[schmidt_norm > 1e-30])))

    # For q=2: ground state has (q-1)=1 multiplet above lev0
    # Classify: lev0 = largest eigenvalue, lev1 = next (q-1) eigenvalues, tail = rest
    lev0_weight = float(schmidt_norm[0])
    lev1_weight = float(np.sum(schmidt_norm[1:q]))  # next (q-1) eigenvalues
    tail_weight = float(np.sum(schmidt_norm[q:]))

    # Entropy contributions
    S_lev0 = float(-schmidt_norm[0] * np.log(schmidt_norm[0])) if schmidt_norm[0] > 1e-30 else 0
    S_lev1 = float(-np.sum(schmidt_norm[1:q] * np.log(np.clip(schmidt_norm[1:q], 1e-30, None))))
    S_tail = S_vN - S_lev0 - S_lev1

    # Entanglement gap: -log(lambda_1/lambda_0)
    if schmidt_norm[1] > 1e-30:
        ent_gap = float(-np.log(schmidt_norm[1] / schmidt_norm[0]))
    else:
        ent_gap = float('inf')

    # Number of significant eigenvalues (> 1e-6)
    n_significant = int(np.sum(schmidt_norm > 1e-6))

    # Participation ratio (inverse purity of spectrum)
    participation = float(1.0 / np.sum(schmidt_norm**2))

    # BW R²
    H_E = entanglement_hamiltonian(rho_A)
    H_BW = build_H_BW(q, n, nA, g)
    R2, alpha = frobenius_R2(H_E, H_BW)
    one_minus_R2 = 1.0 - R2

    dt = time.time() - t1

    print(f"  S_vN = {S_vN:.6f}")
    print(f"  lev0 weight: {lev0_weight:.6f}, lev1: {lev1_weight:.6f}, tail: {tail_weight:.6e}")
    print(f"  %S(lev0): {S_lev0/S_vN:.4f}, %S(lev1): {S_lev1/S_vN:.4f}, %S(tail): {S_tail/S_vN:.4f}")
    print(f"  Entanglement gap: {ent_gap:.6f}")
    print(f"  N_significant: {n_significant}, Participation: {participation:.4f}")
    print(f"  BW R² = {R2:.8f}, 1-R² = {one_minus_R2:.2e}, alpha = {alpha:.4f}")
    print(f"  Time: {dt:.1f}s", flush=True)

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
        'top_10_schmidt': [float(x) for x in schmidt_norm[:10]],
        'time': dt,
    }
    results['data'][f'nA{nA}'] = entry

    record(sprint=96, model='sq_potts', q=q, n=n,
           quantity='bw_R2_spectrum', value=R2,
           method=f'096a_nA{nA}', notes=f'tail={tail_weight:.2e}')
    save()

# Summary table
print(f"\n\n{'='*70}")
print("SUMMARY: Spectrum + BW at nA=3..7, q=2, n=14")
print(f"{'nA':>4} {'nA/n':>6} {'tail_w':>10} {'%S(tail)':>10} {'ent_gap':>10} {'N_sig':>6} {'1-R²':>10}")
print("-" * 62)
for nA in nA_list:
    d = results['data'][f'nA{nA}']
    print(f"{nA:>4} {d['ratio_nA_n']:>6.3f} {d['tail_weight']:>10.2e} "
          f"{d['pct_S_tail']:>10.4f} {d['ent_gap']:>10.4f} "
          f"{d['n_significant']:>6} {d['one_minus_R2']:>10.2e}")

# Log-log correlation: tail_weight vs 1-R²
tails = [results['data'][f'nA{nA}']['tail_weight'] for nA in nA_list]
omr = [results['data'][f'nA{nA}']['one_minus_R2'] for nA in nA_list]
tails_arr = np.array(tails)
omr_arr = np.array(omr)
mask = (tails_arr > 0) & (omr_arr > 0)
if np.sum(mask) >= 3:
    log_t = np.log(tails_arr[mask])
    log_o = np.log(omr_arr[mask])
    coeffs = np.polyfit(log_t, log_o, 1)
    pred = np.polyval(coeffs, log_t)
    ss_res = np.sum((log_o - pred)**2)
    ss_tot = np.sum((log_o - np.mean(log_o))**2)
    fit_R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    print(f"\nCorrelation: log(1-R²) ~ {coeffs[0]:.2f} * log(tail_weight), R² = {fit_R2:.4f}")
    results['correlation_tail_R2'] = {'slope': float(coeffs[0]), 'R2': float(fit_R2)}

# Same for entanglement gap
gaps = [results['data'][f'nA{nA}']['ent_gap'] for nA in nA_list]
gaps_arr = np.array(gaps)
mask2 = (gaps_arr > 0) & (gaps_arr < 100) & (omr_arr > 0)
if np.sum(mask2) >= 3:
    coeffs2 = np.polyfit(gaps_arr[mask2], np.log(omr_arr[mask2]), 1)
    pred2 = np.polyval(coeffs2, gaps_arr[mask2])
    log_o2 = np.log(omr_arr[mask2])
    ss_res2 = np.sum((log_o2 - pred2)**2)
    ss_tot2 = np.sum((log_o2 - np.mean(log_o2))**2)
    fit_R2_2 = 1 - ss_res2 / ss_tot2 if ss_tot2 > 0 else 0
    print(f"Correlation: log(1-R²) ~ {coeffs2[0]:.3f} * ent_gap, R² = {fit_R2_2:.4f}")
    results['correlation_gap_R2'] = {'slope': float(coeffs2[0]), 'R2': float(fit_R2_2)}

save()
print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done!")
