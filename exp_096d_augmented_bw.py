#!/usr/bin/env python3
"""Sprint 096d: Augmented BW with range-2 operators for q=2 n=14.

BW has only range-0,1 operators. 096c showed the residual is dominated by
long-range operators. Test: if we add NNN (range-2) operators to BW
(fitted, not from BW theory), does the threshold shift?

Also: cumulative variance explained by fitting operators up to range r.
This shows EXACTLY how many terms are needed to capture H_E.
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
    'experiment': '096d_augmented_bw',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'n': n, 'g_c': g,
    'data': {},
}

def save():
    with open("results/sprint_096d_augmented_bw.json", "w") as f:
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


# Build H, get ground state
print(f"Sprint 096d: Augmented BW fitting by operator range")
print(f"q={q}, n={n} (periodic)", flush=True)

H = build_sq_potts_periodic(n, q, g)
evals_H, evecs_H = eigsh(H, k=1, which='SA')
psi = evecs_H[:, 0]
E0 = float(evals_H[0])
print(f"E0 = {E0:.8f}, time: {time.time()-t0:.1f}s\n", flush=True)


def build_operator_basis(q, nA, max_range):
    """Build basis of Potts-type operators up to given range.

    Returns list of (label, operator_matrix, range).
    Operators: delta(si,sj) bonds and shift(site) fields at all positions.
    """
    ops = []
    dimA = q**nA

    # Range 0: single-site field operators
    for site in range(nA):
        op = clock_field_site(q, nA, site)
        ops.append((f'F_{site}', op, 0))

    # Range 1+: delta bonds at various distances
    for dist in range(1, max_range + 1):
        for si in range(nA - dist):
            sj = si + dist
            op = potts_delta_bond(q, nA, si, sj)
            ops.append((f'D_{si}_{sj}', op, dist))

    return ops


nA_list = [3, 4, 5, 6, 7]

for nA in nA_list:
    dimA = q**nA
    ratio = nA / n
    print(f"{'='*60}")
    print(f"nA={nA}, dimA={dimA}, nA/n={ratio:.3f}", flush=True)

    rho_A = get_rho_A(psi, n, q, nA)
    H_E = entanglement_hamiltonian(rho_A)

    # Traceless H_E
    H_E_tl = H_E - np.trace(H_E) / dimA * np.eye(dimA)
    total_var = np.sum(H_E_tl**2)

    # Standard BW (1-param)
    H_BW = build_H_BW(q, n, nA, g)
    H_BW_tl = H_BW - np.trace(H_BW) / dimA * np.eye(dimA)
    alpha_bw = np.sum(H_E_tl * H_BW_tl) / np.sum(H_BW_tl**2)
    res_bw = H_E_tl - alpha_bw * H_BW_tl
    var_bw = np.sum(res_bw**2)
    R2_bw = 1 - var_bw / total_var

    print(f"  BW (1-param): R²={R2_bw:.8f}, 1-R²={1-R2_bw:.2e}")

    # Augmented: fit each Potts-type operator independently (not constrained by envelope)
    # Up to each max range, do least-squares fit
    max_range_test = min(nA - 1, 6)

    entry = {
        'nA': nA, 'dimA': dimA,
        'bw_R2': float(R2_bw),
        'bw_1mR2': float(1 - R2_bw),
        'augmented': {},
    }

    for max_r in range(0, max_range_test + 1):
        ops = build_operator_basis(q, nA, max_r)

        # Flatten operators into design matrix
        n_ops = len(ops)
        # Vectorize: flatten each op minus trace, compute coefficients via least squares
        op_vecs = np.zeros((n_ops, dimA * dimA))
        for i, (label, op, r) in enumerate(ops):
            op_tl = op - np.trace(op) / dimA * np.eye(dimA)
            op_vecs[i] = op_tl.flatten()

        target_vec = H_E_tl.flatten()

        # Least squares: target = sum_i c_i * op_i
        # Normal equations: (A^T A) c = A^T b
        ATA = op_vecs @ op_vecs.T
        ATb = op_vecs @ target_vec

        try:
            coeffs = np.linalg.solve(ATA, ATb)
            fit_vec = coeffs @ op_vecs
            res_var = np.sum((target_vec - fit_vec)**2)
            R2_aug = 1 - res_var / total_var
        except np.linalg.LinAlgError:
            # Singular — use pseudoinverse
            coeffs = np.linalg.lstsq(ATA, ATb, rcond=None)[0]
            fit_vec = coeffs @ op_vecs
            res_var = np.sum((target_vec - fit_vec)**2)
            R2_aug = 1 - res_var / total_var

        print(f"  max_range={max_r}: {n_ops} ops, R²={R2_aug:.8f}, 1-R²={1-R2_aug:.2e}")

        entry['augmented'][f'r{max_r}'] = {
            'max_range': max_r,
            'n_ops': n_ops,
            'R2': float(R2_aug),
            'one_minus_R2': float(1 - R2_aug),
        }

    results['data'][f'nA{nA}'] = entry
    save()

# Summary table
print(f"\n\n{'='*70}")
print("SUMMARY: R² with Potts operators up to range r (free coefficients)")
print(f"{'nA':>4} {'BW(1p)':>10}", end="")
for r in range(7):
    print(f"  {'r≤'+str(r):>8}", end="")
print()
print("-" * (16 + 10 * 7))

for nA in nA_list:
    d = results['data'][f'nA{nA}']
    print(f"{nA:>4} {d['bw_1mR2']:>10.2e}", end="")
    for r in range(7):
        aug = d['augmented'].get(f'r{r}')
        if aug:
            print(f"  {aug['one_minus_R2']:>8.2e}", end="")
        else:
            print(f"  {'—':>8}", end="")
    print()

save()
print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done!")
