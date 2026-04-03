#!/usr/bin/env python3
"""Sprint 097c: Operator compactness for q=3 n=10 using clock-shift basis.

For q=3, the complete single-site operator basis is the 9 generalized
Gell-Mann matrices (or equivalently: I, X, X², Z, Z², XZ, XZ², X²Z, X²Z²).
We use the clock-shift basis: X^a Z^b for a,b in {0,...,q-1}.

Full basis for nA sites: (q²)^nA operators. For nA=4: 6561 operators.
Test whether H_E compactness pattern is universal across q.
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record
import json, time
from itertools import product

t0 = time.time()
q = 3
g = 1.0 / q
n = 10

results = {
    'experiment': '097c_q3_compactness',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'n': n, 'g_c': g,
    'data': {},
}

def save():
    with open("results/sprint_097c_q3_compactness.json", "w") as f:
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
    return H_E  # complex for q>=3


# Build single-site clock and shift operators for Z_q
def clock_X(q):
    """Shift operator: X|k> = |k+1 mod q>"""
    X = np.zeros((q, q), dtype=complex)
    for k in range(q):
        X[(k + 1) % q, k] = 1.0
    return X

def clock_Z(q):
    """Clock operator: Z|k> = omega^k |k>"""
    omega = np.exp(2j * np.pi / q)
    Z = np.diag([omega**k for k in range(q)])
    return Z


def build_clock_shift_basis(q):
    """Build q² single-site operators: X^a Z^b for a,b = 0,...,q-1.

    These form an orthogonal basis: Tr((X^a Z^b)† X^c Z^d) = q δ_{ac} δ_{bd}.
    """
    X = clock_X(q)
    Z = clock_Z(q)
    I = np.eye(q, dtype=complex)

    basis = {}
    for a in range(q):
        for b in range(q):
            op = np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b)
            label = f'X{a}Z{b}'
            basis[(a, b)] = (label, op)
    return basis


def classify_cs_operator(q, site_indices):
    """Classify a multi-site clock-shift operator.

    site_indices: list of (a_k, b_k) for each site.
    Returns: body order, range, type.
    Type: 'D' (diagonal, all b only), 'S' (shift, all a only),
          'M' (mixed, has both a>0 and b>0 on different sites),
          'SM' (same-site mixed, a>0 and b>0 on same site).
    """
    nA = len(site_indices)
    non_trivial = [(k, a, b) for k, (a, b) in enumerate(site_indices) if a != 0 or b != 0]
    body = len(non_trivial)
    if body == 0:
        return 0, -1, 'I'
    sites = [k for k, a, b in non_trivial]
    rng = sites[-1] - sites[0] if body >= 2 else 0

    has_shift = any(a > 0 for _, a, b in non_trivial)
    has_clock = any(b > 0 for _, a, b in non_trivial)
    has_same_site_mixed = any(a > 0 and b > 0 for _, a, b in non_trivial)

    if has_same_site_mixed:
        ptype = 'SM'  # same-site mixed (XZ on one site)
    elif has_shift and has_clock:
        ptype = 'M'  # mixed across sites
    elif has_shift:
        ptype = 'S'  # pure shift
    elif has_clock:
        ptype = 'D'  # pure diagonal (clock)
    else:
        ptype = 'I'

    # Is it a BW-type operator?
    # BW has: NN Potts bonds (delta = sum_k Z^k Z^{-k}) and single-site fields (X + X†)
    is_bw = False
    if body == 1:
        _, a, b = non_trivial[0]
        if b == 0 and a > 0:  # pure shift = field type
            is_bw = True
    elif body == 2 and rng == 1:
        _, a1, b1 = non_trivial[0]
        _, a2, b2 = non_trivial[1]
        if a1 == 0 and a2 == 0 and b1 > 0 and b2 > 0:
            if (b1 + b2) % q == 0:  # Z^k Z^{-k} = Potts bond component
                is_bw = True

    return body, rng, ptype, is_bw


# Build H, get ground state
print(f"Sprint 097c: Clock-shift decomposition for q={q} n={n}")
print(f"g_c = {g:.4f} (periodic)", flush=True)

H = build_sq_potts_periodic(n, q, g)
evals_H, evecs_H = eigsh(H, k=1, which='SA')
psi = evecs_H[:, 0]
E0 = float(evals_H[0])
print(f"E0 = {E0:.8f}, time: {time.time()-t0:.1f}s\n", flush=True)

cs_basis_1site = build_clock_shift_basis(q)

# Verify orthogonality
for (a1, b1), (l1, o1) in cs_basis_1site.items():
    for (a2, b2), (l2, o2) in cs_basis_1site.items():
        tr = np.trace(o1.conj().T @ o2)
        expected = q if (a1 == a2 and b1 == b2) else 0
        assert abs(tr - expected) < 1e-10, f"Orthogonality fail: {l1},{l2}: {tr} vs {expected}"
print("Clock-shift basis orthogonality verified.\n")


nA_list = [3, 4]

for nA in nA_list:
    dimA = q**nA
    n_ops_total = q**(2 * nA) - 1  # excluding identity
    t1 = time.time()

    print(f"{'='*70}")
    print(f"nA={nA}, dimA={dimA}, total basis ops={n_ops_total}, nA/n={nA/n:.3f}", flush=True)

    rho_A = get_rho_A(psi, n, q, nA)
    H_E = entanglement_hamiltonian(rho_A)

    # Full clock-shift decomposition
    coeffs = []
    for multi_idx in product(range(q * q), repeat=nA):
        # Convert flat index to (a, b) pairs
        site_indices = [(idx // q, idx % q) for idx in multi_idx]
        if all(a == 0 and b == 0 for a, b in site_indices):
            continue  # skip identity

        # Build multi-site operator
        ops = [cs_basis_1site[(a, b)][1] for a, b in site_indices]
        op = ops[0]
        for k in range(1, nA):
            op = np.kron(op, ops[k])

        # Coefficient: Tr(O† · H_E) / dimA
        c = np.trace(op.conj().T @ H_E) / dimA

        if abs(c) < 1e-14:
            continue

        body, rng, ptype, is_bw = classify_cs_operator(q, site_indices)
        label = '_'.join(cs_basis_1site[(a, b)][0] for a, b in site_indices)

        coeffs.append({
            'label': label,
            'coeff_abs_sq': float(abs(c)**2 * dimA),  # Frobenius contribution
            'coeff_re': float(c.real),
            'coeff_im': float(c.imag),
            'range': rng,
            'body': body,
            'type': ptype,
            'is_bw': is_bw,
        })

    # Sort by Frobenius contribution
    coeffs.sort(key=lambda x: x['coeff_abs_sq'], reverse=True)
    total_frob = sum(c['coeff_abs_sq'] for c in coeffs)

    cum_frob = np.cumsum([c['coeff_abs_sq'] for c in coeffs])
    cum_R2 = cum_frob / total_frob

    # Milestones
    milestones = {}
    for target in [0.90, 0.95, 0.99, 0.999, 0.9999]:
        idx = np.searchsorted(cum_R2, target)
        milestones[f'{target:.4f}'] = int(idx + 1) if idx < len(cum_R2) else len(cum_R2)

    # Participation ratio
    c_sq = np.array([c['coeff_abs_sq'] for c in coeffs])
    PR = float(np.sum(c_sq)**2 / np.sum(c_sq**2))

    # BW stats
    n_bw = sum(1 for c in coeffs if c['is_bw'])
    bw_frac = sum(c['coeff_abs_sq'] for c in coeffs if c['is_bw']) / total_frob

    # Type breakdown
    from collections import defaultdict
    type_frob = defaultdict(float)
    type_count = defaultdict(int)
    for c in coeffs:
        type_frob[c['type']] += c['coeff_abs_sq']
        type_count[c['type']] += 1
    type_frac = {t: v / total_frob for t, v in type_frob.items()}

    dt = time.time() - t1
    print(f"  Non-zero operators: {len(coeffs)} / {n_ops_total}")
    print(f"  Participation ratio: {PR:.1f}")
    print(f"  BW operators: {n_bw}, BW fraction: {bw_frac:.6f}")
    print(f"  Milestones:")
    for target, count in milestones.items():
        print(f"    R²≥{target}: {count} operators ({count/n_ops_total*100:.1f}%)")

    print(f"  Type breakdown:")
    for t in sorted(type_frac.keys(), key=lambda x: type_frac[x], reverse=True):
        print(f"    {t:>4}: {type_frac[t]:.6f} ({type_count[t]} ops)")

    print(f"  Top 10:")
    for i, c in enumerate(coeffs[:10]):
        bw_tag = " [BW]" if c['is_bw'] else ""
        print(f"    {i+1}. {c['label'][:30]:>30} |c|²f={c['coeff_abs_sq']/total_frob:.4f} "
              f"cum={cum_R2[i]:.4f} b={c['body']} r={c['range']} {c['type']}{bw_tag}")

    print(f"  Time: {dt:.1f}s")

    # Non-BW analysis
    non_bw = [c for c in coeffs if not c['is_bw']]
    non_bw_total = sum(c['coeff_abs_sq'] for c in non_bw)
    non_bw_frac = non_bw_total / total_frob
    if non_bw:
        nb_sq = np.array([c['coeff_abs_sq'] for c in non_bw])
        nb_PR = float(np.sum(nb_sq)**2 / np.sum(nb_sq**2))
    else:
        nb_PR = 0

    print(f"  Non-BW: {len(non_bw)} ops, fraction={non_bw_frac:.6f}, PR={nb_PR:.1f}")

    entry = {
        'nA': nA, 'dimA': dimA,
        'n_nonzero': len(coeffs),
        'n_total_basis': n_ops_total,
        'participation_ratio': float(PR),
        'n_bw_ops': n_bw,
        'bw_fraction': float(bw_frac),
        'non_bw_fraction': float(non_bw_frac),
        'non_bw_n': len(non_bw),
        'non_bw_PR': float(nb_PR),
        'milestones': milestones,
        'type_fraction': {t: float(v) for t, v in type_frac.items()},
        'top_10': [
            {'label': c['label'][:40], 'frac': c['coeff_abs_sq']/total_frob,
             'body': c['body'], 'range': c['range'], 'type': c['type'], 'is_bw': c['is_bw']}
            for c in coeffs[:10]
        ],
        'time': dt,
    }
    results['data'][f'nA{nA}'] = entry

    record(sprint=97, model='sq_potts', q=q, n=n,
           quantity='HE_participation_ratio', value=float(PR),
           method=f'097c_nA{nA}')
    record(sprint=97, model='sq_potts', q=q, n=n,
           quantity='HE_bw_fraction', value=float(bw_frac),
           method=f'097c_nA{nA}')
    save()


# Also do nA=5 for q=3 (at the BW threshold from Sprint 096b: 1-R²=0.211)
# dimA=3^5=243, basis=9^5=59049 — too large for full decomposition
# Instead, compute BW R² and non-BW fraction via Frobenius difference
nA = 5
dimA = q**nA
t1 = time.time()
print(f"\n{'='*70}")
print(f"nA={nA} (threshold): BW R² only (full basis too large: {q**(2*nA)} ops)", flush=True)

rho_A = get_rho_A(psi, n, q, nA)
H_E = entanglement_hamiltonian(rho_A)

# Build BW H
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
    op = np.zeros((dimA, dimA), dtype=complex)
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

H_BW = np.zeros((dimA, dimA), dtype=complex)
for x in range(nA - 1):
    beta = bw_envelope_bond(n, nA, x)
    delta = potts_delta_bond(q, nA, x, x + 1)
    H_BW -= beta * delta
for x in range(nA):
    beta = bw_envelope_site(n, nA, x)
    field = clock_field_site(q, nA, x)
    H_BW -= g * beta * field

# BW R²
H_E_tl = H_E - np.trace(H_E) / dimA * np.eye(dimA)
H_BW_tl = H_BW - np.trace(H_BW) / dimA * np.eye(dimA)
alpha = np.real(np.sum(H_E_tl.conj() * H_BW_tl)) / np.real(np.sum(H_BW_tl.conj() * H_BW_tl))
res = H_E_tl - alpha * H_BW_tl
R2 = 1 - np.real(np.sum(res.conj() * res)) / np.real(np.sum(H_E_tl.conj() * H_E_tl))

# Free-fit all Potts operators (NN and beyond)
all_potts_ops = []
for dist in range(1, nA):
    for si in range(nA - dist):
        sj = si + dist
        op = potts_delta_bond(q, nA, si, sj)
        op_tl = op - np.trace(op) / dimA * np.eye(dimA)
        all_potts_ops.append(op_tl.flatten())
for site in range(nA):
    op = clock_field_site(q, nA, site)
    op_tl = op - np.trace(op) / dimA * np.eye(dimA)
    all_potts_ops.append(np.real(op_tl.flatten()))
    all_potts_ops.append(np.imag(op_tl.flatten()))

op_vecs = np.array(all_potts_ops)
target_vec = np.real(H_E_tl.flatten())
ATA = op_vecs @ op_vecs.T
ATb = op_vecs @ target_vec
try:
    coeffs_fit = np.linalg.lstsq(ATA, ATb, rcond=None)[0]
    fit_vec = coeffs_fit @ op_vecs
    res_var = np.sum((target_vec - fit_vec)**2)
    total_var = np.sum(target_vec**2)
    R2_potts = 1 - res_var / total_var
except:
    R2_potts = -1

dt = time.time() - t1
print(f"  BW R² = {R2:.6f}, 1-R² = {1-R2:.2e}")
print(f"  Free-fit all Potts: R² = {R2_potts:.6f}, 1-R² = {1-R2_potts:.2e}")
print(f"  Non-Potts fraction: {1-R2_potts:.4f}")
print(f"  Time: {dt:.1f}s")

results['data']['nA5_threshold'] = {
    'nA': 5, 'dimA': dimA,
    'bw_R2': float(R2),
    'bw_1mR2': float(1-R2),
    'potts_free_R2': float(R2_potts),
    'non_potts_fraction': float(1-R2_potts),
    'time': dt,
}

record(sprint=97, model='sq_potts', q=q, n=n,
       quantity='bw_R2', value=float(R2), method='097c_nA5')
save()

# Summary
print(f"\n\n{'='*70}")
print(f"SUMMARY: q=3 n=10 operator compactness")
print(f"{'nA':>4} {'basis':>8} {'nonzero':>8} {'PR':>8} {'N(99%)':>8} {'BW frac':>10}")
print("-" * 55)
for nA_key in ['nA3', 'nA4']:
    d = results['data'][nA_key]
    print(f"{d['nA']:>4} {d['n_total_basis']:>8} {d['n_nonzero']:>8} "
          f"{d['participation_ratio']:>8.1f} "
          f"{d['milestones'].get('0.9900', '—'):>8} "
          f"{d['bw_fraction']:>10.6f}")
d5 = results['data']['nA5_threshold']
print(f"   5    59048      —        —        — "
      f"  {d5['bw_R2']:>10.6f}")

save()
print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done!")
