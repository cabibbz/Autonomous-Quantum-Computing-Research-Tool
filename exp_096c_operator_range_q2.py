#!/usr/bin/env python3
"""Sprint 096c: BW residual decomposed by operator range for q=2 n=14.

The entanglement spectrum is smooth across the BW threshold (096a),
so the threshold must come from operator structure. Decompose ||H_res||²
by operator range (max distance between non-identity sites).

BW contains only range-0 (field) and range-1 (NN bond). If range≥2
operators have a threshold growth, that's the mechanism.
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record
import json, time
from itertools import product

t0 = time.time()
q = 2
g = 1.0 / q
n = 14

results = {
    'experiment': '096c_operator_range_q2',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'n': n, 'g_c': g,
    'data': {},
}

def save():
    with open("results/sprint_096c_operator_range_q2.json", "w") as f:
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
    return R2, alpha, H_t - alpha * H_m  # return residual too


def pauli_basis_q2(nA):
    """Generate complete Pauli basis for nA qubits (q=2).

    Returns dict: label -> (operator_matrix, range).
    Range = max site index - min site index among non-I Paulis, or -1 for identity.
    """
    I = np.eye(2)
    X = np.array([[0, 1], [1, 0]], dtype=float)
    Y = np.array([[0, -1], [1, 0]], dtype=float)  # real part
    Z = np.array([[1, 0], [0, -1]], dtype=float)
    paulis = {'I': I, 'X': X, 'Y': Y, 'Z': Z}

    basis = {}
    for labels in product('IXYZ', repeat=nA):
        label = ''.join(labels)
        op = paulis[labels[0]]
        for i in range(1, nA):
            op = np.kron(op, paulis[labels[i]])

        # Compute range: distance between first and last non-I site
        non_I_sites = [i for i, l in enumerate(labels) if l != 'I']
        if len(non_I_sites) == 0:
            rng = -1  # identity
        elif len(non_I_sites) == 1:
            rng = 0  # single-site operator
        else:
            rng = non_I_sites[-1] - non_I_sites[0]

        basis[label] = (op, rng)

    return basis


# Build H, get ground state
print(f"Sprint 096c: Operator range decomposition of BW residual")
print(f"q={q}, n={n} (periodic)", flush=True)

H = build_sq_potts_periodic(n, q, g)
evals_H, evecs_H = eigsh(H, k=1, which='SA')
psi = evecs_H[:, 0]
E0 = float(evals_H[0])
print(f"E0 = {E0:.8f}, H+diag: {time.time()-t0:.1f}s\n", flush=True)

nA_list = [3, 4, 5, 6, 7]

for nA in nA_list:
    dimA = q**nA
    ratio = nA / n
    print(f"{'='*60}")
    print(f"nA={nA}, dimA={dimA}, nA/n={ratio:.3f}", flush=True)
    t1 = time.time()

    rho_A = get_rho_A(psi, n, q, nA)
    H_E = entanglement_hamiltonian(rho_A)
    H_BW = build_H_BW(q, n, nA, g)
    R2, alpha, H_res = frobenius_R2(H_E, H_BW)
    one_minus_R2 = 1.0 - R2

    print(f"  BW R²={R2:.8f}, 1-R²={one_minus_R2:.2e}, alpha={alpha:.4f}")

    # Decompose H_res into Pauli basis and sort by range
    basis = pauli_basis_q2(nA)

    # Traceless part of H_E and H_res (remove identity component)
    H_E_tl = H_E - np.trace(H_E) / dimA * np.eye(dimA)
    H_res_tl = H_res  # already traceless from frobenius_R2

    total_norm_sq = np.sum(H_E_tl**2)
    res_norm_sq = np.sum(H_res_tl**2)

    # Pauli decomposition: c_P = Tr(P·M) / dimA
    range_norm_sq = {}  # range -> sum of |c_P|² * dimA
    range_count = {}
    range_top_ops = {}  # range -> list of (|c|, label)

    for label, (op, rng) in basis.items():
        if label == 'I' * nA:
            continue  # skip identity
        c_res = np.sum(H_res_tl * op) / dimA
        c_sq = c_res**2 * dimA  # Frobenius contribution

        if rng not in range_norm_sq:
            range_norm_sq[rng] = 0.0
            range_count[rng] = 0
            range_top_ops[rng] = []
        range_norm_sq[rng] += c_sq
        range_count[rng] += 1
        if abs(c_res) > 1e-10:
            range_top_ops[rng].append((abs(c_res), label))

    # Sort top ops and keep top 3 per range
    for rng in range_top_ops:
        range_top_ops[rng].sort(reverse=True)
        range_top_ops[rng] = range_top_ops[rng][:3]

    # Print range decomposition
    print(f"  ||H_E||² = {total_norm_sq:.4f}, ||H_res||² = {res_norm_sq:.4f}")
    print(f"  Residual fraction: {res_norm_sq/total_norm_sq:.6f}")

    range_data = {}
    print(f"  {'Range':>6} {'||·||²':>10} {'Fraction':>10} {'N_ops':>6}  Top operators")
    for rng in sorted(range_norm_sq.keys()):
        frac = range_norm_sq[rng] / res_norm_sq if res_norm_sq > 0 else 0
        tops = range_top_ops.get(rng, [])
        top_str = ', '.join(f'{lab}({c:.4f})' for c, lab in tops[:3])
        print(f"  {rng:>6} {range_norm_sq[rng]:>10.4e} {frac:>10.4f} {range_count[rng]:>6}  {top_str}")
        range_data[str(rng)] = {
            'norm_sq': float(range_norm_sq[rng]),
            'fraction_of_residual': float(frac),
            'n_ops': range_count[rng],
            'top_ops': [(float(c), lab) for c, lab in tops],
        }

    dt = time.time() - t1
    print(f"  Time: {dt:.1f}s")

    entry = {
        'nA': nA, 'dimA': dimA, 'ratio_nA_n': float(ratio),
        'R2': float(R2), 'one_minus_R2': float(one_minus_R2),
        'alpha': float(alpha),
        'total_norm_sq': float(total_norm_sq),
        'res_norm_sq': float(res_norm_sq),
        'res_fraction': float(res_norm_sq / total_norm_sq),
        'range_decomposition': range_data,
        'time': dt,
    }
    results['data'][f'nA{nA}'] = entry

    record(sprint=96, model='sq_potts', q=q, n=n,
           quantity='bw_res_range', value=float(res_norm_sq / total_norm_sq),
           method=f'096c_nA{nA}')
    save()

# Summary: how does long-range fraction change?
print(f"\n\n{'='*70}")
print("SUMMARY: Range decomposition of BW residual (fraction of ||H_res||²)")
print(f"{'nA':>4} {'1-R²':>10}", end="")
max_range = max(max(int(r) for r in d['range_decomposition'].keys())
                for d in results['data'].values())
for r in range(0, max_range + 1):
    print(f"  {'r='+str(r):>8}", end="")
print()
print("-" * (16 + 10 * (max_range + 1)))

for nA in nA_list:
    d = results['data'][f'nA{nA}']
    print(f"{nA:>4} {d['one_minus_R2']:>10.2e}", end="")
    for r in range(0, max_range + 1):
        rd = d['range_decomposition'].get(str(r), {})
        frac = rd.get('fraction_of_residual', 0)
        print(f"  {frac:>8.4f}", end="")
    print()

# Key question: does long-range fraction grow at threshold?
print(f"\nLong-range (r≥2) fraction of residual:")
for nA in nA_list:
    d = results['data'][f'nA{nA}']
    lr_frac = sum(rd.get('fraction_of_residual', 0)
                  for r_key, rd in d['range_decomposition'].items()
                  if int(r_key) >= 2)
    print(f"  nA={nA}: {lr_frac:.4f}")

save()
print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done!")
