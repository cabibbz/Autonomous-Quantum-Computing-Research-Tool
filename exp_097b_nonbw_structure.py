#!/usr/bin/env python3
"""Sprint 097b: Detailed non-BW operator structure at nA=6 threshold (q=2 n=14).

097a showed H_E goes from 9 operators (nA≤5) to ~1000 (nA=6).
This experiment analyzes the NON-BW content:
- What (body, range) categories dominate?
- Is the non-BW content itself compact or diffuse?
- Can "BW + one correction class" recover most of H_E?
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record
import json, time
from itertools import product
from collections import defaultdict

t0 = time.time()
q = 2
g = 1.0 / q
n = 14

results = {
    'experiment': '097b_nonbw_structure',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'n': n, 'g_c': g,
    'data': {},
}

def save():
    with open("results/sprint_097b_nonbw_structure.json", "w") as f:
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


I2 = np.eye(2)
X2 = np.array([[0, 1], [1, 0]], dtype=float)
Y2 = np.array([[0, -1], [1, 0]], dtype=float)
Z2 = np.array([[1, 0], [0, -1]], dtype=float)
PAULIS = [I2, X2, Y2, Z2]
PAULI_LABELS = ['I', 'X', 'Y', 'Z']


def full_pauli_analysis(H_E, nA):
    """Full Pauli decomposition with detailed classification."""
    dimA = 2**nA
    all_ops = []

    for indices in product(range(4), repeat=nA):
        label = ''.join(PAULI_LABELS[i] for i in indices)
        if all(i == 0 for i in indices):
            continue

        op = PAULIS[indices[0]]
        for k in range(1, nA):
            op = np.kron(op, PAULIS[indices[k]])

        c = np.sum(op * H_E.T) / dimA

        non_I = [k for k in range(nA) if indices[k] != 0]
        body = len(non_I)
        rng = (non_I[-1] - non_I[0]) if body >= 2 else (0 if body == 1 else -1)

        # Classify Pauli content
        pauli_types = tuple(sorted(set(PAULI_LABELS[indices[k]] for k in non_I)))

        # Is BW type?
        is_bw = False
        if body == 1 and indices[non_I[0]] == 1:  # single X
            is_bw = True
        elif body == 2 and rng == 1:
            p1, p2 = indices[non_I[0]], indices[non_I[1]]
            if p1 == 3 and p2 == 3:
                is_bw = True

        all_ops.append({
            'label': label,
            'coeff': float(c),
            'coeff_sq': float(c**2 * dimA),
            'range': rng,
            'body': body,
            'pauli_types': pauli_types,
            'is_bw': is_bw,
        })

    return all_ops


# Build H, get ground state
print(f"Sprint 097b: Non-BW operator structure at threshold")
print(f"q={q}, n={n} (periodic)", flush=True)

H = build_sq_potts_periodic(n, q, g)
evals_H, evecs_H = eigsh(H, k=1, which='SA')
psi = evecs_H[:, 0]
E0 = float(evals_H[0])
print(f"E0 = {E0:.8f}, time: {time.time()-t0:.1f}s\n", flush=True)


for nA in [5, 6]:
    dimA = 2**nA
    t1 = time.time()
    print(f"{'='*70}")
    print(f"nA={nA}", flush=True)

    rho_A = get_rho_A(psi, n, q, nA)
    H_E = entanglement_hamiltonian(rho_A)
    all_ops = full_pauli_analysis(H_E, nA)

    total_frob = sum(o['coeff_sq'] for o in all_ops)
    bw_ops = [o for o in all_ops if o['is_bw']]
    non_bw = [o for o in all_ops if not o['is_bw'] and abs(o['coeff']) > 1e-15]

    bw_frac = sum(o['coeff_sq'] for o in bw_ops) / total_frob
    non_bw_frac = 1 - bw_frac
    non_bw_total = sum(o['coeff_sq'] for o in non_bw)

    print(f"  BW fraction: {bw_frac:.6f}")
    print(f"  Non-BW fraction: {non_bw_frac:.6f} ({len(non_bw)} operators)")

    # Non-BW participation ratio
    if non_bw:
        nb_sq = np.array([o['coeff_sq'] for o in non_bw])
        nb_PR = float(np.sum(nb_sq)**2 / np.sum(nb_sq**2))
    else:
        nb_PR = 0

    print(f"  Non-BW participation ratio: {nb_PR:.1f}")

    # Non-BW cumulative (sorted)
    non_bw_sorted = sorted(non_bw, key=lambda x: x['coeff_sq'], reverse=True)
    if non_bw_sorted:
        nb_cum = np.cumsum([o['coeff_sq'] for o in non_bw_sorted]) / non_bw_total
        for target in [0.50, 0.90, 0.99]:
            idx = np.searchsorted(nb_cum, target)
            if idx < len(nb_cum):
                print(f"  Non-BW N({target:.0%}): {idx+1} operators")

    # Breakdown by (body, range)
    print(f"\n  Non-BW by (body, range):")
    br_frob = defaultdict(float)
    br_count = defaultdict(int)
    for o in non_bw:
        key = (o['body'], o['range'])
        br_frob[key] += o['coeff_sq']
        br_count[key] += 1

    print(f"  {'(body,range)':>14} {'count':>6} {'frac_of_nonBW':>14} {'frac_of_total':>14}")
    for key in sorted(br_frob.keys(), key=lambda k: br_frob[k], reverse=True):
        frac_nb = br_frob[key] / non_bw_total if non_bw_total > 0 else 0
        frac_tot = br_frob[key] / total_frob
        print(f"  {str(key):>14} {br_count[key]:>6} {frac_nb:>14.4f} {frac_tot:>14.6f}")

    # Breakdown by Pauli type tuple
    print(f"\n  Non-BW by Pauli content:")
    pt_frob = defaultdict(float)
    pt_count = defaultdict(int)
    for o in non_bw:
        key = ''.join(o['pauli_types'])
        pt_frob[key] += o['coeff_sq']
        pt_count[key] += 1

    print(f"  {'types':>12} {'count':>6} {'frac_of_nonBW':>14} {'frac_of_total':>14}")
    for key in sorted(pt_frob.keys(), key=lambda k: pt_frob[k], reverse=True)[:15]:
        frac_nb = pt_frob[key] / non_bw_total if non_bw_total > 0 else 0
        frac_tot = pt_frob[key] / total_frob
        print(f"  {key:>12} {pt_count[key]:>6} {frac_nb:>14.4f} {frac_tot:>14.6f}")

    # Top 15 non-BW operators
    print(f"\n  Top 15 non-BW operators:")
    for i, o in enumerate(non_bw_sorted[:15]):
        cum_nb = float(nb_cum[i]) if i < len(nb_cum) else 1.0
        print(f"    {i+1:>3}. {o['label']:>8} c={o['coeff']:+.6f} "
              f"frac_nb={o['coeff_sq']/non_bw_total:.4f} cum={cum_nb:.4f} "
              f"b={o['body']} r={o['range']} types={''.join(o['pauli_types'])}")

    dt = time.time() - t1
    print(f"  Time: {dt:.1f}s")

    entry = {
        'nA': nA,
        'bw_fraction': float(bw_frac),
        'non_bw_fraction': float(non_bw_frac),
        'n_non_bw': len(non_bw),
        'non_bw_PR': float(nb_PR),
        'body_range_breakdown': {
            str(k): {'count': br_count[k], 'frac_nb': br_frob[k]/non_bw_total if non_bw_total > 0 else 0,
                     'frac_total': br_frob[k]/total_frob}
            for k in br_frob
        },
        'pauli_type_breakdown': {
            k: {'count': pt_count[k], 'frac_nb': pt_frob[k]/non_bw_total if non_bw_total > 0 else 0,
                'frac_total': pt_frob[k]/total_frob}
            for k in pt_frob
        },
        'top_15_non_bw': [
            {'label': o['label'], 'coeff': o['coeff'], 'body': o['body'],
             'range': o['range'], 'types': ''.join(o['pauli_types'])}
            for o in non_bw_sorted[:15]
        ],
        'time': dt,
    }
    results['data'][f'nA{nA}'] = entry

    record(sprint=97, model='sq_potts', q=q, n=n,
           quantity='non_bw_PR', value=float(nb_PR),
           method=f'097b_nA{nA}')
    save()


save()
print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done!")
