#!/usr/bin/env python3
"""Sprint 097a: Full Pauli decomposition of H_E for q=2 n=14, nA=3-6.

Decompose H_E into the complete Pauli string basis. Measure:
- Cumulative R² (sorted by |coeff|²): how many operators for 90%, 99%, 99.9%?
- Operator participation ratio PR = (sum |c|²)² / sum |c|⁴
- Compare: BW uses O(nA) operators. Is H_E compact (O(nA²)) or diffuse (O(4^nA))?

Key insight: Pauli strings are orthogonal, so greedy selection = sorted decomposition.
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
    'experiment': '097a_pauli_compactness',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'n': n, 'g_c': g,
    'data': {},
}

def save():
    with open("results/sprint_097a_pauli_compactness.json", "w") as f:
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


# Pauli matrices
I2 = np.eye(2)
X2 = np.array([[0, 1], [1, 0]], dtype=float)
Y2 = np.array([[0, -1], [1, 0]], dtype=float)
Z2 = np.array([[1, 0], [0, -1]], dtype=float)
PAULIS = [I2, X2, Y2, Z2]
PAULI_LABELS = ['I', 'X', 'Y', 'Z']


def pauli_decompose(H_E, nA):
    """Full Pauli decomposition of H_E for nA qubits.

    Returns list of (label, coefficient, range, body_order, pauli_type).
    Pauli type: 'D' (diagonal: I,Z only), 'S' (shift: I,X only), 'M' (mixed).
    """
    dimA = 2**nA
    coeffs = []

    for indices in product(range(4), repeat=nA):
        # Build label
        label = ''.join(PAULI_LABELS[i] for i in indices)
        if all(i == 0 for i in indices):
            continue  # skip identity

        # Build operator via Kronecker product
        op = PAULIS[indices[0]]
        for k in range(1, nA):
            op = np.kron(op, PAULIS[indices[k]])

        # Coefficient: Tr(P · H_E) / dimA
        c = np.sum(op * H_E.T) / dimA  # efficient trace via elementwise

        if abs(c) < 1e-15:
            continue

        # Classify
        non_I = [k for k in range(nA) if indices[k] != 0]
        body_order = len(non_I)
        rng = (non_I[-1] - non_I[0]) if body_order >= 2 else (0 if body_order == 1 else -1)

        # Pauli type classification
        types_present = set(PAULI_LABELS[indices[k]] for k in non_I)
        if types_present <= {'Z'}:
            ptype = 'D'  # diagonal (Potts-like: delta bonds)
        elif types_present <= {'X'}:
            ptype = 'S'  # shift (field-like)
        elif types_present <= {'X', 'Z'}:
            ptype = 'DxS'  # mixed diagonal × shift
        elif 'Y' in types_present and types_present <= {'Y'}:
            ptype = 'M'  # pure Y (momentum-like)
        elif 'Y' in types_present:
            ptype = 'M'  # mixed with Y
        else:
            ptype = 'other'

        # Is this a BW-type operator? BW has: NN ZZ bonds + single-site X fields
        is_bw = False
        if body_order == 1 and indices[non_I[0]] == 1:  # single X
            is_bw = True
        elif body_order == 2 and rng == 1:
            p1, p2 = indices[non_I[0]], indices[non_I[1]]
            if p1 == 3 and p2 == 3:  # ZZ (Potts bond for q=2)
                is_bw = True

        coeffs.append({
            'label': label,
            'coeff': float(c),
            'coeff_sq': float(c**2 * dimA),  # Frobenius contribution
            'range': rng,
            'body': body_order,
            'type': ptype,
            'is_bw': is_bw,
        })

    return coeffs


# Build H, get ground state
print(f"Sprint 097a: Full Pauli decomposition of H_E")
print(f"q={q}, n={n} (periodic)", flush=True)

H = build_sq_potts_periodic(n, q, g)
evals_H, evecs_H = eigsh(H, k=1, which='SA')
psi = evecs_H[:, 0]
E0 = float(evals_H[0])
print(f"E0 = {E0:.8f}, H build+diag: {time.time()-t0:.1f}s\n", flush=True)


nA_list = [3, 4, 5, 6]

for nA in nA_list:
    dimA = 2**nA
    n_ops_total = 4**nA - 1  # total traceless Pauli strings
    ratio = nA / n
    t1 = time.time()

    print(f"{'='*70}")
    print(f"nA={nA}, dimA={dimA}, total Pauli ops={n_ops_total}, nA/n={ratio:.3f}", flush=True)

    rho_A = get_rho_A(psi, n, q, nA)
    H_E = entanglement_hamiltonian(rho_A)

    coeffs = pauli_decompose(H_E, nA)

    # Sort by Frobenius contribution (descending)
    coeffs.sort(key=lambda x: x['coeff_sq'], reverse=True)

    total_frob = sum(c['coeff_sq'] for c in coeffs)

    # Cumulative R²
    cum_frob = np.cumsum([c['coeff_sq'] for c in coeffs])
    cum_R2 = cum_frob / total_frob

    # Find milestones
    milestones = {}
    for target in [0.90, 0.95, 0.99, 0.999, 0.9999]:
        idx = np.searchsorted(cum_R2, target)
        if idx < len(cum_R2):
            milestones[f'{target:.4f}'] = int(idx + 1)
        else:
            milestones[f'{target:.4f}'] = len(cum_R2)

    # Participation ratio
    c_sq = np.array([c['coeff_sq'] for c in coeffs])
    PR = float(np.sum(c_sq)**2 / np.sum(c_sq**2)) if np.sum(c_sq**2) > 0 else 0

    # BW operator count and fraction
    n_bw = sum(1 for c in coeffs if c['is_bw'])
    bw_frob = sum(c['coeff_sq'] for c in coeffs if c['is_bw'])
    bw_frac = bw_frob / total_frob

    # Type breakdown
    type_frob = {}
    type_count = {}
    for c in coeffs:
        t = c['type']
        type_frob[t] = type_frob.get(t, 0) + c['coeff_sq']
        type_count[t] = type_count.get(t, 0) + 1
    type_frac = {t: v / total_frob for t, v in type_frob.items()}

    dt = time.time() - t1

    # Print results
    print(f"  Non-zero operators: {len(coeffs)} / {n_ops_total}")
    print(f"  Participation ratio: {PR:.1f}")
    print(f"  BW operators: {n_bw}, BW fraction: {bw_frac:.6f}")
    print(f"  Milestones (operators for R²):")
    for target, count in milestones.items():
        print(f"    R²≥{target}: {count} operators ({count/n_ops_total*100:.1f}% of basis)")

    print(f"  Type breakdown (fraction of ||H_E||²):")
    for t in sorted(type_frac.keys(), key=lambda x: type_frac[x], reverse=True):
        print(f"    {t:>5}: {type_frac[t]:.6f} ({type_count[t]} ops)")

    # Top 10 operators
    print(f"  Top 10 operators:")
    for i, c in enumerate(coeffs[:10]):
        bw_tag = " [BW]" if c['is_bw'] else ""
        cum = cum_R2[i]
        print(f"    {i+1}. {c['label']:>8} c={c['coeff']:+.6f} "
              f"frac={c['coeff_sq']/total_frob:.4f} cum={cum:.4f} "
              f"r={c['range']} b={c['body']} {c['type']}{bw_tag}")

    print(f"  Time: {dt:.1f}s")

    entry = {
        'nA': nA, 'dimA': dimA,
        'n_nonzero': len(coeffs),
        'n_total_basis': n_ops_total,
        'participation_ratio': float(PR),
        'n_bw_ops': n_bw,
        'bw_fraction': float(bw_frac),
        'milestones': milestones,
        'type_fraction': {t: float(v) for t, v in type_frac.items()},
        'type_count': type_count,
        'top_20': [
            {'label': c['label'], 'coeff': c['coeff'],
             'frac': c['coeff_sq']/total_frob, 'range': c['range'],
             'body': c['body'], 'type': c['type'], 'is_bw': c['is_bw']}
            for c in coeffs[:20]
        ],
        'time': dt,
    }
    results['data'][f'nA{nA}'] = entry

    record(sprint=97, model='sq_potts', q=q, n=n,
           quantity='HE_participation_ratio', value=float(PR),
           method=f'097a_nA{nA}')
    record(sprint=97, model='sq_potts', q=q, n=n,
           quantity='HE_n99pct', value=float(milestones.get('0.9900', 0)),
           method=f'097a_nA{nA}')
    save()


# Summary
print(f"\n\n{'='*70}")
print("SUMMARY: Operator compactness of H_E (q=2 n=14)")
print(f"{'nA':>4} {'dimA':>6} {'basis':>7} {'nonzero':>8} {'PR':>8} "
      f"{'N(90%)':>8} {'N(99%)':>8} {'N(99.9%)':>10} {'BW frac':>10}")
print("-" * 85)
for nA in nA_list:
    d = results['data'][f'nA{nA}']
    print(f"{nA:>4} {d['dimA']:>6} {d['n_total_basis']:>7} {d['n_nonzero']:>8} "
          f"{d['participation_ratio']:>8.1f} "
          f"{d['milestones'].get('0.9000', '—'):>8} "
          f"{d['milestones'].get('0.9900', '—'):>8} "
          f"{d['milestones'].get('0.9990', '—'):>10} "
          f"{d['bw_fraction']:>10.6f}")

save()
print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done!")
