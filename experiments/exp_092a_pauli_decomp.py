#!/usr/bin/env python3
"""Sprint 092a: Pauli decomposition of H_E for q=2 (Ising) at g_c=1/2.

Decompose H_E into Pauli strings on nA qubits. Identify which operators
dominate the BW residual (H_E - alpha*H_BW). Compare nA=3,4,5,6.

BW prediction: H_E ~ alpha * sum_i (distance_from_boundary) * h_i
where h_i are the local Hamiltonian terms restricted to subsystem A.
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from itertools import product
import json, time

t0 = time.time()

results = {
    'experiment': '092a_pauli_decomp',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_092a_pauli_decomp.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# Pauli matrices
I2 = np.eye(2)
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
sz = np.array([[1,0],[0,-1]])
paulis = {'I': I2, 'X': sx, 'Y': sy, 'Z': sz}

def pauli_string(ops):
    """Build tensor product of Pauli operators. ops is a string like 'XZII'."""
    result = np.array([[1.0+0j]])
    for c in ops:
        result = np.kron(result, paulis[c])
    return result

def decompose_pauli(H, nA):
    """Decompose Hermitian H into Pauli basis. Returns dict of coefficients."""
    dim = 2**nA
    labels = [''.join(p) for p in product('IXYZ', repeat=nA)]
    coeffs = {}
    for label in labels:
        P = pauli_string(label)
        c = np.real(np.trace(H @ P)) / dim
        if abs(c) > 1e-12:
            coeffs[label] = float(c)
    return coeffs

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
    return np.real(rho_A)

def entanglement_hamiltonian(rho_A):
    evals, evecs = np.linalg.eigh(rho_A)
    evals = np.clip(evals, 1e-30, None)
    H_E = -evecs @ np.diag(np.log(evals)) @ evecs.T
    return H_E

def classify_pauli(label):
    """Classify a Pauli string by body order, range, and type."""
    non_id = [(i, c) for i, c in enumerate(label) if c != 'I']
    body = len(non_id)
    if body == 0:
        return {'body': 0, 'range': 0, 'type': 'identity'}
    if body == 1:
        return {'body': 1, 'range': 0, 'type': f'1-body-{non_id[0][1]}'}
    positions = [p for p, _ in non_id]
    rng = max(positions) - min(positions)
    ops = ''.join(c for _, c in non_id)
    return {'body': body, 'range': rng, 'type': f'{body}-body-r{rng}-{ops}'}

q = 2
g = 0.5  # g_c = 1/q

for nA in [3, 4, 5, 6]:
    n = 2 * nA  # total system = 2x subsystem for half-chain cut
    print(f"\n=== q=2, n={n}, nA={nA} ===")
    t1 = time.time()

    H = build_sq_potts_periodic(n, q, g)
    evals, evecs = eigsh(H, k=1, which='SA')
    psi = evecs[:, 0]

    rho_A = get_rho_A(psi, n, q, nA)
    H_E = entanglement_hamiltonian(rho_A)

    # Pauli decomposition
    coeffs = decompose_pauli(H_E, nA)

    # Classify and sort
    classified = []
    total_weight = sum(c**2 for c in coeffs.values())
    for label, c in coeffs.items():
        info = classify_pauli(label)
        info['label'] = label
        info['coeff'] = c
        info['weight'] = c**2 / total_weight
        classified.append(info)

    classified.sort(key=lambda x: abs(x['coeff']), reverse=True)

    # Print top 20
    print(f"Total Pauli terms: {len(classified)}")
    print(f"Top 20 by |coeff|:")
    for i, item in enumerate(classified[:20]):
        print(f"  {i+1:2d}. {item['label']:8s} coeff={item['coeff']:+.6f}  "
              f"weight={item['weight']:.4%}  body={item['body']} range={item['range']}")

    # Aggregate by body order
    body_weights = {}
    for item in classified:
        b = item['body']
        body_weights[b] = body_weights.get(b, 0) + item['weight']

    print(f"\nWeight by body order:")
    for b in sorted(body_weights):
        print(f"  {b}-body: {body_weights[b]:.4%}")

    # Aggregate by range (for 2+ body)
    range_weights = {}
    for item in classified:
        if item['body'] >= 2:
            r = item['range']
            range_weights[r] = range_weights.get(r, 0) + item['weight']

    print(f"\nWeight by range (2+ body):")
    for r in sorted(range_weights):
        print(f"  range={r}: {range_weights[r]:.4%}")

    # BW terms: NN ZZ and single-site X (the Ising Hamiltonian terms)
    # For q=2 S_q Potts: H = -delta(s_i,s_j) - g*X = -(1+ZZ)/2 - g*X
    # BW: alpha * sum_site w(site) * [-(1+Z_i Z_{i+1})/2 - g*X_i]
    # where w(site) = (nA/pi) * sin(pi * dist / nA) for periodic subsystem
    bw_labels = set()
    for i in range(nA):
        # ZZ nearest-neighbor within A
        if i + 1 < nA:
            label = list('I' * nA)
            label[i] = 'Z'
            label[i+1] = 'Z'
            bw_labels.add(''.join(label))
        # X single-site
        label = list('I' * nA)
        label[i] = 'X'
        bw_labels.add(''.join(label))
        # Z single-site (from the diagonal part of -(1+ZZ)/2)
        label = list('I' * nA)
        label[i] = 'Z'
        bw_labels.add(''.join(label))

    bw_weight = sum(item['weight'] for item in classified if item['label'] in bw_labels)
    non_bw_weight = sum(item['weight'] for item in classified
                        if item['label'] not in bw_labels and item['body'] > 0)

    # What's outside BW?
    non_bw = [item for item in classified
              if item['label'] not in bw_labels and item['body'] > 0]
    non_bw.sort(key=lambda x: abs(x['coeff']), reverse=True)

    print(f"\nBW terms weight: {bw_weight:.4%}")
    print(f"Non-BW weight: {non_bw_weight:.4%}")
    print(f"Top 10 non-BW terms:")
    for i, item in enumerate(non_bw[:10]):
        print(f"  {i+1:2d}. {item['label']:8s} coeff={item['coeff']:+.6f}  "
              f"weight={item['weight']:.4%}  body={item['body']} range={item['range']}")

    # Identify dominant non-BW operator types
    non_bw_by_type = {}
    for item in non_bw:
        t = item['type']
        non_bw_by_type[t] = non_bw_by_type.get(t, 0) + item['weight']

    top_types = sorted(non_bw_by_type.items(), key=lambda x: x[1], reverse=True)[:10]
    print(f"\nTop non-BW operator types:")
    for t, w in top_types:
        print(f"  {t}: {w:.4%}")

    entry = {
        'n': n, 'nA': nA, 'q': q, 'g': g,
        'total_pauli_terms': len(classified),
        'top20': [{'label': x['label'], 'coeff': x['coeff'], 'weight': x['weight'],
                    'body': x['body'], 'range': x['range']} for x in classified[:20]],
        'body_weights': {str(k): v for k, v in body_weights.items()},
        'range_weights': {str(k): v for k, v in range_weights.items()},
        'bw_weight': bw_weight,
        'non_bw_weight': non_bw_weight,
        'non_bw_top10': [{'label': x['label'], 'coeff': x['coeff'], 'weight': x['weight'],
                          'body': x['body'], 'range': x['range']} for x in non_bw[:10]],
        'non_bw_by_type_top10': top_types,
        'time': time.time() - t1,
    }
    results['data'][f'nA={nA}'] = entry
    save()
    print(f"Time: {time.time()-t1:.1f}s")

print(f"\nTotal time: {time.time()-t0:.1f}s")
save()
print("Done!")
