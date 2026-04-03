#!/usr/bin/env python3
"""Sprint 092b: Generalized operator decomposition of H_E for q=2,3,4,5.

Uses clock-shift (Weyl-Heisenberg) basis at FIXED nA=3 for clean comparison.
nA=3 keeps basis size manageable: q=5 gives 25^3 = 15625 basis elements.

Classify operators by body order, range, and physical type (D/F/M).
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from itertools import product as iproduct
import json, time

t0 = time.time()

results = {
    'experiment': '092b_generalized_decomp',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_092b_generalized_decomp.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def clock_shift_basis(q):
    """Generate q^2 Weyl-Heisenberg operators for q-level system."""
    omega = np.exp(2j * np.pi / q)
    X = np.zeros((q, q), dtype=complex)
    for s in range(q):
        X[(s + 1) % q, s] = 1.0
    Z = np.diag([omega**s for s in range(q)])
    basis = {}
    Id = np.eye(q, dtype=complex)
    Xpow = [Id.copy()]
    for a in range(1, q):
        Xpow.append(Xpow[-1] @ X)
    Zpow = [Id.copy()]
    for b in range(1, q):
        Zpow.append(Zpow[-1] @ Z)
    for a in range(q):
        for b in range(q):
            label = f"X{a}Z{b}"
            basis[label] = Xpow[a] @ Zpow[b]
    return basis

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

def classify_single_op(label, q):
    """Classify a single-site operator X^a Z^b."""
    a = int(label[1])
    b = int(label[3])
    if a == 0 and b == 0:
        return 'I'  # identity
    if a == 0:
        return 'D'  # density (clock/Z-type)
    if b == 0:
        return 'F'  # field (shift/X-type)
    return 'M'      # mixed (XZ-type)

def is_bw_term(labels, q):
    """Check if multi-site operator is a BW term for S_q Potts."""
    non_id = [(i, l) for i, l in enumerate(labels) if l != "X0Z0"]
    body = len(non_id)
    if body == 0:
        return True
    if body == 1:
        _, l = non_id[0]
        a, b = int(l[1]), int(l[3])
        # Field (X^k) or single-site density (Z^b)
        return True  # all 1-body operators are "BW-like"
    if body == 2:
        (p1, l1), (p2, l2) = non_id
        if abs(p2 - p1) != 1:
            return False  # not NN
        a1, b1 = int(l1[1]), int(l1[3])
        a2, b2 = int(l2[1]), int(l2[3])
        # NN density-density (Z^b ⊗ Z^{-b}) from Potts coupling
        if a1 == 0 and a2 == 0:
            return True
        return False
    return False

nA = 3
for q in [2, 3, 4, 5]:
    g = 1.0 / q
    n = 2 * nA  # total = 6
    print(f"\n{'='*60}")
    print(f"=== q={q}, n={n}, nA={nA}, g_c={g:.4f} ===")
    t1 = time.time()

    dim_A = q**nA
    total_basis = (q**2)**nA
    print(f"  dim_A={dim_A}, basis_size={total_basis}")

    H = build_sq_potts_periodic(n, q, g)
    evals, evecs = eigsh(H, k=1, which='SA')
    psi = evecs[:, 0]
    print(f"  E0 = {evals[0]:.6f}, time so far: {time.time()-t1:.1f}s")

    rho_A = get_rho_A(psi, n, q, nA)
    H_E = entanglement_hamiltonian(rho_A)

    basis = clock_shift_basis(q)
    single_labels = list(basis.keys())

    # Full decomposition
    coeffs = {}
    for multi_label in iproduct(single_labels, repeat=nA):
        op = np.array([[1.0+0j]])
        for label in multi_label:
            op = np.kron(op, basis[label])
        c = np.trace(H_E @ op.conj().T) / dim_A
        if abs(c) > 1e-10:
            coeffs[multi_label] = complex(c)

    print(f"  Non-zero coeffs: {len(coeffs)}, time: {time.time()-t1:.1f}s")

    # Classify
    classified = []
    total_weight = sum(abs(c)**2 for c in coeffs.values())
    for labels, c in coeffs.items():
        non_id = [(i, l) for i, l in enumerate(labels) if l != "X0Z0"]
        body = len(non_id)
        positions = [p for p, _ in non_id] if non_id else [0]
        rng = max(positions) - min(positions) if body > 0 else 0
        type_str = ''.join(classify_single_op(l, q) for _, l in non_id) if non_id else 'I'
        is_bw = is_bw_term(labels, q)

        classified.append({
            'labels': labels,
            'coeff_abs': abs(c),
            'weight': abs(c)**2 / total_weight,
            'body': body,
            'range': rng,
            'type': type_str,
            'is_bw': is_bw,
        })

    classified.sort(key=lambda x: x['coeff_abs'], reverse=True)

    # Print top 15
    print(f"\nTop 15 by |coeff|:")
    for i, item in enumerate(classified[:15]):
        bw = "BW" if item['is_bw'] else "  "
        print(f"  {i+1:2d}. [{bw}] |c|={item['coeff_abs']:.6f}  "
              f"w={item['weight']:.4%}  body={item['body']} r={item['range']} "
              f"type={item['type']}")

    # Body-order weights
    body_weights = {}
    for item in classified:
        b = item['body']
        body_weights[b] = body_weights.get(b, 0) + item['weight']
    print(f"\nWeight by body order:")
    for b in sorted(body_weights):
        print(f"  {b}-body: {body_weights[b]:.4%}")

    # BW vs non-BW
    bw_weight = sum(item['weight'] for item in classified if item['is_bw'])
    non_bw = [item for item in classified if not item['is_bw'] and item['body'] > 0]
    non_bw_weight = sum(item['weight'] for item in non_bw)
    print(f"\nBW weight: {bw_weight:.4%}")
    print(f"Non-BW weight: {non_bw_weight:.4%}")

    # Non-BW by type
    type_weights = {}
    for item in non_bw:
        t = item['type']
        type_weights[t] = type_weights.get(t, 0) + item['weight']
    top_types = sorted(type_weights.items(), key=lambda x: x[1], reverse=True)[:10]
    print(f"\nNon-BW by type (D=density, F=field, M=mixed):")
    for t, w in top_types:
        print(f"  {t}: {w:.4%}")

    # Non-BW by body+range
    br_weights = {}
    for item in non_bw:
        key = f"{item['body']}-body r={item['range']}"
        br_weights[key] = br_weights.get(key, 0) + item['weight']
    print(f"\nNon-BW by body+range:")
    for k, w in sorted(br_weights.items(), key=lambda x: x[1], reverse=True)[:8]:
        print(f"  {k}: {w:.4%}")

    entry = {
        'q': q, 'n': n, 'nA': nA, 'g': g,
        'total_nonzero': len(classified),
        'body_weights': {str(k): v for k, v in body_weights.items()},
        'bw_weight': bw_weight,
        'non_bw_weight': non_bw_weight,
        'non_bw_by_type': top_types,
        'non_bw_by_body_range': sorted(br_weights.items(), key=lambda x: x[1], reverse=True)[:8],
        'time': time.time() - t1,
    }
    results['data'][f'q={q}'] = entry
    save()
    print(f"\nTime: {time.time()-t1:.1f}s")

print(f"\nTotal time: {time.time()-t0:.1f}s")
save()
print("Done!")
