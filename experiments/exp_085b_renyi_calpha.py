#!/usr/bin/env python3
"""Sprint 085b: Extract c_α from Rényi entropy using SIZE PAIRS.

CFT formula: S_α = (c/12)(1+1/α) ln(N/π) + const_α
Using pairs: c_α = 12 * (S_α(N2) - S_α(N1)) / ((1+1/α) * ln(N2/N1))

This cancels the non-universal additive constant.
Uses raw data from 085a and computes additional sizes where needed.
"""
import numpy as np
import json, time
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh

results = {
    'experiment': '085b_renyi_calpha',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_085b_renyi_calpha.json", "w") as f:
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

def get_entanglement_spectrum(psi, n, q, nA):
    dimA = q**nA
    dimB = q**(n - nA)
    psi_mat = psi.reshape(dimA, dimB)
    rho_A = psi_mat @ psi_mat.conj().T
    evals = np.linalg.eigvalsh(rho_A)
    evals = evals[::-1]
    evals = np.clip(evals, 1e-30, None)
    return evals

def renyi_entropy(spec, alpha):
    if alpha == 1:
        return -np.sum(spec * np.log(spec + 1e-300))
    elif alpha == np.inf:
        return -np.log(spec[0])
    else:
        return np.log(np.sum(spec**alpha)) / (1.0 - alpha)

def Re_c(q):
    if q <= 4:
        sqrt_Q = np.sqrt(q)
        p = np.pi / np.arccos(sqrt_Q / 2)
        return 1 - 6 / (p * (p - 1))
    else:
        alpha = np.arccosh(np.sqrt(q) / 2)
        return 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)

# Load 085a data
with open("results/sprint_085a_renyi_q28.json") as f:
    data_085a = json.load(f)

alphas_num = [0.5, 1, 2, 3, 5, 10, np.inf]
alpha_labels = ['0.5', '1', '2', '3', '5', '10', 'inf']

# Need additional sizes for better pairs. Compute missing ones.
# q=2: have 10,12,14 -> pairs (10,12), (12,14), (10,14)
# q=3: have 8,10 -> pair (8,10). Add n=6 for second pair.
# q=5: have 6,8 -> pair (6,8). Add n=4 for second pair.
# q=6: have 6,8 -> pair (6,8). Add n=4.
# q=7: have 6,7 -> pair (6,7). Add n=5.
# q=8: have 6,7 -> pair (6,7). Add n=5.

additional = {
    3: [6],   # q=3 n=6 dim=729
    5: [4],   # q=5 n=4 dim=625
    6: [4],   # q=6 n=4 dim=1296
    7: [4, 5],   # q=7 n=4 dim=2401, n=5 dim=16807
    8: [4, 5],   # q=8 n=4 dim=4096, n=5 dim=32768
}

# Compute additional sizes
all_S = {}  # (q, n, alpha_label) -> S_alpha

# First load existing data from 085a
for q_str, qd in data_085a['data'].items():
    q = int(q_str)
    for n_str, nd in qd['sizes'].items():
        n = int(n_str)
        for label in alpha_labels:
            all_S[(q, n, label)] = nd['S_alpha'][label]

print("Sprint 085b: c_α from size pairs (canceling additive constant)")
print("=" * 70)

# Compute additional sizes
for q, ns in additional.items():
    gc = 1.0 / q
    for n in ns:
        nA = n // 2
        dim = q**n
        print(f"Computing q={q} n={n} (dim={dim:,}) ... ", end='', flush=True)
        t0 = time.time()
        H = build_sq_potts_periodic(n, q, gc)
        evals, evecs = eigsh(H, k=1, which='SA')
        psi = evecs[:, 0].real
        psi = psi / np.linalg.norm(psi)
        spec = get_entanglement_spectrum(psi, n, q, nA)
        t = time.time() - t0
        print(f"{t:.1f}s")
        for a, label in zip(alphas_num, alpha_labels):
            S = renyi_entropy(spec, a)
            all_S[(q, n, label)] = float(S)

# Now compute c_α from size pairs
# For each q, use (n1, n2) pair where n2 > n1
size_pairs = {
    2: [(10, 12), (12, 14), (10, 14)],
    3: [(6, 8), (8, 10), (6, 10)],
    5: [(4, 6), (6, 8), (4, 8)],
    6: [(4, 6), (6, 8), (4, 8)],
    7: [(4, 5), (5, 6), (6, 7), (4, 7)],
    8: [(4, 5), (5, 6), (6, 7), (4, 7)],
}

print(f"\n{'='*70}")
print("c_α from size pairs")
print(f"{'='*70}")

for q in [2, 3, 5, 6, 7, 8]:
    rc = Re_c(q)
    print(f"\nq={q}, Re(c)={rc:.4f}")
    q_data = {'q': q, 'Re_c': float(rc), 'pairs': {}}

    for n1, n2 in size_pairs[q]:
        pair_key = f"{n1}_{n2}"
        pair_data = {'n1': n1, 'n2': n2, 'c_alpha': {}, 'c_over_Rec': {}}

        print(f"  Pair ({n1},{n2}):", end='')
        for a, label in zip(alphas_num, alpha_labels):
            S1 = all_S[(q, n1, label)]
            S2 = all_S[(q, n2, label)]
            if a == np.inf:
                prefactor = 1.0 / 12.0
            else:
                prefactor = (1.0 + 1.0/a) / 12.0
            c_alpha = (S2 - S1) / (prefactor * np.log(n2 / n1))
            pair_data['c_alpha'][label] = float(c_alpha)
            pair_data['c_over_Rec'][label] = float(c_alpha / rc)

        q_data['pairs'][pair_key] = pair_data

        # Print the widest pair (most reliable)
        row = ""
        for label in alpha_labels:
            row += f"  α={label}:{pair_data['c_over_Rec'][label]:.3f}"
        print(row)

    results['data'][str(q)] = q_data

save()

# Best pair summary (widest baseline for each q)
best_pairs = {2: '10_14', 3: '6_10', 5: '4_8', 6: '4_8', 7: '4_7', 8: '4_7'}
print(f"\n{'='*70}")
print("SUMMARY: c_α/Re(c) from best (widest) size pair")
print(f"{'='*70}")
header = f"{'q':>3} {'pair':>7}"
for label in alpha_labels:
    header += f" {'α='+label:>8}"
print(header)

for q in [2, 3, 5, 6, 7, 8]:
    bp = best_pairs[q]
    d = results['data'][str(q)]['pairs'][bp]
    n1, n2 = d['n1'], d['n2']
    row = f"{q:3d} ({n1},{n2:2d})"
    for label in alpha_labels:
        row += f" {d['c_over_Rec'][label]:8.4f}"
    print(row)

# Key question: does c_∞/Re(c) → 1 at large q?
print(f"\n{'='*70}")
print("KEY: c_α/Re(c) convergence across α")
print(f"{'='*70}")
print("If walking breakdown is entropy-only, c_∞/Re(c) should be CLOSER to 1 than c_1/Re(c)")
print()
for q in [2, 3, 5, 6, 7, 8]:
    bp = best_pairs[q]
    d = results['data'][str(q)]['pairs'][bp]
    c1 = d['c_over_Rec']['1']
    cinf = d['c_over_Rec']['inf']
    print(f"  q={q}: c_1/Re(c) = {c1:.4f}, c_∞/Re(c) = {cinf:.4f}, "
          f"Δ = {abs(cinf - 1) - abs(c1 - 1):+.4f} ({'∞ closer' if abs(cinf-1) < abs(c1-1) else '1 closer'})")

save()
print(f"\nSaved to results/sprint_085b_renyi_calpha.json")

from db_utils import record
for q in [2, 3, 5, 6, 7, 8]:
    bp = best_pairs[q]
    d = results['data'][str(q)]['pairs'][bp]
    for label in alpha_labels:
        record(sprint=85, model='sq_potts', q=q,
               quantity=f'c_renyi_{label}', value=d['c_alpha'][label],
               method='renyi_size_pair',
               notes=f'pair=({d["n1"]},{d["n2"]}), c/Re(c)={d["c_over_Rec"][label]:.4f}')
print("Recorded to DB.")
