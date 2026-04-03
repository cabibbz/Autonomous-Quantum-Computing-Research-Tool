#!/usr/bin/env python3
"""Sprint 085a: Rényi entropies S_α for q=2-8 at g_c=1/q (periodic S_q Potts).

Compute ground state, get entanglement spectrum, then compute S_α for multiple α.
Extract effective c_α from CFT formula: S_α = (c/12)(1+1/α) ln[(N/π)sin(πℓ/N)] + const
For half-chain: S_α = (c_α/12)(1+1/α) ln(N/π) + const_α

Sizes: q=2 n=10,12,14; q=3 n=8,10; q=5 n=6,8; q=6 n=6,8; q=7 n=6,7; q=8 n=6,7
"""
import numpy as np
import json, time
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh

results = {
    'experiment': '085a_renyi_q28',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_085a_renyi_q28.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def build_sq_potts_periodic(n, q, g):
    """Build S_q Potts Hamiltonian on periodic chain."""
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
    """Get entanglement spectrum for half-chain bipartition."""
    dimA = q**nA
    dimB = q**(n - nA)
    psi_mat = psi.reshape(dimA, dimB)
    rho_A = psi_mat @ psi_mat.conj().T
    evals = np.linalg.eigvalsh(rho_A)
    evals = evals[::-1]  # descending
    evals = np.clip(evals, 1e-30, None)
    return evals

def renyi_entropy(spec, alpha):
    """Compute Rényi entropy S_α = (1/(1-α)) ln(Σ λ^α).
    For α=1 returns von Neumann entropy. α='inf' returns min-entropy."""
    if alpha == 1:
        return -np.sum(spec * np.log(spec + 1e-300))
    elif alpha == np.inf or alpha == 'inf':
        return -np.log(spec[0])
    else:
        return np.log(np.sum(spec**alpha)) / (1.0 - alpha)

def c_from_renyi(S_alpha, alpha, n):
    """Extract c_α from Rényi entropy using CFT formula.
    S_α = (c/12)(1 + 1/α) ln(N/π) + const
    For pairs of sizes: c_α = 12(S2-S1) / ((1+1/α)(ln(n2/π) - ln(n1/π)))
    For single size: c_α = S_α / ((1+1/α)/12 * ln(n/π))
    """
    if alpha == np.inf or alpha == 'inf':
        prefactor = 1.0 / 12.0  # (1+1/∞)/12 = 1/12
    else:
        prefactor = (1.0 + 1.0/alpha) / 12.0
    return S_alpha / (prefactor * np.log(n / np.pi))

# Configuration
configs = {
    2: [10, 12, 14],
    3: [8, 10],
    5: [6, 8],
    6: [6, 8],
    7: [6, 7],
    8: [6, 7],
}

# Coulomb gas Re(c)
def Re_c(q):
    if q <= 4:
        sqrt_Q = np.sqrt(q)
        p = np.pi / np.arccos(sqrt_Q / 2)
        return 1 - 6 / (p * (p - 1))
    else:
        alpha = np.arccosh(np.sqrt(q) / 2)
        return 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)

Re_c_vals = {q: Re_c(q) for q in [2, 3, 4, 5, 6, 7, 8]}

alphas = [0.5, 1, 2, 3, 5, 10, 'inf']
alpha_labels = ['0.5', '1', '2', '3', '5', '10', 'inf']

print("Sprint 085a: Rényi entropies for S_q Potts (q=2-8)")
print("=" * 70, flush=True)

for q in [2, 3, 5, 6, 7, 8]:
    gc = 1.0 / q
    sizes = configs[q]
    rc = Re_c_vals[q]
    print(f"\n{'='*70}")
    print(f"q={q}, g_c={gc:.4f}, Re(c)={rc:.4f}")
    print(f"{'='*70}", flush=True)

    q_data = {'q': q, 'gc': gc, 'Re_c': float(rc), 'sizes': {}}

    for n in sizes:
        nA = n // 2
        dim = q**n
        print(f"\n  n={n}, dim={dim:,} ... ", end='', flush=True)

        t0 = time.time()
        H = build_sq_potts_periodic(n, q, gc)
        evals, evecs = eigsh(H, k=1, which='SA')
        psi = evecs[:, 0].real
        psi = psi / np.linalg.norm(psi)
        spec = get_entanglement_spectrum(psi, n, q, nA)
        t_total = time.time() - t0
        print(f"{t_total:.1f}s")

        # Compute Rényi entropies
        S_vals = {}
        c_vals = {}
        for a, label in zip(alphas, alpha_labels):
            a_num = np.inf if a == 'inf' else a
            S = renyi_entropy(spec, a_num)
            c_a = c_from_renyi(S, a_num, n)
            S_vals[label] = float(S)
            c_vals[label] = float(c_a)

        # Print results
        print(f"  {'α':>5} {'S_α':>10} {'c_α':>10} {'c_α/Re(c)':>10}")
        for label in alpha_labels:
            ratio = c_vals[label] / rc
            print(f"  {label:>5} {S_vals[label]:10.6f} {c_vals[label]:10.4f} {ratio:10.4f}")

        q_data['sizes'][str(n)] = {
            'n': n, 'nA': nA, 'dim': dim,
            'S_alpha': S_vals,
            'c_alpha': c_vals,
            'c_alpha_over_Rec': {label: c_vals[label]/rc for label in alpha_labels},
            'spectrum_top5': [float(x) for x in spec[:5]],
            'time': t_total,
        }

    results['data'][str(q)] = q_data
    save()

# Summary: c_α/Re(c) for all q at largest size, for each α
print(f"\n{'='*70}")
print("SUMMARY: c_α / Re(c) at largest size for each q")
print(f"{'='*70}")
header = f"{'q':>3} {'n':>3}"
for label in alpha_labels:
    header += f" {'α='+label:>8}"
print(header)
for q in [2, 3, 5, 6, 7, 8]:
    sizes = list(results['data'][str(q)]['sizes'].keys())
    n_max = max(int(s) for s in sizes)
    d = results['data'][str(q)]['sizes'][str(n_max)]
    row = f"{q:3d} {n_max:3d}"
    for label in alpha_labels:
        row += f" {d['c_alpha_over_Rec'][label]:8.4f}"
    print(row)

save()
print(f"\nSaved to results/sprint_085a_renyi_q28.json")

from db_utils import record
for q in [2, 3, 5, 6, 7, 8]:
    sizes = list(results['data'][str(q)]['sizes'].keys())
    n_max = max(int(s) for s in sizes)
    d = results['data'][str(q)]['sizes'][str(n_max)]
    for label in alpha_labels:
        record(sprint=85, model='sq_potts', q=q, n=int(n_max),
               quantity=f'c_alpha_{label}', value=d['c_alpha'][label],
               method='renyi_periodic',
               notes=f'S_{label}={d["S_alpha"][label]:.6f}, c/Re(c)={d["c_alpha_over_Rec"][label]:.4f}')
print("Recorded to DB.")
