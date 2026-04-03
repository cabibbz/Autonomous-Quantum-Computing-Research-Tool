#!/usr/bin/env python3
"""Sprint 084a: Entanglement spectrum for q=2,3,5 at g_c=1/q (periodic S_q Potts).

Compute ground state, trace out half chain, get full eigenvalue spectrum of rho_A.
Extract: entanglement energies xi_i = -ln(lambda_i), entanglement gap, per-eigenvalue
entropy contributions s_i = -lambda_i * ln(lambda_i).

Sizes: q=2 n=10,12,14; q=3 n=8,10; q=5 n=6,8
"""
import numpy as np
import json, time
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh

results = {
    'experiment': '084a_entspec_q235',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_084a_entspec_q235.json", "w") as f:
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
        # Diagonal: coupling
        diag = 0.0
        for site in range(n):
            nxt = (site + 1) % n
            if c[site] == c[nxt]:
                diag -= 1.0
        H[idx, idx] += diag
        # Off-diagonal: transverse field (S_q: all cyclic shifts)
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
    """Get entanglement spectrum for bipartition into nA sites (A) and n-nA sites (B).
    Returns eigenvalues of rho_A sorted descending."""
    dimA = q**nA
    dimB = q**(n - nA)
    # Reshape state vector into (dimA, dimB) matrix
    psi_mat = psi.reshape(dimA, dimB)
    # rho_A = Tr_B(|psi><psi|) = psi_mat @ psi_mat^dagger
    rho_A = psi_mat @ psi_mat.conj().T
    # Eigenvalues
    evals = np.linalg.eigvalsh(rho_A)
    evals = evals[::-1]  # descending
    # Clip tiny negatives
    evals = np.clip(evals, 1e-30, None)
    return evals

# Configuration
configs = {
    2: {'sizes': [10, 12, 14], 'nA_frac': 0.5},
    3: {'sizes': [8, 10], 'nA_frac': 0.5},
    5: {'sizes': [6, 8], 'nA_frac': 0.5},
}

# Known Re(c) values
Re_c = {2: 0.500, 3: 0.800, 5: 1.138}

print("Sprint 084a: Entanglement spectrum for S_q Potts (q=2,3,5)")
print("=" * 70, flush=True)

for q in [2, 3, 5]:
    gc = 1.0 / q
    sizes = configs[q]['sizes']
    print(f"\n{'='*70}")
    print(f"q={q}, g_c={gc:.4f}, Re(c)={Re_c[q]}")
    print(f"{'='*70}", flush=True)

    q_data = {'q': q, 'gc': gc, 'Re_c': Re_c[q], 'sizes': {}}

    for n in sizes:
        nA = n // 2
        dim = q**n
        dimA = q**nA
        print(f"\n  n={n}, nA={nA}, dim={dim:,}, dimA={dimA:,} ... ", flush=True)

        t0 = time.time()
        H = build_sq_potts_periodic(n, q, gc)
        t_build = time.time() - t0

        t0 = time.time()
        evals, evecs = eigsh(H, k=1, which='SA')
        t_eig = time.time() - t0

        psi = evecs[:, 0].real  # Ground state (real for this Hamiltonian)
        # Normalize
        psi = psi / np.linalg.norm(psi)

        t0 = time.time()
        spec = get_entanglement_spectrum(psi, n, q, nA)
        t_spec = time.time() - t0

        # Entanglement energies
        xi = -np.log(spec)
        # Entropy
        S_total = -np.sum(spec * np.log(spec))
        # Per-eigenvalue entropy contribution
        s_i = -spec * np.log(spec)
        # Entanglement gap
        ent_gap = xi[1] - xi[0] if len(xi) > 1 else 0
        # Number of significant eigenvalues (contribute >1% of entropy)
        n_sig = int(np.sum(s_i > 0.01 * S_total))
        # Fraction of entropy from top k eigenvalues
        cumfrac = np.cumsum(s_i) / S_total

        # CFT prediction: S = (c/6) * ln(N/pi * sin(pi*nA/N)) + const
        # For half chain: S = (c/6) * ln(N/pi) + const
        c_eff = S_total / (np.log(n / np.pi) / 6) if n > 2 else 0

        print(f"  build={t_build:.1f}s, eig={t_eig:.1f}s, spec={t_spec:.1f}s")
        print(f"  S = {S_total:.6f}, c_eff = {c_eff:.4f} (Re(c)={Re_c[q]})")
        print(f"  Entanglement gap: Δξ = {ent_gap:.4f}")
        print(f"  Significant eigenvalues: {n_sig} of {dimA}")
        print(f"  Top-1 λ = {spec[0]:.6f} (frac of S: {s_i[0]/S_total:.3f})")
        print(f"  Top-5 cumulative S fraction: {cumfrac[min(4,len(cumfrac)-1)]:.4f}")
        print(f"  Top-10 cumulative S fraction: {cumfrac[min(9,len(cumfrac)-1)]:.4f}")
        print(f"  Entanglement energies (first 10): {[f'{x:.3f}' for x in xi[:10]]}")

        size_data = {
            'n': n, 'nA': nA, 'dim': dim, 'dimA': dimA,
            'S_total': float(S_total), 'c_eff': float(c_eff),
            'ent_gap': float(ent_gap),
            'n_significant': n_sig,
            'top_eigenvalues': [float(x) for x in spec[:20]],
            'top_ent_energies': [float(x) for x in xi[:20]],
            'top_entropy_contributions': [float(x) for x in s_i[:20]],
            'cumulative_S_frac': [float(x) for x in cumfrac[:20]],
            'lambda_max': float(spec[0]),
            'times': {'build': t_build, 'eig': t_eig, 'spec': t_spec},
        }
        q_data['sizes'][str(n)] = size_data

    results['data'][str(q)] = q_data
    save()

# Summary table
print(f"\n{'='*70}")
print("SUMMARY: Entanglement spectrum characteristics")
print(f"{'='*70}")
print(f"{'q':>3} {'n':>3} {'S':>8} {'c_eff':>7} {'Δξ':>7} {'λ_max':>8} {'n_sig':>6} {'top5%S':>7}")
for q in [2, 3, 5]:
    for n_str, d in results['data'][str(q)]['sizes'].items():
        top5 = d['cumulative_S_frac'][min(4, len(d['cumulative_S_frac'])-1)]
        print(f"{q:3d} {d['n']:3d} {d['S_total']:8.4f} {d['c_eff']:7.4f} "
              f"{d['ent_gap']:7.4f} {d['lambda_max']:8.5f} {d['n_significant']:6d} {top5:7.4f}")

save()
print(f"\nSaved to results/sprint_084a_entspec_q235.json")

from db_utils import record
for q in [2, 3, 5]:
    for n_str, d in results['data'][str(q)]['sizes'].items():
        record(sprint=84, model='sq_potts', q=q, n=d['n'],
               quantity='ent_gap', value=d['ent_gap'],
               method='exact_diag_periodic',
               notes=f'half-chain bipartition, S={d["S_total"]:.4f}')
        record(sprint=84, model='sq_potts', q=q, n=d['n'],
               quantity='c_eff', value=d['c_eff'],
               method='entanglement_spectrum',
               notes=f'n_sig={d["n_significant"]}, lambda_max={d["lambda_max"]:.5f}')
print("Recorded to DB.")
