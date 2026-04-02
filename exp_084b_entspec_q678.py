#!/usr/bin/env python3
"""Sprint 084b: Entanglement spectrum for q=6,7,8 at g_c=1/q (periodic S_q Potts).

These are the walking-broken cases. Compare spectrum structure to q=2,3,5 baseline.
Sizes: q=6 n=6,7; q=7 n=5,6,7; q=8 n=5,6
"""
import numpy as np
import json, time
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh

results = {
    'experiment': '084b_entspec_q678',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_084b_entspec_q678.json", "w") as f:
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
    """Get entanglement spectrum for bipartition into nA sites (A) and n-nA sites (B)."""
    dimA = q**nA
    dimB = q**(n - nA)
    psi_mat = psi.reshape(dimA, dimB)
    rho_A = psi_mat @ psi_mat.conj().T
    evals = np.linalg.eigvalsh(rho_A)
    evals = evals[::-1]
    evals = np.clip(evals, 1e-30, None)
    return evals

# Configuration
configs = {
    6: {'sizes': [6, 7]},
    7: {'sizes': [5, 6, 7]},
    8: {'sizes': [5, 6]},
}

# Known Re(c) values
Re_c = {6: 1.253, 7: 1.351, 8: 1.438}

print("Sprint 084b: Entanglement spectrum for S_q Potts (q=6,7,8)")
print("=" * 70, flush=True)

for q in [6, 7, 8]:
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

        psi = evecs[:, 0].real
        psi = psi / np.linalg.norm(psi)

        t0 = time.time()
        spec = get_entanglement_spectrum(psi, n, q, nA)
        t_spec = time.time() - t0

        xi = -np.log(spec)
        S_total = -np.sum(spec * np.log(spec))
        s_i = -spec * np.log(spec)
        ent_gap = xi[1] - xi[0] if len(xi) > 1 else 0
        n_sig = int(np.sum(s_i > 0.01 * S_total))
        cumfrac = np.cumsum(s_i) / S_total

        # Degeneracy detection: find groups of eigenvalues within 1% relative
        deg_groups = []
        i = 0
        while i < min(20, len(xi)):
            group = [xi[i]]
            j = i + 1
            while j < min(20, len(xi)) and abs(xi[j] - xi[i]) < 0.01 * max(abs(xi[i]), 0.1):
                group.append(xi[j])
                j += 1
            deg_groups.append(len(group))
            i = j

        print(f"  build={t_build:.1f}s, eig={t_eig:.1f}s, spec={t_spec:.1f}s")
        print(f"  S = {S_total:.6f}")
        print(f"  Entanglement gap: Δξ = {ent_gap:.4f}")
        print(f"  Significant eigenvalues: {n_sig} of {dimA}")
        print(f"  Top-1 λ = {spec[0]:.6f} (frac of S: {s_i[0]/S_total:.3f})")
        print(f"  Top-5 cumulative S fraction: {cumfrac[min(4,len(cumfrac)-1)]:.4f}")
        print(f"  Top-10 cumulative S fraction: {cumfrac[min(9,len(cumfrac)-1)]:.4f}")
        print(f"  Degeneracy pattern (first 20 levels): {deg_groups}")
        print(f"  Ent energies (first 15): {[f'{x:.3f}' for x in xi[:15]]}")

        size_data = {
            'n': n, 'nA': nA, 'dim': dim, 'dimA': dimA,
            'S_total': float(S_total),
            'ent_gap': float(ent_gap),
            'n_significant': n_sig,
            'top_eigenvalues': [float(x) for x in spec[:30]],
            'top_ent_energies': [float(x) for x in xi[:30]],
            'top_entropy_contributions': [float(x) for x in s_i[:30]],
            'cumulative_S_frac': [float(x) for x in cumfrac[:30]],
            'lambda_max': float(spec[0]),
            'degeneracy_pattern': deg_groups,
            'times': {'build': t_build, 'eig': t_eig, 'spec': t_spec},
        }
        q_data['sizes'][str(n)] = size_data

    results['data'][str(q)] = q_data
    save()

# Summary
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")
print(f"{'q':>3} {'n':>3} {'S':>8} {'Δξ':>7} {'λ_max':>8} {'n_sig':>6} {'top5%S':>7} {'degs':>20}")
for q in [6, 7, 8]:
    for n_str, d in results['data'][str(q)]['sizes'].items():
        top5 = d['cumulative_S_frac'][min(4, len(d['cumulative_S_frac'])-1)]
        degs = str(d['degeneracy_pattern'][:6])
        print(f"{q:3d} {d['n']:3d} {d['S_total']:8.4f} "
              f"{d['ent_gap']:7.4f} {d['lambda_max']:8.5f} {d['n_significant']:6d} {top5:7.4f} {degs}")

save()
print(f"\nSaved to results/sprint_084b_entspec_q678.json")

from db_utils import record
for q in [6, 7, 8]:
    for n_str, d in results['data'][str(q)]['sizes'].items():
        record(sprint=84, model='sq_potts', q=q, n=d['n'],
               quantity='ent_gap', value=d['ent_gap'],
               method='exact_diag_periodic',
               notes=f'half-chain, S={d["S_total"]:.4f}, n_sig={d["n_significant"]}')
print("Recorded to DB.")
