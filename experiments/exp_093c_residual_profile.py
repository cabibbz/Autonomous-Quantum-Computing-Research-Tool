#!/usr/bin/env python3
"""Sprint 093c: Spatial profile of BW residual within subsystem A.

After subtracting the best BW fit, where does the residual concentrate?
Is it at the entanglement boundary (edges of A) or uniformly spread?

Method: Compute H_res = H_E - alpha*H_BW - c*I, then measure
the "local Frobenius weight" at each position within A by projecting
onto all 2-body operators at bond (x, x+1).

Test at q=2,3,5 with nA=4,5 and q=2 nA=6.
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record
import json, time

t0 = time.time()

results = {
    'experiment': '093c_residual_profile',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_093c_residual_profile.json", "w") as f:
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


def bw_envelope_bond(N, nA, bond_x):
    return (N / np.pi) * np.sin(np.pi * (bond_x + 1) / N) * \
           np.sin(np.pi * (nA - bond_x - 1) / N) / np.sin(np.pi * nA / N)


def bw_envelope_site(N, nA, site_x):
    return (N / np.pi) * np.sin(np.pi * (site_x + 0.5) / N) * \
           np.sin(np.pi * (nA - site_x - 0.5) / N) / np.sin(np.pi * nA / N)


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


def local_weight_at_site(H, q, nA, site):
    """Frobenius weight of H restricted to operators acting non-trivially at 'site'.

    Partial trace approach: for each pair of local states at 'site',
    sum |H_{..a..;..b..}|^2 where a != b (off-diagonal in site index).
    Also include diagonal variations.

    Simpler: build projectors onto site subspace and measure.
    """
    dimA = q**nA
    # Reshape H as tensor with legs for each site
    # H[i1..iN, j1..jN] where i,j are q-valued site indices
    # Weight at site s = ||H - Tr_s(H) x I_s / q||^2
    # This is expensive. Use a simpler proxy: sum of |H_ij|^2 where
    # the local state at 'site' differs between i and j.

    total = 0.0
    acting_on_site = 0.0

    for i in range(dimA):
        ci = []
        tmp = i
        for _ in range(nA):
            ci.append(tmp % q)
            tmp //= q
        for j in range(dimA):
            cj = []
            tmp = j
            for _ in range(nA):
                cj.append(tmp % q)
                tmp //= q
            val = H[i, j]**2
            total += val
            if ci[site] != cj[site]:
                acting_on_site += val

    return acting_on_site, total


def local_weight_at_bond(H_res, q, nA, s1, s2):
    """Frobenius weight of residual that acts non-trivially at bond (s1, s2).

    Count |H_ij|^2 where local states at s1 or s2 differ.
    """
    dimA = q**nA
    acting = 0.0

    for i in range(dimA):
        ci = []
        tmp = i
        for _ in range(nA):
            ci.append(tmp % q)
            tmp //= q
        for j in range(i, dimA):  # use symmetry
            cj = []
            tmp = j
            for _ in range(nA):
                cj.append(tmp % q)
                tmp //= q

            if ci[s1] != cj[s1] or ci[s2] != cj[s2]:
                val = H_res[i, j]**2
                if i == j:
                    acting += val
                else:
                    acting += 2 * val  # symmetry factor

    return acting


# ============= Main experiment =============
configs = [
    (2, 4), (2, 5), (2, 6),
    (3, 4), (3, 5),
    (5, 4),
]

for q, nA in configs:
    n = 2 * nA
    g = 1.0 / q
    dimA = q**nA
    key = f"q{q}_nA{nA}"

    if q**n > 600000:
        print(f"SKIP q={q} nA={nA}: dim={q**n} too large")
        continue

    print(f"\n=== q={q}, n={n}, nA={nA}, dim_A={dimA} ===")
    t1 = time.time()

    H = build_sq_potts_periodic(n, q, g)
    evals, evecs = eigsh(H, k=1, which='SA')
    psi = evecs[:, 0]

    rho_A = get_rho_A(psi, n, q, nA)
    H_E = entanglement_hamiltonian(rho_A)

    H_BW = build_H_BW(q, n, nA, g)

    # Best fit
    dim = dimA
    t_mean = np.trace(H_E) / dim
    m_mean = np.trace(H_BW) / dim
    H_E_tl = H_E - t_mean * np.eye(dim)
    H_BW_tl = H_BW - m_mean * np.eye(dim)
    alpha = np.sum(H_E_tl * H_BW_tl) / np.sum(H_BW_tl * H_BW_tl)

    H_res = H_E_tl - alpha * H_BW_tl
    total_res = np.sum(H_res**2)
    total_HE = np.sum(H_E_tl**2)

    print(f"  alpha = {alpha:.4f}")
    print(f"  ||residual||²/||H_E||² = {total_res/total_HE:.6f}")

    # Measure where the residual concentrates: per-site weight
    # For efficiency, use trace of (H_res restricted to site) approach:
    # Weight at site s = sum_{i,j: i_s != j_s} |H_res[i,j]|^2

    # For small enough dimA, iterate directly
    if dimA <= 1300:
        site_weights = []
        for s in range(nA):
            w, _ = local_weight_at_site(H_res, q, nA, s)
            site_weights.append(w)

        # Normalize
        total_w = sum(site_weights)
        site_fracs = [w / total_w for w in site_weights]
        uniform_frac = 1.0 / nA

        print(f"  Residual site profile (fraction):")
        for s in range(nA):
            deviation = (site_fracs[s] - uniform_frac) / uniform_frac * 100
            pos = "boundary" if s == 0 or s == nA - 1 else "bulk"
            print(f"    site {s} ({pos}): {site_fracs[s]:.4f} (dev from uniform: {deviation:+.1f}%)")

        # Boundary vs bulk concentration
        boundary_w = site_fracs[0] + site_fracs[nA - 1]
        bulk_w = sum(site_fracs[1:nA-1])
        expected_boundary = 2.0 / nA
        expected_bulk = (nA - 2.0) / nA
        boundary_enrichment = (boundary_w / expected_boundary) if expected_boundary > 0 else 0

        print(f"  Boundary fraction: {boundary_w:.4f} (expected {expected_boundary:.4f}, enrichment={boundary_enrichment:.3f}x)")
        print(f"  Bulk fraction: {bulk_w:.4f} (expected {expected_bulk:.4f})")

        entry = {
            'q': q, 'n': n, 'nA': nA,
            'alpha': float(alpha),
            'residual_fraction': float(total_res / total_HE),
            'site_weights': [float(w) for w in site_weights],
            'site_fractions': [float(f) for f in site_fracs],
            'boundary_fraction': float(boundary_w),
            'boundary_enrichment': float(boundary_enrichment),
            'time': time.time() - t1,
        }
    else:
        # For large dimA, use a sampling approach or just report global residual
        print(f"  dimA={dimA} too large for per-site analysis, reporting global only")
        entry = {
            'q': q, 'n': n, 'nA': nA,
            'alpha': float(alpha),
            'residual_fraction': float(total_res / total_HE),
            'time': time.time() - t1,
        }

    results['data'][key] = entry

    # Record boundary enrichment to DB
    if 'boundary_enrichment' in entry:
        record(sprint=93, model='sq_potts', q=q, n=n,
               quantity='bw_boundary_enrichment', value=entry['boundary_enrichment'],
               method='residual_profile', notes=f'nA={nA}')

    save()
    print(f"  Time: {time.time()-t1:.1f}s")

print(f"\nTotal time: {time.time()-t0:.1f}s")
save()
print("Done!")
