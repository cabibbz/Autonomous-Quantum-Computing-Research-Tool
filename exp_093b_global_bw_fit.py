#!/usr/bin/env python3
"""Sprint 093b: Global BW Hamiltonian fit across q=2-7.

Construct H_BW = sum_x beta(x) * h_x where h_x is the full local Hamiltonian
density (coupling + field) at bond/site x, with BW envelope beta(x).
Find best alpha: H_E = alpha * H_BW + c*I.
Compare: (1) BW with single alpha (coupling=field share same beta)
         (2) BW with separate alpha for coupling vs field
Quantify residual across q and check if deviations correlate with walking.

Fixed nA=4 for all q=2,3,4,5,6,7 to compare cleanly.
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record
import json, time

t0 = time.time()

results = {
    'experiment': '093b_global_bw_fit',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_093b_global_bw_fit.json", "w") as f:
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
    """delta(s_i, s_j) on nA-site Hilbert space."""
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
    """sum_{k=1}^{q-1} X^k at given site."""
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
    """BW envelope at bond position x (between site x and x+1)."""
    return (N / np.pi) * np.sin(np.pi * (bond_x + 1) / N) * \
           np.sin(np.pi * (nA - bond_x - 1) / N) / np.sin(np.pi * nA / N)


def bw_envelope_site(N, nA, site_x):
    """BW envelope at site position x. Use midpoint convention."""
    return (N / np.pi) * np.sin(np.pi * (site_x + 0.5) / N) * \
           np.sin(np.pi * (nA - site_x - 0.5) / N) / np.sin(np.pi * nA / N)


def build_H_BW(q, n, nA, g):
    """Build the BW Hamiltonian with proper envelope.
    H_BW = sum_bond beta(bond) * [-delta(si,sj)] + sum_site beta(site) * [-g * field]
    Using bond-centered envelope for coupling, site-centered for field.
    For sites at boundary of A: field acts only within A, so we include it.
    """
    dimA = q**nA
    H_BW = np.zeros((dimA, dimA))

    # Bond terms: -delta(si, si+1) * beta(bond)
    for x in range(nA - 1):
        beta = bw_envelope_bond(n, nA, x)
        delta = potts_delta_bond(q, nA, x, x + 1)
        H_BW -= beta * delta

    # Field terms: -g * field(site) * beta(site)
    for x in range(nA):
        beta = bw_envelope_site(n, nA, x)
        field = clock_field_site(q, nA, x)
        H_BW -= g * beta * field

    return H_BW


def build_H_BW_2param(q, n, nA, g):
    """Build separate coupling and field parts for 2-parameter fit.
    Returns (H_coup, H_field) each with BW envelope baked in.
    """
    dimA = q**nA
    H_coup = np.zeros((dimA, dimA))
    H_field = np.zeros((dimA, dimA))

    for x in range(nA - 1):
        beta = bw_envelope_bond(n, nA, x)
        delta = potts_delta_bond(q, nA, x, x + 1)
        H_coup -= beta * delta

    for x in range(nA):
        beta = bw_envelope_site(n, nA, x)
        field = clock_field_site(q, nA, x)
        H_field -= beta * field

    return H_coup, H_field


def frobenius_R2(H_target, H_model):
    """R² in Frobenius norm. Remove identity component from both first."""
    dim = H_target.shape[0]
    t_mean = np.trace(H_target) / dim
    m_mean = np.trace(H_model) / dim
    H_t = H_target - t_mean * np.eye(dim)
    H_m = H_model - m_mean * np.eye(dim)

    # Best scalar fit: alpha = Tr(H_t @ H_m) / Tr(H_m @ H_m)
    alpha = np.sum(H_t * H_m) / np.sum(H_m * H_m)
    residual = H_t - alpha * H_m
    ss_res = np.sum(residual**2)
    ss_tot = np.sum(H_t**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return R2, alpha


# ============= Main experiment =============
nA = 4

for q in [2, 3, 4, 5, 6, 7]:
    n = 2 * nA
    g = 1.0 / q
    dimA = q**nA
    key = f"q{q}"
    print(f"\n=== q={q}, n={n}, nA={nA}, g_c={g:.4f}, dim={q**n} ===")
    t1 = time.time()

    # Skip if Hilbert space too large
    if q**n > 600000:
        print(f"  SKIP: dim={q**n} too large for exact diag")
        continue

    H = build_sq_potts_periodic(n, q, g)
    evals, evecs = eigsh(H, k=1, which='SA')
    psi = evecs[:, 0]
    E0 = float(evals[0])

    rho_A = get_rho_A(psi, n, q, nA)
    H_E = entanglement_hamiltonian(rho_A)

    # Method 1: Single-parameter BW fit (coupling+field share same alpha)
    H_BW = build_H_BW(q, n, nA, g)
    R2_1param, alpha_1param = frobenius_R2(H_E, H_BW)
    print(f"  1-param BW: alpha={alpha_1param:.6f}, R²={R2_1param:.8f}")

    # Method 2: Two-parameter fit (separate alpha for coupling and field)
    H_coup, H_field = build_H_BW_2param(q, n, nA, g)
    # Least squares: H_E_traceless = a1*H_coup + a2*H_field
    dim = dimA
    H_E_tl = H_E - (np.trace(H_E) / dim) * np.eye(dim)
    H_c_flat = H_coup.flatten()
    H_f_flat = H_field.flatten()
    A_mat = np.column_stack([H_c_flat, H_f_flat])
    params, _, _, _ = np.linalg.lstsq(A_mat, H_E_tl.flatten(), rcond=None)
    a_coup, a_field = params[0], params[1]
    H_fit_2p = a_coup * H_coup + a_field * H_field
    H_fit_2p_tl = H_fit_2p - (np.trace(H_fit_2p) / dim) * np.eye(dim)
    ss_res = np.sum((H_E_tl - H_fit_2p_tl)**2)
    ss_tot = np.sum(H_E_tl**2)
    R2_2param = 1 - ss_res / ss_tot
    ratio_cf = a_field / a_coup if a_coup != 0 else float('nan')
    print(f"  2-param BW: a_coup={a_coup:.6f}, a_field={a_field:.6f}, "
          f"ratio={ratio_cf:.6f}, R²={R2_2param:.8f}")
    print(f"  Expected ratio (g_c) = {g:.4f}")
    print(f"  Actual / expected = {ratio_cf/g:.4f}")

    # Method 3: Per-bond and per-site coefficients (unconstrained)
    # Extract individual bond coefficients
    bond_coeffs = []
    for x in range(nA - 1):
        delta = potts_delta_bond(q, nA, x, x + 1)
        # Project: Tr(H_E @ delta) / Tr(delta @ delta)
        c = np.sum(H_E * delta) / np.sum(delta * delta)
        bond_coeffs.append(float(c))

    field_coeffs = []
    for x in range(nA):
        field = clock_field_site(q, nA, x)
        c = np.sum(H_E * field) / np.sum(field * field)
        field_coeffs.append(float(c))

    # BW envelope values
    bw_bonds = [bw_envelope_bond(n, nA, x) for x in range(nA - 1)]
    bw_sites = [bw_envelope_site(n, nA, x) for x in range(nA)]

    # Normalized profile comparison
    bond_norm = np.array(bond_coeffs) / np.mean(bond_coeffs)
    bw_bond_norm = np.array(bw_bonds) / np.mean(bw_bonds)
    bond_profile_dev = np.max(np.abs(bond_norm - bw_bond_norm))

    field_norm = np.array(field_coeffs) / np.mean(field_coeffs)
    bw_site_norm = np.array(bw_sites) / np.mean(bw_sites)
    field_profile_dev = np.max(np.abs(field_norm - bw_site_norm))

    print(f"  Bond profile: {bond_coeffs}")
    print(f"  BW bond env:  {[f'{x:.3f}' for x in bw_bonds]}")
    print(f"  Bond shape deviation: {bond_profile_dev:.6f}")
    print(f"  Field profile: {field_coeffs}")
    print(f"  BW site env:   {[f'{x:.3f}' for x in bw_sites]}")
    print(f"  Field shape deviation: {field_profile_dev:.6f}")

    # Frobenius norm of residual (what's NOT BW)
    H_residual_1p = H_E_tl - alpha_1param * (H_BW - (np.trace(H_BW)/dim)*np.eye(dim))
    frac_residual_1p = np.sum(H_residual_1p**2) / np.sum(H_E_tl**2)

    H_residual_2p = H_E_tl - H_fit_2p_tl
    frac_residual_2p = np.sum(H_residual_2p**2) / np.sum(H_E_tl**2)

    print(f"  Residual fraction: 1-param={frac_residual_1p:.6f}, 2-param={frac_residual_2p:.6f}")
    print(f"  Improvement from 2-param: {(frac_residual_1p - frac_residual_2p)/frac_residual_1p*100:.1f}%")

    entry = {
        'q': q, 'n': n, 'nA': nA, 'g_c': g, 'E0': E0,
        'R2_1param': float(R2_1param),
        'alpha_1param': float(alpha_1param),
        'R2_2param': float(R2_2param),
        'a_coup': float(a_coup),
        'a_field': float(a_field),
        'ratio_field_coup': float(ratio_cf),
        'ratio_over_gc': float(ratio_cf / g),
        'bond_coeffs': bond_coeffs,
        'field_coeffs': field_coeffs,
        'bw_bond_envelope': bw_bonds,
        'bw_site_envelope': bw_sites,
        'bond_profile_deviation': float(bond_profile_dev),
        'field_profile_deviation': float(field_profile_dev),
        'residual_1param': float(frac_residual_1p),
        'residual_2param': float(frac_residual_2p),
        'time': time.time() - t1,
    }
    results['data'][key] = entry

    # Record to DB
    record(sprint=93, model='sq_potts', q=q, n=n,
           quantity='bw_R2_1param', value=R2_1param, method='global_bw_fit')
    record(sprint=93, model='sq_potts', q=q, n=n,
           quantity='bw_R2_2param', value=R2_2param, method='global_bw_fit')
    record(sprint=93, model='sq_potts', q=q, n=n,
           quantity='bw_alpha', value=alpha_1param, method='global_bw_fit')
    record(sprint=93, model='sq_potts', q=q, n=n,
           quantity='bw_field_coup_ratio', value=ratio_cf, method='global_bw_fit')

    save()
    print(f"  Time: {time.time()-t1:.1f}s")

print(f"\nTotal time: {time.time()-t0:.1f}s")
save()
print("Done!")
