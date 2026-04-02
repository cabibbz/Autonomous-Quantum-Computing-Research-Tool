#!/usr/bin/env python3
"""Sprint 093a: Position-dependent BW entanglement temperature profile.

For periodic S_q Potts at g_c=1/q, extract H_E = -log(rho_A) and project
onto local Potts operators at each bond/site within subsystem A.
Compare the position-dependent coefficients to BW envelope prediction:
  beta_BW(x) ~ (N/pi) * sin(pi*(x+1)/N) * sin(pi*(nA-x)/N) / sin(pi*nA/N)

Test q=2 (nA=3,4,5,6), q=3 (nA=3,4), q=5 (nA=3,4).
"""
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from gpu_utils import eigsh
import json, time

t0 = time.time()

results = {
    'experiment': '093a_bw_envelope',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_093a_bw_envelope.json", "w") as f:
        json.dump(results, f, indent=2, default=str)


def build_sq_potts_periodic(n, q, g):
    """S_q Potts: H = -sum delta(s_i,s_{i+1}) - g * sum_{k=1}^{q-1} X^k."""
    dim = q**n
    H = lil_matrix((dim, dim), dtype=complex)
    for idx in range(dim):
        c = []
        tmp = idx
        for _ in range(n):
            c.append(tmp % q)
            tmp //= q
        # Potts coupling
        diag = 0.0
        for site in range(n):
            nxt = (site + 1) % n
            if c[site] == c[nxt]:
                diag -= 1.0
        H[idx, idx] += diag
        # S_q transverse field: sum_{k=1}^{q-1} X^k at each site
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
    """Reduced density matrix of first nA sites (contiguous)."""
    dimA = q**nA
    dimB = q**(n - nA)
    psi_mat = psi.reshape(dimA, dimB)
    rho_A = psi_mat @ psi_mat.conj().T
    return np.array(rho_A)


def entanglement_hamiltonian(rho_A):
    """H_E = -log(rho_A)."""
    evals, evecs = np.linalg.eigh(rho_A)
    evals = np.clip(evals, 1e-30, None)
    H_E = -evecs @ np.diag(np.log(evals)) @ evecs.conj().T
    return np.real(H_E)


def potts_delta_operator(q, nA, site_i, site_j):
    """Build delta(s_i, s_j) operator on nA-site Hilbert space (dim q^nA)."""
    dimA = q**nA
    op = np.zeros((dimA, dimA))
    for idx in range(dimA):
        c = []
        tmp = idx
        for _ in range(nA):
            c.append(tmp % q)
            tmp //= q
        if c[site_i] == c[site_j]:
            op[idx, idx] = 1.0
    return op


def clock_field_operator(q, nA, site):
    """Build sum_{k=1}^{q-1} X^k operator at given site on nA-site space."""
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


def project_onto_operator(H_E, op):
    """Project H_E onto operator: coeff = Tr(H_E @ op) / Tr(op @ op)."""
    return np.trace(H_E @ op) / np.trace(op @ op)


def bw_envelope(N, nA, bond_x):
    """BW entanglement temperature for bond x (0-indexed) within A.

    For periodic chain of N sites, subsystem A = sites [0, nA-1]:
    beta(x) ~ (N/pi) * sin(pi*(x+1)/N) * sin(pi*(nA-x-1)/N) / sin(pi*nA/N)

    bond_x = 0 means bond between site 0 and site 1 (near left boundary of A).
    """
    return (N / np.pi) * np.sin(np.pi * (bond_x + 1) / N) * \
           np.sin(np.pi * (nA - bond_x - 1) / N) / np.sin(np.pi * nA / N)


def bw_envelope_site(N, nA, site_x):
    """BW entanglement temperature for site x (0-indexed) within A.

    For field term at site x, average of adjacent bond positions.
    Actually: beta(x) at site position x+0.5.
    Use midpoint: sin(pi*(x+0.5)/N) * sin(pi*(nA-x-0.5)/N) / sin(pi*nA/N)
    """
    return (N / np.pi) * np.sin(np.pi * (site_x + 0.5) / N) * \
           np.sin(np.pi * (nA - site_x - 0.5) / N) / np.sin(np.pi * nA / N)


# ============= Main experiment =============
configs = [
    # (q, nA)
    (2, 3), (2, 4), (2, 5), (2, 6),
    (3, 3), (3, 4),
    (5, 3), (5, 4),
]

for q, nA in configs:
    n = 2 * nA  # half-chain cut
    g = 1.0 / q
    key = f"q{q}_nA{nA}"
    print(f"\n=== q={q}, n={n}, nA={nA}, g_c={g:.4f} ===")
    t1 = time.time()

    # Build H, get ground state
    H = build_sq_potts_periodic(n, q, g)
    evals, evecs = eigsh(H, k=1, which='SA')
    psi = evecs[:, 0]
    E0 = evals[0]
    print(f"  E0 = {E0:.6f}")

    # Get rho_A and H_E
    rho_A = get_rho_A(psi, n, q, nA)
    H_E = entanglement_hamiltonian(rho_A)

    # Project onto Potts coupling at each bond within A
    bond_coeffs = []
    for x in range(nA - 1):
        delta_op = potts_delta_operator(q, nA, x, x + 1)
        coeff = project_onto_operator(H_E, delta_op)
        bond_coeffs.append(float(np.real(coeff)))

    # Project onto field operator at each site within A
    field_coeffs = []
    for x in range(nA):
        field_op = clock_field_operator(q, nA, x)
        coeff = project_onto_operator(H_E, field_op)
        field_coeffs.append(float(np.real(coeff)))

    # BW envelope prediction (unnormalized)
    bw_bond = [bw_envelope(n, nA, x) for x in range(nA - 1)]
    bw_site = [bw_envelope_site(n, nA, x) for x in range(nA)]

    # Fit: actual = alpha * BW_envelope + const
    # For bonds: fit alpha from least squares
    if len(bond_coeffs) >= 2:
        bw_arr = np.array(bw_bond)
        actual_arr = np.array(bond_coeffs)
        # alpha = (actual . bw) / (bw . bw) if we assume zero offset
        # With offset: [alpha, c] = argmin |actual - alpha*bw - c|^2
        A_mat = np.column_stack([bw_arr, np.ones(len(bw_arr))])
        params_bond, res_bond, _, _ = np.linalg.lstsq(A_mat, actual_arr, rcond=None)
        alpha_bond = params_bond[0]
        offset_bond = params_bond[1]
        fit_bond = alpha_bond * bw_arr + offset_bond
        ss_res = np.sum((actual_arr - fit_bond)**2)
        ss_tot = np.sum((actual_arr - np.mean(actual_arr))**2)
        R2_bond = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    else:
        alpha_bond = bond_coeffs[0] / bw_bond[0] if bw_bond[0] != 0 else 0
        offset_bond = 0
        R2_bond = 1.0
        fit_bond = [alpha_bond * bw_bond[0]]

    # For sites
    if len(field_coeffs) >= 2:
        bw_s_arr = np.array(bw_site)
        actual_s_arr = np.array(field_coeffs)
        A_mat_s = np.column_stack([bw_s_arr, np.ones(len(bw_s_arr))])
        params_site, res_site, _, _ = np.linalg.lstsq(A_mat_s, actual_s_arr, rcond=None)
        alpha_site = params_site[0]
        offset_site = params_site[1]
        fit_site = alpha_site * bw_s_arr + offset_site
        ss_res_s = np.sum((actual_s_arr - fit_site)**2)
        ss_tot_s = np.sum((actual_s_arr - np.mean(actual_s_arr))**2)
        R2_site = 1 - ss_res_s / ss_tot_s if ss_tot_s > 0 else 1.0
    else:
        alpha_site = field_coeffs[0] / bw_site[0] if bw_site[0] != 0 else 0
        offset_site = 0
        R2_site = 1.0
        fit_site = [alpha_site * bw_site[0]]

    # Print results
    print(f"  Bond coefficients (Potts coupling delta):")
    for x in range(nA - 1):
        dev = (bond_coeffs[x] - fit_bond[x]) / abs(bond_coeffs[x]) * 100 if bond_coeffs[x] != 0 else 0
        print(f"    bond {x}-{x+1}: actual={bond_coeffs[x]:.6f}  "
              f"BW_fit={fit_bond[x]:.6f}  dev={dev:.2f}%")
    print(f"  Bond: alpha={alpha_bond:.6f}, offset={offset_bond:.6f}, R²={R2_bond:.6f}")

    print(f"  Field coefficients (clock field):")
    for x in range(nA):
        dev = (field_coeffs[x] - fit_site[x]) / abs(field_coeffs[x]) * 100 if field_coeffs[x] != 0 else 0
        print(f"    site {x}: actual={field_coeffs[x]:.6f}  "
              f"BW_fit={fit_site[x]:.6f}  dev={dev:.2f}%")
    print(f"  Field: alpha={alpha_site:.6f}, offset={offset_site:.6f}, R²={R2_site:.6f}")

    # Ratio of field to bond alphas (should be g_c for BW)
    if alpha_bond != 0:
        alpha_ratio = alpha_site / alpha_bond
        print(f"  alpha_field/alpha_bond = {alpha_ratio:.6f} (expected g_c={g:.4f})")

    # Max deviation from BW envelope
    if len(bond_coeffs) >= 2:
        max_dev_bond = np.max(np.abs(np.array(bond_coeffs) - np.array(fit_bond))) / np.max(np.abs(bond_coeffs))
    else:
        max_dev_bond = 0

    entry = {
        'q': q, 'n': n, 'nA': nA, 'g_c': g, 'E0': float(E0),
        'bond_coeffs': bond_coeffs,
        'field_coeffs': field_coeffs,
        'bw_bond_envelope': bw_bond,
        'bw_site_envelope': bw_site,
        'bond_alpha': float(alpha_bond),
        'bond_offset': float(offset_bond),
        'bond_R2': float(R2_bond),
        'site_alpha': float(alpha_site),
        'site_offset': float(offset_site),
        'site_R2': float(R2_site),
        'alpha_ratio': float(alpha_site / alpha_bond) if alpha_bond != 0 else None,
        'max_dev_bond': float(max_dev_bond),
        'time': time.time() - t1,
    }
    results['data'][key] = entry
    save()
    print(f"  Time: {time.time()-t1:.1f}s")

print(f"\nTotal time: {time.time()-t0:.1f}s")
save()
print("Done!")
