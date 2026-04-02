#!/usr/bin/env python3
"""Sprint 075b: Entanglement entropy on q=2 cylinders at criticality.

Measure half-chain entanglement entropy S(Lx/2) on cylinders at g_c(Ly)
for Ly=2,3,4,5. Compare with 1D (Ly=1) entropy.

For 1D CFT: S = (c/3)*ln(n) + const. For cylinder: S should scale with
the cut area (Ly sites) plus logarithmic corrections from the 1D direction.

Key question: How does c_eff on cylinders interpolate between 1D (c=0.5) and 2D?

Method: Full diagonalization for small dim, eigsh for larger. Compute ground
state, trace out left half of cylinder, compute von Neumann entropy.
"""
import numpy as np
import json
import time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from gpu_utils import eigsh


def build_H_cylinder(Lx, Ly, q, g):
    """Build hybrid Hamiltonian on Lx x Ly cylinder (open x, periodic y)."""
    n = Lx * Ly
    dim = q ** n
    potts_2site = np.zeros(q ** 2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q ** 2, q ** 2), format='csr')
    X = np.zeros((q, q))
    for s in range(q):
        X[(s + 1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)
    H = csr_matrix((dim, dim))
    bonds = set()
    for y in range(Ly):
        for x in range(Lx):
            site = y * Lx + x
            if x + 1 < Lx:
                nbr = y * Lx + (x + 1)
                bonds.add((min(site, nbr), max(site, nbr)))
            nbr_v = ((y + 1) % Ly) * Lx + x
            if nbr_v != site:
                bonds.add((min(site, nbr_v), max(site, nbr_v)))
    for (i, j) in bonds:
        if j == i + 1:
            left = q ** i
            right = q ** (n - j - 1)
            op = potts_op
            if left > 1:
                op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1:
                op = sp_kron(op, sp_eye(right), format='csr')
        else:
            indices = np.arange(dim)
            si = (indices // (q ** i)) % q
            sj = (indices // (q ** j)) % q
            diag_vals = (si == sj).astype(float)
            op = diags(diag_vals, 0, shape=(dim, dim), format='csr')
        H = H - op
    for i in range(n):
        left = q ** i
        right = q ** (n - i - 1)
        op = XpXd.copy()
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H


def get_ground_state(Lx, Ly, q, g):
    """Get ground state vector."""
    H = build_H_cylinder(Lx, Ly, q, g)
    dim = H.shape[0]
    if dim < 2000:
        evals, evecs = np.linalg.eigh(H.toarray())
        return evecs[:, 0], evals[0]
    evals, evecs = eigsh(H, k=2, which='SA')
    idx = np.argsort(evals)
    return evecs[:, idx[0]], evals[idx[0]]


def entanglement_entropy(psi, q, n, subsys_sites):
    """Compute von Neumann entropy for a subsystem of specified sites.

    Sites are labeled 0..n-1. State is in the product basis |s_0 s_1 ... s_{n-1}>.
    Reshape into (dim_A, dim_B) and compute SVD.

    For cylinder with site ordering (y*Lx + x), a vertical cut at x=l
    means subsystem A = all sites with x < l.
    """
    n_A = len(subsys_sites)
    n_B = n - n_A
    dim_A = q ** n_A
    dim_B = q ** n_B

    # Reorder sites: A sites first, then B sites
    all_sites = list(range(n))
    B_sites = [s for s in all_sites if s not in subsys_sites]

    # Build permuted state
    # psi is indexed by sum_i s_i * q^i
    # We need to reshape into (dim_A, dim_B)
    psi_matrix = np.zeros((dim_A, dim_B))

    for idx in range(len(psi)):
        # Decode index into site values
        site_vals = np.zeros(n, dtype=int)
        tmp = idx
        for s in range(n):
            site_vals[s] = tmp % q
            tmp //= q

        # Compute A and B indices
        a_idx = 0
        for ki, s in enumerate(subsys_sites):
            a_idx += site_vals[s] * (q ** ki)
        b_idx = 0
        for ki, s in enumerate(B_sites):
            b_idx += site_vals[s] * (q ** ki)

        psi_matrix[a_idx, b_idx] += psi[idx]

    # SVD
    sv = np.linalg.svd(psi_matrix, compute_uv=False)
    sv = sv[sv > 1e-14]
    sv2 = sv ** 2
    return -np.sum(sv2 * np.log(sv2))


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 075b: Entanglement entropy on q=2 cylinders", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '075b_cylinder_entropy',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': 2,
}

# Critical points for each Ly
gc_values = {
    1: 0.250,   # 1D chain (exact)
    2: 0.451,   # Ly=2 cylinder
    3: 0.655,   # Ly=3 cylinder
    4: 0.688,   # Ly=4 cylinder
    5: 0.701,   # Ly=5 cylinder (from 075a)
}

# For each Ly, compute S at g_c for several Lx values
# Subsystem A = left half of cylinder (sites with x < Lx/2)
entropy_data = {}

for Ly in [1, 2, 3, 4]:
    gc = gc_values[Ly]
    print(f"\n--- Ly={Ly}, g_c={gc:.3f} ---", flush=True)

    # Choose Lx values (even, so half-cut is clean)
    if Ly == 1:
        Lx_list = [4, 6, 8, 10, 12]
    elif Ly == 2:
        Lx_list = [4, 6, 8]
    elif Ly == 3:
        Lx_list = [4, 6]
    elif Ly == 4:
        Lx_list = [4]

    ly_data = []
    for Lx in Lx_list:
        n = Lx * Ly
        dim = 2 ** n
        if dim > 2_000_000:
            print(f"  Lx={Lx} (n={n}, dim={dim:,}) — skipping (too large for eigenvector)", flush=True)
            continue

        t0 = time.time()
        psi, E0 = get_ground_state(Lx, Ly, 2, gc)
        dt_gs = time.time() - t0

        # Subsystem A: left half (x < Lx/2)
        half = Lx // 2
        A_sites = []
        for y in range(Ly):
            for x in range(half):
                A_sites.append(y * Lx + x)

        t1 = time.time()
        S = entanglement_entropy(psi, 2, n, A_sites)
        dt_ent = time.time() - t1

        n_A = len(A_sites)
        print(f"  Lx={Lx} (n={n}, dim={dim:,}): S={S:.4f}, n_A={n_A}, {dt_gs:.1f}s GS + {dt_ent:.1f}s entropy", flush=True)

        ly_data.append({
            'Lx': Lx, 'n': n, 'dim': dim,
            'S': float(S), 'n_A': n_A,
            'E0': float(E0),
            'time_gs': float(dt_gs), 'time_ent': float(dt_ent)
        })

    entropy_data[Ly] = ly_data

results['entropy_data'] = {str(k): v for k, v in entropy_data.items()}

# Analysis: extract c_eff for each Ly
# For 1D CFT: S = (c/3)*ln((2n/pi)*sin(pi*l/n)) + const
# For open BC (cylinder x-direction): S = (c/6)*ln(n_x) + f(Ly)
# With open x-BC, use S = (c_eff/6)*ln((2Lx/pi)*sin(pi/2)) + const = (c_eff/6)*ln(2Lx/pi) + const
print(f"\n{'='*60}", flush=True)
print("Analysis: c_eff extraction from S vs ln(Lx)", flush=True)
print(f"{'='*60}", flush=True)

c_eff_results = {}
for Ly in sorted(entropy_data.keys()):
    data = entropy_data[Ly]
    if len(data) < 2:
        print(f"  Ly={Ly}: Only {len(data)} point(s) — cannot extract c_eff", flush=True)
        if len(data) == 1:
            c_eff_results[Ly] = {'S': data[0]['S'], 'Lx': data[0]['Lx'], 'c_eff': None}
        continue

    Lx_arr = np.array([d['Lx'] for d in data])
    S_arr = np.array([d['S'] for d in data])

    # For open BC: S = (c/6)*ln(Lx) + const  (half-chain cut on open chain)
    # Actually for periodic x: S = (c/3)*ln(Lx) + const
    # Our cylinder has OPEN x, so S ~ (c_eff/6)*ln(Lx) + const
    # But we also have Ly cuts (area law contribution)
    # S = (c_eff/6)*ln(Lx) + alpha*Ly + const
    # With fixed Ly, fit S = a*ln(Lx) + b, then c_eff = 6*a

    ln_Lx = np.log(Lx_arr)
    # Linear fit: S = a*ln(Lx) + b
    if len(data) >= 2:
        A_mat = np.vstack([ln_Lx, np.ones_like(ln_Lx)]).T
        result_fit = np.linalg.lstsq(A_mat, S_arr, rcond=None)
        a, b = result_fit[0]
        c_eff = 6 * a  # Open BC: coefficient is c/6

        # Pairwise c_eff from consecutive pairs
        c_pairs = []
        for i in range(len(data) - 1):
            dS = S_arr[i+1] - S_arr[i]
            dln = ln_Lx[i+1] - ln_Lx[i]
            c_pairs.append(6 * dS / dln)

        print(f"  Ly={Ly}: c_eff = {c_eff:.4f} (from fit), pairwise: {[f'{c:.3f}' for c in c_pairs]}", flush=True)
        c_eff_results[Ly] = {
            'c_eff_fit': float(c_eff),
            'c_eff_pairwise': [float(c) for c in c_pairs],
            'S_values': [float(s) for s in S_arr],
            'Lx_values': [int(l) for l in Lx_arr],
        }

results['c_eff_results'] = {str(k): v for k, v in c_eff_results.items()}

# Save
with open("results/sprint_075b_cylinder_entropy.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_075b_cylinder_entropy.json", flush=True)

from db_utils import record
for Ly in sorted(c_eff_results.keys()):
    data = c_eff_results[Ly]
    if isinstance(data, dict) and data.get('c_eff_fit') is not None:
        record(sprint=75, model='hybrid_cyl', q=2, n=Ly,
               quantity='c_eff', value=data['c_eff_fit'],
               method=f'entropy_halfcut_Ly{Ly}',
               notes=f'Ly={Ly} cylinder at g_c, open-BC fit')

print(f"\n{'='*60}", flush=True)
print("SUMMARY", flush=True)
print(f"{'='*60}", flush=True)
for Ly in sorted(entropy_data.keys()):
    data = entropy_data[Ly]
    S_vals = [d['S'] for d in data]
    Lx_vals = [d['Lx'] for d in data]
    ceff = c_eff_results.get(Ly, {})
    c_str = f"c_eff={ceff.get('c_eff_fit', 'N/A'):.3f}" if isinstance(ceff, dict) and ceff.get('c_eff_fit') is not None else "c_eff=N/A"
    print(f"  Ly={Ly}: S={[f'{s:.3f}' for s in S_vals]} at Lx={Lx_vals}, {c_str}", flush=True)
