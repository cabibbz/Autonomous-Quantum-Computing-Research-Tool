#!/usr/bin/env python3
"""Sprint 069a: 2D entanglement entropy profiles at q=2 (validation)

Compute ground state on L×L periodic torus at g_c(2D)=0.771,
then measure entanglement entropy for strip partitions of width w=1,...,L-1.

Strip partition: subsystem A = all sites with x-coordinate in [0, w).
Boundary length = 2*L (two cuts across the periodic y-direction).

Lattices: L=2 (4 sites, dim=16), L=3 (9 sites, dim=512), L=4 (16 sites, dim=65536)
All CPU-feasible for q=2.

Also measure at g=0.1 (ordered) and g=3.0 (disordered) for comparison.
"""
import numpy as np
import json
import time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigvalsh

def build_H_2d(Lx, Ly, q, g):
    """Build 2D hybrid Hamiltonian on Lx x Ly periodic square lattice."""
    n = Lx * Ly
    dim = q**n

    potts_2site = np.zeros(q**2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q**2, q**2), format='csr')

    X = np.zeros((q, q))
    for s in range(q):
        X[(s + 1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)

    H = csr_matrix((dim, dim))

    # Enumerate bonds (periodic square lattice)
    bonds = []
    for y in range(Ly):
        for x in range(Lx):
            site = y * Lx + x
            nbr_h = y * Lx + (x + 1) % Lx
            if nbr_h != site:
                bonds.append((min(site, nbr_h), max(site, nbr_h)))
            nbr_v = ((y + 1) % Ly) * Lx + x
            if nbr_v != site:
                bonds.append((min(site, nbr_v), max(site, nbr_v)))
    bonds = list(set(bonds))

    for (i, j) in bonds:
        if j == i + 1:
            left = q**i
            right = q**(n - i - 2)
            op = potts_op
            if left > 1:
                op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1:
                op = sp_kron(op, sp_eye(right), format='csr')
        else:
            diag_vals = np.zeros(dim)
            for idx in range(dim):
                si = (idx // (q**i)) % q
                sj = (idx // (q**j)) % q
                diag_vals[idx] = 1.0 if si == sj else 0.0
            op = diags(diag_vals, 0, shape=(dim, dim), format='csr')
        H = H - op

    for i in range(n):
        left = q**i
        right = q**(n - i - 1)
        op = XpXd.copy()
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op

    return H


def get_ground_state(H, dim):
    """Get ground state vector."""
    if dim <= 2000:
        Hd = H.toarray()
        evals, evecs = np.linalg.eigh(Hd)
        return evals[0], evecs[:, 0]
    else:
        evals, evecs = eigsh(H, k=1, which='SA')
        return evals[0], evecs[:, 0]


def entanglement_entropy(psi, q, n, subsystem_A):
    """Compute von Neumann entropy of subsystem A.

    psi: ground state vector of length q^n
    subsystem_A: list of site indices in subsystem A
    """
    nA = len(subsystem_A)
    nB = n - nA
    subsystem_B = [i for i in range(n) if i not in subsystem_A]

    dimA = q**nA
    dimB = q**nB

    # Reshape psi into a matrix (A x B) by reordering indices
    # Site ordering: index = sum_i s_i * q^i
    # We need to group A sites together and B sites together

    # Build the reshaped state
    psi_matrix = np.zeros((dimA, dimB))

    for idx in range(q**n):
        # Extract site values
        site_vals = []
        tmp = idx
        for i in range(n):
            site_vals.append(tmp % q)
            tmp //= q

        # Compute A index and B index
        a_idx = 0
        for ki, site in enumerate(subsystem_A):
            a_idx += site_vals[site] * (q**ki)

        b_idx = 0
        for ki, site in enumerate(subsystem_B):
            b_idx += site_vals[site] * (q**ki)

        psi_matrix[a_idx, b_idx] += psi[idx]

    # SVD to get entanglement spectrum
    s = np.linalg.svd(psi_matrix, compute_uv=False)
    s = s[s > 1e-15]
    probs = s**2
    probs = probs / probs.sum()  # normalize

    # Von Neumann entropy
    S = -np.sum(probs * np.log(probs))
    return S


def strip_partition(Lx, Ly, w):
    """Return list of sites in a strip of width w (x in [0, w))."""
    sites = []
    for y in range(Ly):
        for x in range(w):
            sites.append(y * Lx + x)
    return sites


# ---- MAIN ----
results = {
    'experiment': '069a_2d_entropy_q2',
    'description': '2D entanglement entropy profiles at q=2 critical point',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': 2,
}

q = 2
gc_2d = 0.771  # From Sprint 068

print("=" * 60)
print("Sprint 069a: 2D Entanglement Entropy at q=2")
print("=" * 60)

# g values: critical, ordered, disordered
g_values = {'critical': gc_2d, 'ordered': 0.1, 'disordered': 3.0}

for L in [2, 3, 4]:
    n = L * L
    dim = q**n
    print(f"\n{'='*40}")
    print(f"L={L} ({n} sites, dim={dim})")
    print(f"{'='*40}")

    L_data = {}

    for phase, g in g_values.items():
        print(f"\n  g={g:.3f} ({phase}):")
        t0 = time.time()
        H = build_H_2d(L, L, q, g)
        E0, psi = get_ground_state(H, dim)
        dt_H = time.time() - t0
        print(f"    Ground state: E0/n = {E0/n:.6f} ({dt_H:.1f}s)")

        # Entropy profile: strip of width w=1,...,L-1
        entropy_profile = []
        for w in range(1, L):
            A_sites = strip_partition(L, L, w)
            nA = len(A_sites)
            boundary_length = 2 * L  # Two cuts across periodic y-direction

            t1 = time.time()
            S = entanglement_entropy(psi, q, n, A_sites)
            dt_S = time.time() - t1

            print(f"    w={w}: |A|={nA}, boundary={boundary_length}, S={S:.6f} ({dt_S:.2f}s)")
            entropy_profile.append({
                'w': w,
                'nA': nA,
                'boundary_length': boundary_length,
                'S': round(float(S), 8),
                'S_per_boundary': round(float(S / boundary_length), 8),
            })

        L_data[phase] = {
            'g': g,
            'E0_per_site': round(float(E0/n), 8),
            'entropy_profile': entropy_profile,
        }

    results[f'L{L}'] = L_data

# ---- Analysis ----
print("\n" + "=" * 60)
print("ANALYSIS")
print("=" * 60)

# 1. Area law check: S vs boundary length at criticality
print("\n1. Entropy at critical point vs strip width:")
print(f"   {'L':>3} {'w':>3} {'|A|':>4} {'bdry':>5} {'S':>10} {'S/bdry':>10}")
for L in [2, 3, 4]:
    for entry in results[f'L{L}']['critical']['entropy_profile']:
        print(f"   {L:3d} {entry['w']:3d} {entry['nA']:4d} {entry['boundary_length']:5d} "
              f"{entry['S']:10.6f} {entry['S_per_boundary']:10.6f}")

# 2. Critical vs gapped: entropy should be much larger at g_c
print("\n2. Critical vs gapped entropy (half-system or w=L//2):")
for L in [2, 3, 4]:
    w_half = max(1, L // 2)
    idx = w_half - 1  # w starts at 1
    if idx < len(results[f'L{L}']['critical']['entropy_profile']):
        S_crit = results[f'L{L}']['critical']['entropy_profile'][idx]['S']
        S_ord = results[f'L{L}']['ordered']['entropy_profile'][idx]['S']
        S_dis = results[f'L{L}']['disordered']['entropy_profile'][idx]['S']
        print(f"   L={L}, w={w_half}: S_ordered={S_ord:.6f}, S_critical={S_crit:.6f}, S_disordered={S_dis:.6f}")

# 3. S/boundary at criticality across sizes (should be roughly constant = area law)
print("\n3. S/boundary_length at criticality (area law test):")
print("   Fixed w=1 strip (boundary=2L):")
for L in [2, 3, 4]:
    entry = results[f'L{L}']['critical']['entropy_profile'][0]  # w=1
    print(f"   L={L}: S={entry['S']:.6f}, boundary={entry['boundary_length']}, S/bdry={entry['S_per_boundary']:.6f}")

# Save
with open("results/sprint_069a_2d_entropy_q2.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_069a_2d_entropy_q2.json")
