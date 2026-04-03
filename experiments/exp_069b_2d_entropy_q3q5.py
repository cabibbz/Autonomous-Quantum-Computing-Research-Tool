#!/usr/bin/env python3
"""Sprint 069b: 2D entanglement entropy profiles at q=3,5 critical points

q=3: L=2 (dim=81), L=3 (dim=19683) — CPU
q=5: L=2 (dim=625), L=3 (dim=1953125) — GPU for L=3

Measure strip entropy S(w) at 2D g_c, plus ordered and disordered phases.
"""
import numpy as np
import json
import time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh

try:
    import cupy as cp
    from cupyx.scipy.sparse import csr_matrix as cp_csr
    from cupyx.scipy.sparse.linalg import eigsh as cp_eigsh
    HAS_GPU = True
    print("GPU available")
except ImportError:
    HAS_GPU = False
    print("GPU not available")

GPU_THRESHOLD = 50000


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
    """Get ground state vector, using GPU for large systems."""
    if dim <= 2000:
        Hd = H.toarray()
        evals, evecs = np.linalg.eigh(Hd)
        return evals[0], evecs[:, 0]

    if HAS_GPU and dim > GPU_THRESHOLD:
        try:
            H_gpu = cp_csr(H)
            evals_gpu, evecs_gpu = cp_eigsh(H_gpu, k=1, which='SA')
            E0 = float(cp.asnumpy(evals_gpu)[0])
            psi = cp.asnumpy(evecs_gpu[:, 0])
            return E0, psi
        except Exception as e:
            print(f"    GPU failed ({e}), falling back to CPU")

    evals, evecs = eigsh(H, k=1, which='SA')
    return evals[0], evecs[:, 0]


def entanglement_entropy(psi, q, n, subsystem_A):
    """Compute von Neumann entropy of subsystem A."""
    nA = len(subsystem_A)
    subsystem_B = [i for i in range(n) if i not in subsystem_A]
    dimA = q**nA
    dimB = q**(n - nA)

    psi_matrix = np.zeros((dimA, dimB))
    for idx in range(q**n):
        site_vals = []
        tmp = idx
        for i in range(n):
            site_vals.append(tmp % q)
            tmp //= q

        a_idx = sum(site_vals[site] * (q**ki) for ki, site in enumerate(subsystem_A))
        b_idx = sum(site_vals[site] * (q**ki) for ki, site in enumerate(subsystem_B))
        psi_matrix[a_idx, b_idx] += psi[idx]

    s = np.linalg.svd(psi_matrix, compute_uv=False)
    s = s[s > 1e-15]
    probs = s**2
    probs = probs / probs.sum()
    S = -np.sum(probs * np.log(probs))
    return S


def strip_partition(Lx, Ly, w):
    """Sites in a strip of width w (x in [0, w))."""
    sites = []
    for y in range(Ly):
        for x in range(w):
            sites.append(y * Lx + x)
    return sites


# ---- MAIN ----
results = {
    'experiment': '069b_2d_entropy_q3q5',
    'description': '2D entanglement entropy profiles at q=3,5 critical points',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

configs = [
    {'q': 3, 'gc_2d': 1.267, 'sizes': [2, 3], 'g_ord': 0.1, 'g_dis': 5.0},
    {'q': 5, 'gc_2d': 1.588, 'sizes': [2, 3], 'g_ord': 0.1, 'g_dis': 5.0},
]

for cfg in configs:
    q = cfg['q']
    gc = cfg['gc_2d']
    print(f"\n{'='*60}")
    print(f"q={q}, g_c(2D)={gc}")
    print(f"{'='*60}")

    q_data = {}

    for L in cfg['sizes']:
        n = L * L
        dim = q**n
        print(f"\n  L={L} ({n} sites, dim={dim})")

        g_values = {'critical': gc, 'ordered': cfg['g_ord'], 'disordered': cfg['g_dis']}
        L_data = {}

        for phase, g in g_values.items():
            print(f"\n    g={g:.3f} ({phase}):")
            t0 = time.time()
            H = build_H_2d(L, L, q, g)
            dt_build = time.time() - t0
            print(f"      H built ({dt_build:.1f}s)")

            t0 = time.time()
            E0, psi = get_ground_state(H, dim)
            dt_gs = time.time() - t0
            print(f"      GS: E0/n={E0/n:.6f} ({dt_gs:.1f}s)")

            entropy_profile = []
            for w in range(1, L):
                A_sites = strip_partition(L, L, w)
                nA = len(A_sites)
                boundary = 2 * L

                t1 = time.time()
                S = entanglement_entropy(psi, q, n, A_sites)
                dt_S = time.time() - t1

                print(f"      w={w}: |A|={nA}, bdry={boundary}, S={S:.6f} ({dt_S:.1f}s)")
                entropy_profile.append({
                    'w': w, 'nA': nA, 'boundary_length': boundary,
                    'S': round(float(S), 8),
                    'S_per_boundary': round(float(S / boundary), 8),
                })

            L_data[phase] = {
                'g': g,
                'E0_per_site': round(float(E0/n), 8),
                'entropy_profile': entropy_profile,
            }

        q_data[f'L{L}'] = L_data

    results[f'q{q}'] = q_data

# ---- Analysis ----
print("\n" + "=" * 60)
print("ANALYSIS")
print("=" * 60)

print("\n1. Critical entropy by q and L:")
print(f"   {'q':>3} {'L':>3} {'w':>3} {'bdry':>5} {'S':>10} {'S/bdry':>10}")
for q_val in [3, 5]:
    key = f'q{q_val}'
    for L in [2, 3]:
        Lkey = f'L{L}'
        if Lkey in results[key]:
            for entry in results[key][Lkey]['critical']['entropy_profile']:
                print(f"   {q_val:3d} {L:3d} {entry['w']:3d} {entry['boundary_length']:5d} "
                      f"{entry['S']:10.6f} {entry['S_per_boundary']:10.6f}")

# Compare with q=2 (from 069a)
print("\n2. S/boundary at w=1 across q values (area law coefficient):")
q2_data = {2: 0.087, 3: 0.392, 4: 0.477}  # L=2,3,4 from 069a
print(f"   q=2: L=3 S/bdry=0.065, L=4 S/bdry=0.060")
for q_val in [3, 5]:
    key = f'q{q_val}'
    for L in [2, 3]:
        Lkey = f'L{L}'
        if Lkey in results[key]:
            prof = results[key][Lkey]['critical']['entropy_profile']
            if prof:
                print(f"   q={q_val}: L={L} S/bdry={prof[0]['S_per_boundary']:.6f}")

print("\n3. Ordered phase entropy (GHZ-like, should be ~ln(q)):")
for q_val in [3, 5]:
    key = f'q{q_val}'
    for L in [2, 3]:
        Lkey = f'L{L}'
        if Lkey in results[key]:
            prof = results[key][Lkey]['ordered']['entropy_profile']
            if prof:
                S_ord = prof[0]['S']
                print(f"   q={q_val}, L={L}: S_ord={S_ord:.6f}, ln(q)={np.log(q_val):.6f}")

# Save
with open("results/sprint_069b_2d_entropy_q3q5.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_069b_2d_entropy_q3q5.json")
