"""Sprint 057a: Validate CFT spectrum extraction on q=2 Ising (TFIM) with periodic BC.

At criticality on periodic chain of length L, CFT predicts:
  E_n - E_0 = (2*pi*v_s/L) * x_n
where x_n are scaling dimensions. Ratios R_n = x_n/x_1 are universal.

Ising CFT (c=1/2) primary scaling dimensions:
  Identity sector: 0, 2, 4, ...  (descendants)
  Spin field sigma: 1/8, 1+1/8, ...
  Energy field epsilon: 1, 2, 3, ...

On periodic chain, the low-lying spectrum should give ratios matching these.
We use the standard TFIM: H = -sum_i Z_i Z_{i+1} - g sum_i X_i, with g_c = 1.
"""
import numpy as np
from scipy.sparse import csr_matrix, eye as sp_eye, kron as sp_kron
from scipy.sparse.linalg import eigsh
import json, time

def tfim_hamiltonian_periodic(n, g):
    """Build TFIM Hamiltonian with periodic BC: H = -ZZ - g*X"""
    dim = 2**n
    sz = csr_matrix(np.array([[1,0],[0,-1]], dtype=float))
    sx = csr_matrix(np.array([[0,1],[1,0]], dtype=float))

    H = csr_matrix((dim, dim))

    # ZZ interaction (periodic)
    for i in range(n):
        j = (i + 1) % n
        if i < j:
            ops = [sp_eye(2, format='csr')] * n
            ops[i] = sz
            ops[j] = sz
        else:  # wrap-around: i=n-1, j=0
            ops = [sp_eye(2, format='csr')] * n
            ops[i] = sz
            ops[j] = sz

        result = ops[0]
        for k in range(1, n):
            result = sp_kron(result, ops[k], format='csr')
        H = H - result

    # Transverse field
    for i in range(n):
        ops = [sp_eye(2, format='csr')] * n
        ops[i] = sx
        result = ops[0]
        for k in range(1, n):
            result = sp_kron(result, ops[k], format='csr')
        H = H - g * result

    return H

# Test at g_c = 1 for several chain lengths
results = {}
for n in [8, 10, 12]:
    t0 = time.time()
    H = tfim_hamiltonian_periodic(n, g=1.0)

    # Get lowest 12 eigenvalues
    n_eig = min(12, 2**n - 2)
    evals, _ = eigsh(H, k=n_eig, which='SA')
    evals = np.sort(evals)

    dt = time.time() - t0

    # Energy gaps
    gaps = evals - evals[0]
    # Ratios relative to first gap
    ratios = gaps[1:] / gaps[1]  # R_1 = 1 by definition

    # Ground state energy per site (for v_s extraction)
    e0_per_site = evals[0] / n

    # Sound velocity: v_s = (E_1 - E_0) * L / (2*pi*x_1)
    # For Ising, x_1 = 1/8 (spin field in Ramond sector) or x_1 = 1 (energy field)
    # The lowest gap on periodic chain depends on the sector
    gap1 = gaps[1]

    print(f"\nn={n} (dim={2**n}), time={dt:.1f}s")
    print(f"  E0 = {evals[0]:.6f}, E0/n = {e0_per_site:.6f}")
    print(f"  First gap = {gap1:.6f}")
    print(f"  Gap ratios (E_n-E_0)/(E_1-E_0):")
    for i, r in enumerate(ratios):
        print(f"    R_{i+1} = {r:.4f}")

    results[f"n={n}"] = {
        "n": n, "dim": 2**n, "time_s": round(dt, 2),
        "E0": round(float(evals[0]), 8),
        "gaps": [round(float(g), 8) for g in gaps],
        "ratios": [round(float(r), 6) for r in ratios],
    }

# Ising CFT expected scaling dimensions (periodic chain, both sectors):
# The spectrum depends on boundary conditions and sector.
# For periodic chain (torus), the partition function gives towers from:
#   (h, h_bar) = (0,0), (1/2, 1/2), (1/16, 1/16), (1/2, 0), (0, 1/2)
# Scaling dimension x = h + h_bar:
#   x = 0 (identity), 1/8 (sigma), 1 (epsilon), 1/2 (fermion)
# On a spin chain, the spectrum includes states from all Z2 sectors.
# Expected first few gaps sorted: 1/8, 1, 1+1/8, 2, ...
# But actual gap ordering depends on chain parity and sector structure.

print("\n\nIsing CFT expected scaling dimensions (x values):")
print("  Identity: 0, 2, 2, 4, 4, ...")
print("  Sigma: 1/8 = 0.125, 1+1/8 = 1.125, ...")
print("  Epsilon: 1, 2, 3, ...")
print("  Fermion: 1/2 = 0.5, 3/2 = 1.5, ...")

# Save
with open("results/sprint_057a_ising_spectrum.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nResults saved.")
