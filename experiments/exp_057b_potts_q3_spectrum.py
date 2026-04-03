"""Sprint 057b: CFT spectrum for q=3 Potts at g_c=1/3 with periodic BC.

q=3 Potts CFT has c=4/5 with W_3 algebra. Known primary scaling dimensions:
  (h, h_bar): (0,0), (2/5,2/5), (1/15,1/15), (2/3,2/3), (7/5,7/5), ...
  x = h + h_bar: 0, 4/5, 2/15, 4/3, 14/5, ...
  Sorted non-zero: 2/15 ≈ 0.133, 4/5 = 0.8, 4/3 ≈ 1.333, ...

Hamiltonian: H = -J sum_i delta(s_i, s_{i+1}) - g sum_i (X_i + X_i†)
g_c = 1/3 (exact, self-dual).
"""
import numpy as np
from scipy.sparse import csr_matrix, eye as sp_eye, kron as sp_kron
from scipy.sparse.linalg import eigsh
import json, time

def potts_hamiltonian_periodic(n, q, g):
    """Build q-state Potts Hamiltonian with periodic BC."""
    d = q
    dim = d**n

    # X operator: cyclic shift |s> -> |s+1 mod q>
    X = np.zeros((d, d))
    for s in range(d):
        X[(s+1) % d, s] = 1.0
    Xd = X.T  # X†
    XpXd = csr_matrix(X + Xd)

    # delta(s_i, s_j) as d^2 x d^2 matrix: sum_a |a,a><a,a|
    delta_op = np.zeros((d*d, d*d))
    for a in range(d):
        idx = a * d + a
        delta_op[idx, idx] = 1.0
    delta_sp = csr_matrix(delta_op)

    H = csr_matrix((dim, dim))

    # Interaction: -J delta(s_i, s_{i+1}), periodic BC
    for i in range(n):
        j = (i + 1) % n
        # Build full operator via kronecker products
        # Sites are ordered: 0, 1, ..., n-1
        if j == i + 1:  # adjacent (non-wrapping)
            left = sp_eye(d**i, format='csr') if i > 0 else csr_matrix(np.array([[1.0]]))
            right = sp_eye(d**(n-i-2), format='csr') if i + 2 < n else csr_matrix(np.array([[1.0]]))
            H = H - sp_kron(sp_kron(left, delta_sp, format='csr'), right, format='csr')
        else:  # wrapping: i=n-1, j=0
            # Need to handle non-adjacent kronecker product
            # delta(s_{n-1}, s_0) = sum_a |a>_{n-1}<a| ⊗ |a>_0<a|
            for a in range(d):
                # |a><a| on site 0
                proj0 = np.zeros((d, d))
                proj0[a, a] = 1.0
                # |a><a| on site n-1
                projn = np.zeros((d, d))
                projn[a, a] = 1.0
                # Full: proj0 ⊗ I^{n-2} ⊗ projn
                op = csr_matrix(proj0)
                middle = sp_eye(d**(n-2), format='csr')
                op = sp_kron(sp_kron(op, middle, format='csr'), csr_matrix(projn), format='csr')
                H = H - op

    # Transverse field: -g sum_i (X_i + X_i†)
    for i in range(n):
        left = sp_eye(d**i, format='csr') if i > 0 else csr_matrix(np.array([[1.0]]))
        right = sp_eye(d**(n-i-1), format='csr') if i < n-1 else csr_matrix(np.array([[1.0]]))
        H = H - g * sp_kron(sp_kron(left, XpXd, format='csr'), right, format='csr')

    return H

# q=3 Potts at g_c = 1/3
q = 3
g_c = 1.0 / 3.0
results = {}

for n in [6, 8, 10]:
    t0 = time.time()
    dim = q**n
    print(f"\nq={q}, n={n}, dim={dim}")

    H = potts_hamiltonian_periodic(n, q, g_c)

    n_eig = min(16, dim - 2)
    evals, _ = eigsh(H, k=n_eig, which='SA')
    evals = np.sort(evals)
    dt = time.time() - t0

    gaps = evals - evals[0]
    # Filter out near-zero gaps (degeneracies of ground state)
    nonzero = gaps[gaps > 1e-8]
    ratios = nonzero / nonzero[0]

    print(f"  Time: {dt:.1f}s")
    print(f"  E0 = {evals[0]:.6f}")
    print(f"  First {min(6, len(gaps))} gaps: {[f'{g:.6f}' for g in gaps[:6]]}")
    print(f"  Non-zero gaps: {[f'{g:.6f}' for g in nonzero[:10]]}")
    print(f"  Gap ratios:")
    for i, r in enumerate(ratios[:10]):
        print(f"    R_{i+1} = {r:.4f}")

    results[f"n={n}"] = {
        "n": n, "q": q, "g_c": g_c, "dim": dim, "time_s": round(dt, 2),
        "E0": round(float(evals[0]), 8),
        "eigenvalues": [round(float(e), 8) for e in evals],
        "gaps": [round(float(g), 8) for g in gaps],
        "nonzero_gaps": [round(float(g), 8) for g in nonzero[:12]],
        "ratios": [round(float(r), 6) for r in ratios[:12]],
    }

# Expected q=3 Potts CFT scaling dimensions
# Primaries: x = 2/15 ≈ 0.133, 4/5 = 0.8, 4/3 ≈ 1.333
# Ratios relative to x_1 = 2/15: x_2/x_1 = (4/5)/(2/15) = 6, x_3/x_1 = (4/3)/(2/15) = 10
print("\n\nq=3 Potts CFT (c=4/5) expected scaling dimensions:")
print(f"  x_1 = 2/15 = {2/15:.4f} (spin field sigma)")
print(f"  x_2 = 4/5 = {4/5:.4f} (energy field epsilon)")
print(f"  x_3 = 4/3 = {4/3:.4f} (subleading)")
print(f"  Expected ratios: R_2 = {(4/5)/(2/15):.2f}, R_3 = {(4/3)/(2/15):.2f}")

with open("results/sprint_057b_q3_spectrum.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nResults saved.")
