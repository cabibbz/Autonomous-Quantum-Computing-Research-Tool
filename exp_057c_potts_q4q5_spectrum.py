"""Sprint 057c: CFT spectrum for q=4 and q=5 Potts with periodic BC.

q=4 (c=1, Ashkin-Teller/marginal): g_c = 0.392
q=5 (c≈1.10, NOVEL): g_c = 0.441

For q=4, max feasible: n=6 (dim = 4^6 = 4096)
For q=5, max feasible: n=6 (dim = 5^6 = 15625)
"""
import numpy as np
from scipy.sparse import csr_matrix, eye as sp_eye, kron as sp_kron
from scipy.sparse.linalg import eigsh
import json, time

def potts_hamiltonian_periodic(n, q, g):
    """Build q-state Potts Hamiltonian with periodic BC."""
    d = q
    dim = d**n

    X = np.zeros((d, d))
    for s in range(d):
        X[(s+1) % d, s] = 1.0
    XpXd = csr_matrix(X + X.T)

    delta_op = np.zeros((d*d, d*d))
    for a in range(d):
        idx = a * d + a
        delta_op[idx, idx] = 1.0
    delta_sp = csr_matrix(delta_op)

    H = csr_matrix((dim, dim))

    for i in range(n):
        j = (i + 1) % n
        if j == i + 1:
            left = sp_eye(d**i, format='csr') if i > 0 else csr_matrix(np.array([[1.0]]))
            right = sp_eye(d**(n-i-2), format='csr') if i + 2 < n else csr_matrix(np.array([[1.0]]))
            H = H - sp_kron(sp_kron(left, delta_sp, format='csr'), right, format='csr')
        else:
            for a in range(d):
                proj0 = np.zeros((d, d)); proj0[a, a] = 1.0
                projn = np.zeros((d, d)); projn[a, a] = 1.0
                op = sp_kron(sp_kron(csr_matrix(proj0), sp_eye(d**(n-2), format='csr'), format='csr'),
                             csr_matrix(projn), format='csr')
                H = H - op

    for i in range(n):
        left = sp_eye(d**i, format='csr') if i > 0 else csr_matrix(np.array([[1.0]]))
        right = sp_eye(d**(n-i-1), format='csr') if i < n-1 else csr_matrix(np.array([[1.0]]))
        H = H - g * sp_kron(sp_kron(left, XpXd, format='csr'), right, format='csr')

    return H

results = {}

for q, g_c in [(4, 0.392), (5, 0.441)]:
    for n in [4, 6]:
        dim = q**n
        if dim > 20000:
            if n > 6:
                print(f"Skipping q={q}, n={n} (dim={dim} too large)")
                continue

        t0 = time.time()
        print(f"\nq={q}, n={n}, dim={dim}, g_c={g_c}")

        H = potts_hamiltonian_periodic(n, q, g_c)
        n_eig = min(20, dim - 2)
        evals, _ = eigsh(H, k=n_eig, which='SA')
        evals = np.sort(evals)
        dt = time.time() - t0

        gaps = evals - evals[0]
        nonzero = gaps[gaps > 1e-8]
        ratios = nonzero / nonzero[0] if len(nonzero) > 0 else []

        print(f"  Time: {dt:.1f}s")
        print(f"  E0 = {evals[0]:.6f}")
        print(f"  Ground state degeneracy: {np.sum(gaps < 1e-6)}")
        print(f"  First 8 non-zero gaps: {[f'{g:.6f}' for g in nonzero[:8]]}")
        print(f"  Gap ratios:")
        for i, r in enumerate(ratios[:12]):
            print(f"    R_{i+1} = {r:.4f}")

        key = f"q={q}_n={n}"
        results[key] = {
            "q": q, "n": n, "g_c": g_c, "dim": dim, "time_s": round(dt, 2),
            "E0": round(float(evals[0]), 8),
            "eigenvalues": [round(float(e), 8) for e in evals],
            "gaps": [round(float(g), 8) for g in gaps],
            "nonzero_gaps": [round(float(g), 8) for g in nonzero[:15]],
            "ratios": [round(float(r), 6) for r in ratios[:15]],
            "gs_degeneracy": int(np.sum(gaps < 1e-6)),
        }

with open("results/sprint_057c_q4q5_spectrum.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nResults saved.")
