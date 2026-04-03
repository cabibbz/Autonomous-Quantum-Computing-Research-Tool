"""Sprint 057d: CFT spectrum for q=7 and q=4(n=8) Potts with periodic BC.

q=7 at n=4 (dim=2401) and n=6 (dim=117649 — might be too large, test timing)
q=4 at n=8 (dim=65536) for better convergence of ratios

Key question: does q=7 show THREE pairs of spin fields (σ,σ†), (σ²,σ⁵), (σ³,σ⁴)?
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

# q=4, n=8 for better convergence
print("=== q=4, n=8 ===")
t0 = time.time()
H = potts_hamiltonian_periodic(8, 4, 0.392)
print(f"  Build time: {time.time()-t0:.1f}s, dim={4**8}")
t1 = time.time()
evals, _ = eigsh(H, k=20, which='SA')
evals = np.sort(evals)
dt = time.time() - t0
gaps = evals - evals[0]
nonzero = gaps[gaps > 1e-8]
ratios = nonzero / nonzero[0]
print(f"  Diag time: {time.time()-t1:.1f}s")
print(f"  Gap ratios:")
for i, r in enumerate(ratios[:12]):
    print(f"    R_{i+1} = {r:.4f}")
results["q=4_n=8"] = {
    "q": 4, "n": 8, "g_c": 0.392, "dim": 4**8, "time_s": round(dt, 2),
    "eigenvalues": [round(float(e), 8) for e in evals],
    "nonzero_gaps": [round(float(g), 8) for g in nonzero[:15]],
    "ratios": [round(float(r), 6) for r in ratios[:15]],
}

# q=7, n=4
print("\n=== q=7, n=4 ===")
t0 = time.time()
H = potts_hamiltonian_periodic(4, 7, 0.535)
evals, _ = eigsh(H, k=20, which='SA')
evals = np.sort(evals)
dt = time.time() - t0
gaps = evals - evals[0]
nonzero = gaps[gaps > 1e-8]
ratios = nonzero / nonzero[0]
print(f"  Time: {dt:.1f}s, dim={7**4}")
print(f"  Gap ratios:")
for i, r in enumerate(ratios[:15]):
    print(f"    R_{i+1} = {r:.4f}")
results["q=7_n=4"] = {
    "q": 7, "n": 4, "g_c": 0.535, "dim": 7**4, "time_s": round(dt, 2),
    "eigenvalues": [round(float(e), 8) for e in evals],
    "nonzero_gaps": [round(float(g), 8) for g in nonzero[:18]],
    "ratios": [round(float(r), 6) for r in ratios[:18]],
}

# q=7, n=6 — test if feasible
print("\n=== q=7, n=6 (testing feasibility) ===")
dim76 = 7**6
print(f"  dim = {dim76}")
if dim76 <= 150000:
    t0 = time.time()
    H = potts_hamiltonian_periodic(6, 7, 0.535)
    print(f"  Build time: {time.time()-t0:.1f}s")
    t1 = time.time()
    evals, _ = eigsh(H, k=20, which='SA')
    evals = np.sort(evals)
    dt = time.time() - t0
    gaps = evals - evals[0]
    nonzero = gaps[gaps > 1e-8]
    ratios = nonzero / nonzero[0]
    print(f"  Diag time: {time.time()-t1:.1f}s")
    print(f"  Gap ratios:")
    for i, r in enumerate(ratios[:15]):
        print(f"    R_{i+1} = {r:.4f}")
    results["q=7_n=6"] = {
        "q": 7, "n": 6, "g_c": 0.535, "dim": dim76, "time_s": round(dt, 2),
        "eigenvalues": [round(float(e), 8) for e in evals],
        "nonzero_gaps": [round(float(g), 8) for g in nonzero[:18]],
        "ratios": [round(float(r), 6) for r in ratios[:18]],
    }
else:
    print("  Too large, skipping")

with open("results/sprint_057d_q7_spectrum.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nResults saved.")
