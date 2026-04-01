"""Sprint 058b: q=5 Potts spectrum at n=8 periodic (dim=390,625).

Extract E0 and first gap for absolute x1 extraction.
"""
import numpy as np
from scipy.sparse import csr_matrix, eye as sp_eye, kron as sp_kron
from scipy.sparse.linalg import eigsh
import json, time

def potts_hamiltonian_periodic(n, q, g):
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

print("=== Sprint 058b: q=5 Potts, n=8 periodic ===")
print(f"Hilbert space dim = 5^8 = {5**8}")

q = 5
g_c = 0.441  # from Sprint 052
n = 8

t0 = time.time()
print(f"Building Hamiltonian (dim={q**n})...")
H = potts_hamiltonian_periodic(n, q, g_c)
t_build = time.time() - t0
print(f"  Build time: {t_build:.1f}s, nnz={H.nnz}")

print("Running eigsh for 8 lowest eigenvalues...")
t1 = time.time()
evals, _ = eigsh(H, k=8, which='SA')
evals = np.sort(evals)
t_eig = time.time() - t1
print(f"  Eigsh time: {t_eig:.1f}s")

gaps = evals - evals[0]
print(f"\n  E0 = {evals[0]:.10f}")
print(f"  E0/n = {evals[0]/n:.10f}")
print(f"  Eigenvalues: {evals}")
print(f"  Gaps: {gaps}")

nonzero = gaps[gaps > 1e-8]
if len(nonzero) > 0:
    print(f"\n  First gap (Delta_1) = {nonzero[0]:.10f}")
    print(f"  Delta_1 * n = {nonzero[0] * n:.6f}")
    if len(nonzero) > 1:
        ratios = nonzero / nonzero[0]
        print(f"  Gap ratios: {[f'{r:.4f}' for r in ratios]}")

total_time = time.time() - t0
print(f"\nTotal time: {total_time:.1f}s")

results = {
    "q": q, "n": n, "g_c": g_c,
    "dim": q**n,
    "E0": float(evals[0]),
    "E0_per_site": float(evals[0] / n),
    "eigenvalues": [float(e) for e in evals],
    "gaps": [float(g) for g in gaps],
    "nonzero_gaps": [float(g) for g in nonzero],
    "build_time_s": round(t_build, 2),
    "eigsh_time_s": round(t_eig, 2),
    "total_time_s": round(total_time, 2)
}

with open("results/sprint_058b_q5_n8.json", "w") as f:
    json.dump(results, f, indent=2)
print("Results saved to results/sprint_058b_q5_n8.json")
