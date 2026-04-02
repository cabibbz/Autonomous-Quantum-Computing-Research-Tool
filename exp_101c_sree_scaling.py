"""Sprint 101c: Size scaling of SREE — does equipartition improve with n?

For S_q Potts at g_c=1/q: measure p(0), equipartition CV, and S_number/S_total
as functions of system size n, at fixed nA/n=0.5 (midchain bipartition).

Key question: does p(0)*q → 1 as n → ∞ (equipartition restored)?
Or does the charge-0 enrichment persist?

Also test: nA=n//2 vs small fixed nA to disentangle subsystem vs system size effects.
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from gpu_utils import eigsh
import json, time

def potts_hamiltonian_periodic(n, q, g, model='S_q'):
    dim = q**n
    eye_q = np.eye(q)
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    if model == 'hybrid':
        field_op = X + X.T
    else:
        field_op = np.zeros((q, q))
        Xk = np.eye(q)
        for k in range(1, q):
            Xk = Xk @ X
            field_op += Xk
    delta_2 = np.zeros((q**2, q**2))
    for s in range(q):
        delta_2[s*q+s, s*q+s] = 1.0
    H = csr_matrix((dim, dim))
    for i in range(n):
        j = (i + 1) % n
        if j == i + 1:
            left = q**i if i > 0 else 1
            right = q**(n-j-1) if n-j-1 > 0 else 1
            op = csr_matrix(delta_2)
            if i > 0:
                op = sp_kron(sp_eye(left), op, format='csr')
            if n-j-1 > 0:
                op = sp_kron(op, sp_eye(right), format='csr')
            H = H - op
        else:
            for s in range(q):
                proj = np.zeros((q, q)); proj[s, s] = 1.0
                ops = [eye_q]*n; ops[0] = proj; ops[n-1] = proj
                r = csr_matrix(ops[0])
                for op in ops[1:]:
                    r = sp_kron(r, csr_matrix(op), format='csr')
                H = H - r
    for i in range(n):
        ops = [eye_q]*n; ops[i] = field_op
        r = csr_matrix(ops[0])
        for op in ops[1:]:
            r = sp_kron(r, csr_matrix(op), format='csr')
        H = H - g * r
    return H

def reduced_density_matrix(psi, n, q, nA):
    dimA = q**nA
    dimB = q**(n - nA)
    psi_mat = psi.reshape((dimA, dimB), order='F')
    return psi_mat @ psi_mat.conj().T

def compute_sree(rho_A, nA, q):
    dimA = q**nA
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    GA = X.copy()
    for _ in range(nA - 1):
        GA = np.kron(GA, X)

    omega = np.exp(2j * np.pi / q)
    GA_k = np.eye(dimA)
    GA_powers = [GA_k.copy()]
    for k in range(1, q):
        GA_k = GA_k @ GA
        GA_powers.append(GA_k.copy())

    # Total entropy
    evals_total = np.linalg.eigvalsh(rho_A)
    evals_total = evals_total[evals_total > 1e-15]
    S_total = -np.sum(evals_total * np.log(evals_total))

    p_alpha = []
    S_alpha = []
    for alpha in range(q):
        P = np.zeros((dimA, dimA), dtype=complex)
        for k in range(q):
            P += omega**(-alpha * k) * GA_powers[k]
        P /= q
        rho_alpha = P @ rho_A @ P.conj().T
        p = np.real(np.trace(rho_alpha))
        p_alpha.append(float(p))
        if p > 1e-15:
            rho_norm = rho_alpha / p
            evals_sec = np.linalg.eigvalsh(rho_norm)
            evals_sec = evals_sec[evals_sec > 1e-15]
            S_alpha.append(float(-np.sum(evals_sec * np.log(evals_sec))))
        else:
            S_alpha.append(0.0)

    p_arr = np.array(p_alpha)
    p_pos = p_arr[p_arr > 1e-15]
    S_number = float(-np.sum(p_pos * np.log(p_pos)))
    equi_cv = float(np.std(p_alpha) / np.mean(p_alpha)) if np.mean(p_alpha) > 0 else 0

    return {
        'S_total': float(S_total), 'S_number': float(S_number),
        'p_alpha': p_alpha, 'S_alpha': S_alpha,
        'p0': p_alpha[0], 'p0_times_q': p_alpha[0] * q,
        'equi_cv': equi_cv,
        'sn_ratio': float(S_number / S_total) if S_total > 0 else 0,
    }

# ============================================================
# Part 1: S_q Potts size scaling with midchain bipartition
# ============================================================
print("=" * 70)
print("Sprint 101c: SREE Size Scaling")
print("=" * 70)

all_results = {}

# q=2: n=4,6,8,10,12,14 (dim up to 16384)
# q=3: n=4,6,8,10 (dim up to 59049)
# q=5: n=4,6,8 (dim up to 390625), n=10 with GPU (9.8M)
# q=7: n=4,6 (dim up to 117649)

configs = [
    # (q, n_list) — midchain bipartition nA = n//2
    (2, [4, 6, 8, 10, 12, 14]),
    (3, [4, 6, 8, 10]),
    (5, [4, 6, 8]),
    (7, [4, 6]),
]

for q, n_list in configs:
    g_c = 1.0 / q
    print(f"\n--- q={q}, g_c={g_c:.4f} (S_q Potts, midchain bipartition) ---")
    for n in n_list:
        nA = n // 2
        dim = q**n
        dimA = q**nA
        if dim > 10_000_000:
            print(f"  n={n}: dim={dim} > 10M, skipping")
            continue
        if dimA > 5000:
            print(f"  n={n}: dimA={dimA} > 5000, projectors too expensive, skipping")
            continue

        t0 = time.time()
        H = potts_hamiltonian_periodic(n, q, g_c, model='S_q')
        E0, psi = eigsh(H, k=1, which='SA')
        E0 = E0[0]; psi = psi[:, 0]
        rho_A = reduced_density_matrix(psi, n, q, nA)
        sree = compute_sree(rho_A, nA, q)
        elapsed = time.time() - t0

        print(f"  n={n:>2} nA={nA}: p0*q={sree['p0_times_q']:.4f}, "
              f"CV={sree['equi_cv']:.5f}, S_n/S_t={sree['sn_ratio']:.4f} ({elapsed:.1f}s)")

        key = f"Sq_q{q}_n{n}_nA{nA}"
        all_results[key] = {'q': q, 'n': n, 'nA': nA, 'model': 'S_q', 'g_c': g_c,
                            'dim': dim, 'time_s': round(elapsed, 2), **sree}

# Part 2: Fixed small nA=2 with varying n — isolate system size effect
print("\n\n--- Fixed nA=2, varying n ---")
for q in [2, 3, 5]:
    g_c = 1.0 / q
    n_list = list(range(4, 14, 2)) if q == 2 else list(range(4, 10, 2))
    if q == 5:
        n_list = [4, 6, 8]
    for n in n_list:
        nA = 2
        dim = q**n
        if dim > 10_000_000:
            continue

        t0 = time.time()
        H = potts_hamiltonian_periodic(n, q, g_c, model='S_q')
        E0, psi = eigsh(H, k=1, which='SA')
        E0 = E0[0]; psi = psi[:, 0]
        rho_A = reduced_density_matrix(psi, n, q, nA)
        sree = compute_sree(rho_A, nA, q)
        elapsed = time.time() - t0

        print(f"  q={q} n={n:>2} nA={nA}: p0*q={sree['p0_times_q']:.4f}, "
              f"CV={sree['equi_cv']:.5f}, S_n/S_t={sree['sn_ratio']:.4f} ({elapsed:.1f}s)")

        key = f"Sq_q{q}_n{n}_nA{nA}"
        all_results[key] = {'q': q, 'n': n, 'nA': nA, 'model': 'S_q', 'g_c': g_c,
                            'dim': dim, 'time_s': round(elapsed, 2), **sree}

# Summary: convergence towards equipartition
print("\n\n" + "=" * 70)
print("SIZE SCALING SUMMARY: p(0)*q (=1 is equipartition)")
print("=" * 70)

for q in [2, 3, 5, 7]:
    print(f"\n  q={q} (midchain bipartition nA=n/2):")
    entries = [(k, v) for k, v in all_results.items()
               if v['q'] == q and v['nA'] == v['n'] // 2]
    entries.sort(key=lambda x: x[1]['n'])
    for key, r in entries:
        print(f"    n={r['n']:>2}, nA={r['nA']}: p0*q={r['p0_times_q']:.4f}, "
              f"CV={r['equi_cv']:.5f}")

for q in [2, 3, 5]:
    print(f"\n  q={q} (fixed nA=2):")
    entries = [(k, v) for k, v in all_results.items()
               if v['q'] == q and v['nA'] == 2]
    entries.sort(key=lambda x: x[1]['n'])
    for key, r in entries:
        print(f"    n={r['n']:>2}, nA=2: p0*q={r['p0_times_q']:.4f}, "
              f"CV={r['equi_cv']:.5f}")

with open("results/sprint_101c_sree_scaling.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print("\nResults saved to results/sprint_101c_sree_scaling.json")
