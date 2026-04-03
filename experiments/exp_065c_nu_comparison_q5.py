#!/usr/bin/env python3
"""Sprint 065c: ν comparison between hybrid and clock at q=5.

Method: energy gap slope d(Δ·N)/dg ~ N^{1/ν} at g_c.
Use n=4,6,8 for both models. Same method as Sprint 053.
If ν matches → strong evidence for same universality class.
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh
import json, time

q = 5

def hybrid_hamiltonian_periodic(n, q, g):
    dim = q**n
    potts_2site = np.zeros(q**2)
    for a in range(q):
        potts_2site[a*q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q**2, q**2), format='csr')
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)
    H = csr_matrix((dim, dim))
    for i in range(n - 1):
        left = q**i; right = q**(n - i - 2)
        op = potts_op
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op
    diag_boundary = np.zeros(dim)
    for idx in range(dim):
        s0 = idx % q; sn = (idx // (q**(n-1))) % q
        diag_boundary[idx] = 1.0 if s0 == sn else 0.0
    H = H - diags(diag_boundary, 0, shape=(dim, dim), format='csr')
    for i in range(n):
        left = q**i; right = q**(n - i - 1)
        op = XpXd
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H

def clock_hamiltonian_periodic(n, q, g):
    dim = q**n
    clock_2site = np.zeros(q**2)
    for a in range(q):
        for b in range(q):
            clock_2site[a*q + b] = np.cos(2*np.pi*(a-b)/q)
    clock_op = diags(clock_2site, 0, shape=(q**2, q**2), format='csr')
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)
    H = csr_matrix((dim, dim))
    for i in range(n - 1):
        left = q**i; right = q**(n - i - 2)
        op = clock_op
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op
    diag_boundary = np.zeros(dim)
    for idx in range(dim):
        s0 = idx % q; sn = (idx // (q**(n-1))) % q
        diag_boundary[idx] = np.cos(2*np.pi*(sn - s0)/q)
    H = H - diags(diag_boundary, 0, shape=(dim, dim), format='csr')
    for i in range(n):
        left = q**i; right = q**(n - i - 1)
        op = XpXd
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H

def gap_slope(ham_fn, q, gc, n, dg=0.005):
    """Compute d(Δ·N)/dg via central difference at g_c."""
    gaps = {}
    for g in [gc - dg, gc, gc + dg]:
        H = ham_fn(n, q, g)
        k_eig = min(4, q**n - 2)
        evals, _ = eigsh(H, k=k_eig, which='SA')
        evals = np.sort(evals)
        gaps[g] = (evals[1] - evals[0]) * n
    slope = (gaps[gc + dg] - gaps[gc - dg]) / (2 * dg)
    return slope, gaps[gc]

# ============================================================
# ν EXTRACTION
# ============================================================

results = {"q": q}

for model_name, ham_fn, gc in [("HYBRID", hybrid_hamiltonian_periodic, 0.441),
                                 ("CLOCK", clock_hamiltonian_periodic, 0.52)]:
    print(f"\n{'='*60}")
    print(f"{model_name} q={q}: ν EXTRACTION at g_c={gc}")
    print(f"{'='*60}")

    slopes = {}
    gaps_at_gc = {}
    sizes = [4, 6, 8]

    for n in sizes:
        dim = q**n
        if dim > 500_000:
            print(f"  n={n}: dim={dim} too large, skip")
            continue
        t0 = time.time()
        slope, gap = gap_slope(ham_fn, q, gc, n)
        dt = time.time() - t0
        slopes[n] = abs(slope)
        gaps_at_gc[n] = gap
        print(f"  n={n}: |d(Δ·N)/dg| = {abs(slope):.4f}, Δ·N = {gap:.6f} ({dt:.1f}s)")

    # Pairwise ν from slope ratios: slope ~ N^{1/ν}
    # slope(n2)/slope(n1) = (n2/n1)^{1/ν}
    # 1/ν = log(slope2/slope1) / log(n2/n1)
    ns = sorted(slopes.keys())
    nu_pairs = {}
    print(f"\n  Pairwise ν:")
    for i in range(len(ns)):
        for j in range(i+1, len(ns)):
            n1, n2 = ns[i], ns[j]
            if slopes[n1] > 0 and slopes[n2] > 0:
                inv_nu = np.log(slopes[n2] / slopes[n1]) / np.log(n2 / n1)
                if inv_nu > 0:
                    nu = 1.0 / inv_nu
                    nu_pairs[(n1, n2)] = nu
                    print(f"    ({n1},{n2}): ν = {nu:.4f}")
                else:
                    print(f"    ({n1},{n2}): inv_nu = {inv_nu:.4f} (negative, slope non-monotonic)")

    # Corrected ν (Sprint 053 method: b=0.86 correction)
    # slope(N) = A * N^{1/ν} * (1 + b/N), b = 0.86
    b_corr = 0.86
    print(f"\n  Corrected ν (b={b_corr}):")
    for i in range(len(ns)):
        for j in range(i+1, len(ns)):
            n1, n2 = ns[i], ns[j]
            if slopes[n1] > 0 and slopes[n2] > 0:
                # corrected_slope(N) = slope(N) / (1 + b/N)
                s1_corr = slopes[n1] / (1 + b_corr/n1)
                s2_corr = slopes[n2] / (1 + b_corr/n2)
                inv_nu = np.log(s2_corr / s1_corr) / np.log(n2 / n1)
                if inv_nu > 0:
                    nu = 1.0 / inv_nu
                    print(f"    ({n1},{n2}): ν_corr = {nu:.4f}")
                    nu_pairs[f"({n1},{n2})_corr"] = nu

    results[model_name.lower()] = {
        "gc": gc,
        "slopes": {str(n): float(v) for n, v in slopes.items()},
        "gaps_at_gc": {str(n): float(v) for n, v in gaps_at_gc.items()},
        "nu_pairs": {str(k): float(v) for k, v in nu_pairs.items()},
    }

# ============================================================
# SUMMARY
# ============================================================

print(f"\n{'='*60}")
print("ν COMPARISON SUMMARY: q=5")
print(f"{'='*60}")

h_data = results.get("hybrid", {})
c_data = results.get("clock", {})

h_nu = h_data.get("nu_pairs", {})
c_nu = c_data.get("nu_pairs", {})

print(f"\n{'Pair':<20} {'Hybrid ν':>12} {'Clock ν':>12}")
print("-" * 46)
for key in sorted(set(list(h_nu.keys()) + list(c_nu.keys()))):
    h_v = h_nu.get(key, None)
    c_v = c_nu.get(key, None)
    h_str = f"{h_v:.4f}" if h_v else "—"
    c_str = f"{c_v:.4f}" if c_v else "—"
    print(f"{key:<20} {h_str:>12} {c_str:>12}")

with open("results/sprint_065c_nu_comparison_q5.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_065c_nu_comparison_q5.json")
