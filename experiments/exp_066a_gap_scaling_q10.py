#!/usr/bin/env python3
"""Sprint 066a: Δ·N scaling at q=10 for n=4,5,6,7.

Key diagnostic: if the transition is truly continuous, Δ·N converges to a constant.
If weakly first-order, Δ·N eventually diverges (gap stays finite at large N).

Uses GPU (CuPy) for n=6 (dim=10^6) and n=7 (dim=10^7).
Also compute Δ₂/Δ₁ ratio — universal at continuous transitions, → 1 at first-order.
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh
import json, time

q = 10
gc = 0.684  # Sprint 052 formula + Sprint 051 gap crossing

def hybrid_hamiltonian_periodic(n, q, g, use_gpu=False):
    """Build H = -Σ δ(s_i,s_j) - g Σ (X+X†) with periodic BC."""
    dim = q**n

    # Potts coupling: δ(s_i, s_j)
    potts_2site = np.zeros(q**2)
    for a in range(q):
        potts_2site[a*q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q**2, q**2), format='csr')

    # Clock field: X + X†
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)

    H = csr_matrix((dim, dim))

    # Nearest-neighbor Potts coupling (bulk)
    for i in range(n - 1):
        left = q**i; right = q**(n - i - 2)
        op = potts_op
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op

    # Periodic boundary: sites 0 and n-1
    diag_boundary = np.zeros(dim)
    for idx in range(dim):
        s0 = idx % q
        sn = (idx // (q**(n-1))) % q
        diag_boundary[idx] = 1.0 if s0 == sn else 0.0
    H = H - diags(diag_boundary, 0, shape=(dim, dim), format='csr')

    # Transverse field
    for i in range(n):
        left = q**i; right = q**(n - i - 1)
        op = XpXd
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op

    if use_gpu:
        import cupy as cp
        from cupyx.scipy.sparse import csr_matrix as cp_csr
        return cp_csr(H)
    return H

def get_lowest_eigenvalues(H, k=6, use_gpu=False):
    """Get k lowest eigenvalues, using GPU if requested."""
    if use_gpu:
        import cupy as cp
        from cupyx.scipy.sparse.linalg import eigsh as cp_eigsh
        evals, _ = cp_eigsh(H, k=k, which='SA')
        return np.sort(cp.asnumpy(evals))
    else:
        evals, _ = eigsh(H, k=k, which='SA')
        return np.sort(evals)

# ============================================================
# MAIN: Δ·N at g_c for multiple sizes
# ============================================================

print(f"Sprint 066a: Δ·N scaling at q={q}, g_c={gc}")
print(f"{'='*60}")

sizes = [4, 5, 6, 7]
results = {"q": q, "gc": gc, "sizes": {}}

for n in sizes:
    dim = q**n
    use_gpu = dim > 500_000
    k_eig = min(6, dim - 2)

    print(f"\nn={n}: dim={dim:,}, {'GPU' if use_gpu else 'CPU'}, k={k_eig}")

    t0 = time.time()
    try:
        H = hybrid_hamiltonian_periodic(n, q, gc, use_gpu=use_gpu)
        evals = get_lowest_eigenvalues(H, k=k_eig, use_gpu=use_gpu)
        dt = time.time() - t0

        E0 = evals[0]
        E1 = evals[1]
        gap1 = E1 - E0
        gap1_x_N = gap1 * n

        # Second gap if available
        if len(evals) >= 3:
            E2 = evals[2]
            gap2 = E2 - E0
            gap2_x_N = gap2 * n
            ratio_21 = gap2 / gap1
        else:
            gap2 = gap2_x_N = ratio_21 = None

        # E0/N for Casimir energy
        e0_per_site = E0 / n

        print(f"  E0 = {E0:.8f}, E0/N = {e0_per_site:.8f}")
        print(f"  Δ₁ = {gap1:.8f}, Δ₁·N = {gap1_x_N:.6f}")
        if gap2 is not None:
            print(f"  Δ₂ = {gap2:.8f}, Δ₂·N = {gap2_x_N:.6f}")
            print(f"  Δ₂/Δ₁ = {ratio_21:.4f}")
        print(f"  Spectrum: {evals[:6]}")
        print(f"  Time: {dt:.1f}s")

        results["sizes"][str(n)] = {
            "dim": dim,
            "gpu": use_gpu,
            "E0": float(E0),
            "E0_per_site": float(e0_per_site),
            "gap1": float(gap1),
            "gap1_x_N": float(gap1_x_N),
            "gap2": float(gap2) if gap2 is not None else None,
            "gap2_x_N": float(gap2_x_N) if gap2_x_N is not None else None,
            "ratio_21": float(ratio_21) if ratio_21 is not None else None,
            "evals": [float(e) for e in evals[:6]],
            "time_s": dt,
        }
    except Exception as e:
        dt = time.time() - t0
        print(f"  FAILED: {e} ({dt:.1f}s)")
        results["sizes"][str(n)] = {"error": str(e), "time_s": dt}

# ============================================================
# ANALYSIS: Convergence vs divergence of Δ·N
# ============================================================

print(f"\n{'='*60}")
print("ANALYSIS: Δ·N trend")
print(f"{'='*60}")

ns = []
gap_x_N = []
ratio_21_vals = []

for n_str, data in sorted(results["sizes"].items(), key=lambda x: int(x[0])):
    if "error" not in data:
        n = int(n_str)
        ns.append(n)
        gap_x_N.append(data["gap1_x_N"])
        if data.get("ratio_21") is not None:
            ratio_21_vals.append((n, data["ratio_21"]))

print(f"\n{'n':>4} {'Δ₁·N':>10} {'Δ₂/Δ₁':>10}")
print("-" * 26)
for i, n in enumerate(ns):
    r21_str = ""
    for nn, r in ratio_21_vals:
        if nn == n:
            r21_str = f"{r:.4f}"
    print(f"{n:>4} {gap_x_N[i]:>10.6f} {r21_str:>10}")

# Check for acceleration: second differences
if len(gap_x_N) >= 3:
    print(f"\nSuccessive differences in Δ₁·N:")
    diffs = [gap_x_N[i+1] - gap_x_N[i] for i in range(len(gap_x_N)-1)]
    for i, d in enumerate(diffs):
        print(f"  n={ns[i]}→{ns[i+1]}: Δ(Δ·N) = {d:+.6f}")

    if len(diffs) >= 2:
        accel = diffs[-1] - diffs[-2]
        print(f"\n  Acceleration (2nd diff): {accel:+.6f}")
        if accel > 0:
            print("  → Δ·N ACCELERATING (suggests first-order)")
        else:
            print("  → Δ·N DECELERATING (suggests continuous)")

# Save
with open("results/sprint_066a_gap_scaling_q10.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nResults saved to results/sprint_066a_gap_scaling_q10.json")

# Record to DB
from db_utils import record
for n_str, data in results["sizes"].items():
    if "error" not in data:
        n = int(n_str)
        record(sprint=66, model='hybrid', q=q, n=n, quantity='gap1_x_N',
               value=data['gap1_x_N'], method='exact_diag_periodic',
               notes=f'at g_c={gc}')
        if data.get('ratio_21') is not None:
            record(sprint=66, model='hybrid', q=q, n=n, quantity='gap_ratio_21',
                   value=data['ratio_21'], method='exact_diag_periodic',
                   notes=f'Δ₂/Δ₁ at g_c={gc}')
print("Results recorded to DB.")
