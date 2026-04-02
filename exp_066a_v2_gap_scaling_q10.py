#!/usr/bin/env python3
"""Sprint 066a v2: Δ·N scaling at q=10 for n=4,5,6 (CPU).

CuPy unavailable, so use scipy eigsh for all sizes.
n=6 dim=10^6 should work in ~60s with sparse eigsh.
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh
import json, time

q = 10
gc = 0.684

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

print(f"Sprint 066a v2: Δ·N scaling at q={q}, g_c={gc}")
print(f"{'='*60}")

sizes = [4, 5, 6]
results = {"q": q, "gc": gc, "sizes": {}}

for n in sizes:
    dim = q**n
    k_eig = min(6, dim - 2)
    print(f"\nn={n}: dim={dim:,}, k={k_eig}")
    t0 = time.time()
    try:
        H = hybrid_hamiltonian_periodic(n, q, gc)
        evals, _ = eigsh(H, k=k_eig, which='SA')
        evals = np.sort(evals)
        dt = time.time() - t0

        E0, E1 = evals[0], evals[1]
        gap1 = E1 - E0
        gap1_x_N = gap1 * n
        E2 = evals[2] if len(evals) >= 3 else None
        gap2 = (E2 - E0) if E2 is not None else None
        gap2_x_N = gap2 * n if gap2 is not None else None
        ratio_21 = gap2 / gap1 if gap2 is not None else None

        print(f"  E0 = {E0:.8f}, E0/N = {E0/n:.8f}")
        print(f"  Δ₁ = {gap1:.8f}, Δ₁·N = {gap1_x_N:.6f}")
        if gap2 is not None:
            print(f"  Δ₂ = {gap2:.8f}, Δ₂·N = {gap2_x_N:.6f}")
            print(f"  Δ₂/Δ₁ = {ratio_21:.4f}")
        print(f"  Spectrum: {evals[:6]}")
        print(f"  Time: {dt:.1f}s")

        results["sizes"][str(n)] = {
            "dim": dim, "E0": float(E0), "E0_per_site": float(E0/n),
            "gap1": float(gap1), "gap1_x_N": float(gap1_x_N),
            "gap2": float(gap2) if gap2 else None,
            "gap2_x_N": float(gap2_x_N) if gap2_x_N else None,
            "ratio_21": float(ratio_21) if ratio_21 else None,
            "evals": [float(e) for e in evals[:6]], "time_s": dt,
        }
    except Exception as e:
        dt = time.time() - t0
        print(f"  FAILED: {e} ({dt:.1f}s)")
        results["sizes"][str(n)] = {"error": str(e), "time_s": dt}

# ============================================================
# ANALYSIS
# ============================================================
print(f"\n{'='*60}")
print("ANALYSIS: Δ·N trend at q=10")
print(f"{'='*60}")

ns, gap_x_N = [], []
print(f"\n{'n':>4} {'dim':>12} {'Δ₁·N':>10} {'Δ₂/Δ₁':>10}")
print("-" * 40)
for n_str in sorted(results["sizes"].keys(), key=int):
    data = results["sizes"][n_str]
    if "error" not in data:
        n = int(n_str)
        ns.append(n)
        gap_x_N.append(data["gap1_x_N"])
        r21 = f"{data['ratio_21']:.4f}" if data.get('ratio_21') else "—"
        print(f"{n:>4} {data['dim']:>12,} {data['gap1_x_N']:>10.6f} {r21:>10}")

if len(gap_x_N) >= 2:
    print(f"\nSuccessive changes in Δ₁·N:")
    for i in range(len(gap_x_N)-1):
        d = gap_x_N[i+1] - gap_x_N[i]
        pct = d / gap_x_N[i] * 100
        print(f"  n={ns[i]}→{ns[i+1]}: Δ(Δ·N) = {d:+.6f} ({pct:+.2f}%)")

    if len(gap_x_N) >= 3:
        d1 = gap_x_N[1] - gap_x_N[0]
        d2 = gap_x_N[2] - gap_x_N[1]
        accel = d2 - d1
        print(f"\n  1st diff: {d1:+.6f}, {d2:+.6f}")
        print(f"  Acceleration: {accel:+.6f}")
        if accel > 0:
            print("  → ACCELERATING: consistent with weakly first-order")
        elif abs(accel) < abs(d1) * 0.1:
            print("  → ROUGHLY LINEAR: inconclusive")
        else:
            print("  → DECELERATING: consistent with continuous transition")

# Save
with open("results/sprint_066a_gap_scaling_q10.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to results/sprint_066a_gap_scaling_q10.json")

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
print("Recorded to DB.")
