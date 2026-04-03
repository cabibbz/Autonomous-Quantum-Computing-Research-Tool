#!/usr/bin/env python3
"""Sprint 066c: Cross-q comparison of Δ·N trend.

Compare Δ·N at g_c for q=3,5,7,10 at multiple sizes.
q=3 is the control (known continuous, exact ν=5/6).
If q=10 anomalous FSS is just a correction, it should match q=3 pattern (decelerate).
If it's first-order, it should accelerate (diverge).
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh
import json, time

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

# Critical points from KNOWLEDGE.md
q_data = {
    3:  {"gc": 0.333, "sizes": [4, 6, 8, 10, 12]},  # dim=12^q... check limits
    5:  {"gc": 0.441, "sizes": [4, 6, 8]},
    7:  {"gc": 0.535, "sizes": [4, 6]},
    10: {"gc": 0.684, "sizes": [4, 5, 6]},
}

# Enforce dim limit ~1M
for q_val, info in q_data.items():
    info["sizes"] = [n for n in info["sizes"] if q_val**n <= 1_500_000]

print(f"Sprint 066c: Cross-q Δ·N comparison")
print(f"{'='*60}")

results = {}

for q_val, info in sorted(q_data.items()):
    gc = info["gc"]
    sizes = info["sizes"]
    print(f"\nq={q_val}, gc={gc}, sizes={sizes}")

    q_results = {"gc": gc, "sizes": {}}

    for n in sizes:
        dim = q_val**n
        t0 = time.time()
        H = hybrid_hamiltonian_periodic(n, q_val, gc)
        k_eig = min(4, dim - 2)
        evals, _ = eigsh(H, k=k_eig, which='SA')
        evals = np.sort(evals)
        dt = time.time() - t0

        gap = evals[1] - evals[0]
        gap_x_N = gap * n
        print(f"  n={n} (dim={dim:,}): Δ·N = {gap_x_N:.6f} ({dt:.1f}s)")

        q_results["sizes"][str(n)] = {
            "dim": dim, "gap": float(gap), "gap_x_N": float(gap_x_N),
            "E0": float(evals[0]), "time_s": dt
        }

    results[str(q_val)] = q_results

# ============================================================
# ANALYSIS: trend comparison
# ============================================================
print(f"\n{'='*60}")
print("CROSS-q COMPARISON")
print(f"{'='*60}")

for q_str, q_res in sorted(results.items(), key=lambda x: int(x[0])):
    q_val = int(q_str)
    ns = sorted([int(k) for k in q_res["sizes"].keys()])
    gaps = [q_res["sizes"][str(n)]["gap_x_N"] for n in ns]

    print(f"\nq={q_val} (gc={q_res['gc']}):")
    for i, n in enumerate(ns):
        print(f"  n={n}: Δ·N = {gaps[i]:.6f}")

    if len(gaps) >= 2:
        diffs = [gaps[i+1] - gaps[i] for i in range(len(gaps)-1)]
        for i, d in enumerate(diffs):
            print(f"  Δ(Δ·N) n={ns[i]}→{ns[i+1]}: {d:+.6f}")

        if len(diffs) >= 2:
            accel = diffs[-1] - diffs[-2]
            if accel > 0:
                trend = "ACCELERATING"
            elif accel < -abs(diffs[0]) * 0.1:
                trend = "DECELERATING"
            else:
                trend = "ROUGHLY LINEAR"
            print(f"  Trend: {trend} (accel={accel:+.6f})")

# Normalize: Δ·N / Δ·N(n_min) to compare shapes
print(f"\n{'='*60}")
print("NORMALIZED Δ·N (relative to n=4)")
print(f"{'='*60}")

print(f"\n{'q':>3} | ", end="")
all_ns = sorted(set(int(n) for q_res in results.values() for n in q_res["sizes"]))
for n in all_ns:
    print(f"{'n='+str(n):>10}", end="")
print()
print("-" * (5 + 10*len(all_ns)))

for q_str in sorted(results.keys(), key=int):
    q_res = results[q_str]
    # Normalize to n=4 if available
    base = q_res["sizes"].get("4", {}).get("gap_x_N", None)
    print(f"{q_str:>3} | ", end="")
    for n in all_ns:
        val = q_res["sizes"].get(str(n), {}).get("gap_x_N", None)
        if val is not None and base is not None:
            print(f"{val/base:>10.4f}", end="")
        elif val is not None:
            print(f"{'—':>10}", end="")
        else:
            print(f"{'':>10}", end="")
    print()

# Save
with open("results/sprint_066c_cross_q_gap_trend.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to results/sprint_066c_cross_q_gap_trend.json")

from db_utils import record
for q_str, q_res in results.items():
    q_val = int(q_str)
    for n_str, data in q_res["sizes"].items():
        n = int(n_str)
        record(sprint=66, model='hybrid', q=q_val, n=n, quantity='gap1_x_N',
               value=data['gap_x_N'], method='exact_diag_periodic',
               notes=f'at g_c={q_res["gc"]}')
print("Recorded to DB.")
