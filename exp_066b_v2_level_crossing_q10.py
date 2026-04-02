#!/usr/bin/env python3
"""Sprint 066b v2: Level crossing scan at q=10 (lean version).

Fewer g-points, flush output, focus on gap minimum scaling.
"""
import numpy as np, sys
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh
import json, time

q = 10; gc = 0.684

def build_H(n, q, g):
    dim = q**n
    potts_2site = np.zeros(q**2)
    for a in range(q): potts_2site[a*q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q**2, q**2), format='csr')
    X = np.zeros((q, q))
    for s in range(q): X[(s+1) % q, s] = 1.0
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

def scan_gaps(n, q, g_values):
    data = []
    for g in g_values:
        H = build_H(n, q, g)
        k = min(4, q**n - 2)
        evals, _ = eigsh(H, k=k, which='SA')
        evals = np.sort(evals)
        gap = evals[1] - evals[0]
        data.append({"g": float(g), "E0": float(evals[0]), "gap": float(gap),
                      "gap_x_N": float(gap*n)})
    return data

results = {"q": q, "gc": gc, "scans": {}}

# n=4: full scan, fast
print(f"n=4: 21 points", flush=True)
t0 = time.time()
g4 = np.linspace(0.50, 0.90, 21)
d4 = scan_gaps(4, q, g4)
dt = time.time() - t0
gaps4 = [d["gap"] for d in d4]
imin4 = np.argmin(gaps4)
print(f"  gap_min={gaps4[imin4]:.6f} at g={g4[imin4]:.3f} ({dt:.1f}s)", flush=True)
results["scans"]["4"] = {"data": d4, "g_min": float(g4[imin4]), "gap_min": float(gaps4[imin4])}

# n=5: medium scan
print(f"n=5: 15 points", flush=True)
t0 = time.time()
g5 = np.linspace(0.55, 0.85, 15)
d5 = scan_gaps(5, q, g5)
dt = time.time() - t0
gaps5 = [d["gap"] for d in d5]
imin5 = np.argmin(gaps5)
print(f"  gap_min={gaps5[imin5]:.6f} at g={g5[imin5]:.3f} ({dt:.1f}s)", flush=True)
results["scans"]["5"] = {"data": d5, "g_min": float(g5[imin5]), "gap_min": float(gaps5[imin5])}

# n=6: sparse scan near gc
print(f"n=6: 11 points near gc", flush=True)
t0 = time.time()
g6 = np.linspace(0.60, 0.78, 11)
d6 = scan_gaps(6, q, g6)
dt = time.time() - t0
gaps6 = [d["gap"] for d in d6]
imin6 = np.argmin(gaps6)
print(f"  gap_min={gaps6[imin6]:.6f} at g={g6[imin6]:.3f} ({dt:.1f}s)", flush=True)
results["scans"]["6"] = {"data": d6, "g_min": float(g6[imin6]), "gap_min": float(gaps6[imin6])}

# ============================================================
# ANALYSIS
# ============================================================
print(f"\n{'='*60}", flush=True)
print("Gap minimum scaling", flush=True)
ns = [4, 5, 6]
gap_mins = [results["scans"][str(n)]["gap_min"] for n in ns]
g_mins = [results["scans"][str(n)]["g_min"] for n in ns]

print(f"\n{'n':>4} {'g_min':>8} {'Δ_min':>12} {'Δ_min·N':>10} {'ln(Δ)':>10}")
print("-" * 44)
for i, n in enumerate(ns):
    print(f"{n:>4} {g_mins[i]:>8.3f} {gap_mins[i]:>12.6f} {gap_mins[i]*n:>10.6f} {np.log(gap_mins[i]):>10.4f}")

# Power-law fit: Δ_min ~ N^{-z}
log_n = np.log(ns)
log_gap = np.log(gap_mins)
z, c = np.polyfit(log_n, log_gap, 1)
print(f"\nPower law: Δ_min ~ N^{{{z:.3f}}}  (z={-z:.3f})")

# Exponential fit: log(Δ) = -α·N + b
alpha, b = np.polyfit(ns, log_gap, 1)
resid_pow = np.sum((log_gap - (z * log_n + c))**2)
resid_exp = np.sum((log_gap - (alpha * np.array(ns) + b))**2)
print(f"Exp fit: Δ_min ~ exp({alpha:.4f}·N)  (α={-alpha:.4f})")
print(f"  Power residual: {resid_pow:.8f}")
print(f"  Exp residual: {resid_exp:.8f}")

if resid_exp < resid_pow * 0.5:
    print("  → Exponential BETTER: first-order signature")
elif resid_pow < resid_exp * 0.5:
    print("  → Power-law BETTER: continuous signature")
else:
    print("  → Comparable: inconclusive")

# Pseudo-critical drift
print(f"\nPseudo-critical drift: g_min converges toward gc={gc}")
for i, n in enumerate(ns):
    print(f"  n={n}: g_min={g_mins[i]:.3f}, offset={g_mins[i]-gc:+.3f}")

# Save
with open("results/sprint_066b_level_crossing_q10.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to results/sprint_066b_level_crossing_q10.json", flush=True)

from db_utils import record
for n_str in results["scans"]:
    n = int(n_str)
    d = results["scans"][n_str]
    record(sprint=66, model='hybrid', q=q, n=n, quantity='gap_min',
           value=d['gap_min'], method='gap_scan',
           notes=f'g_min={d["g_min"]:.3f}')
    record(sprint=66, model='hybrid', q=q, n=n, quantity='g_pseudo_critical',
           value=d['g_min'], method='gap_minimum')
print("Recorded to DB.", flush=True)
