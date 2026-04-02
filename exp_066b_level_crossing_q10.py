#!/usr/bin/env python3
"""Sprint 066b: Level crossing scan at q=10.

Diagnostics for first-order vs continuous:
1. Gap minimum structure: single smooth minimum (continuous) vs sharp dip (first-order)
2. Gap minimum value scaling: Δ_min ~ 1/N^z (continuous) vs ~ exp(-αN) (first-order)
3. E₀/N curvature: smooth (continuous) vs kink (first-order)
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

print(f"Sprint 066b: Level crossing scan at q={q}")
print(f"{'='*60}")

# Scan range: g from 0.5 to 0.9 (centered on gc=0.684)
g_values = np.linspace(0.50, 0.90, 41)  # 41 points, step=0.01
sizes = [4, 5, 6]

results = {"q": q, "gc": gc, "scan": {}}

for n in sizes:
    dim = q**n
    if dim > 1_500_000:
        # For n=6, use coarser grid near gc
        g_scan = np.concatenate([
            np.linspace(0.55, 0.63, 5),
            np.linspace(0.64, 0.74, 21),  # Fine near gc
            np.linspace(0.75, 0.85, 5)
        ])
    else:
        g_scan = g_values

    print(f"\nn={n}: dim={dim:,}, {len(g_scan)} g-points")
    t0 = time.time()

    scan_data = {"g": [], "E0": [], "E1": [], "E2": [], "gap": [], "gap_x_N": []}

    for gi, g in enumerate(g_scan):
        H = hybrid_hamiltonian_periodic(n, q, g)
        k_eig = min(4, dim - 2)
        evals, _ = eigsh(H, k=k_eig, which='SA')
        evals = np.sort(evals)

        gap = evals[1] - evals[0]
        scan_data["g"].append(float(g))
        scan_data["E0"].append(float(evals[0]))
        scan_data["E1"].append(float(evals[1]))
        scan_data["E2"].append(float(evals[2]) if len(evals) > 2 else None)
        scan_data["gap"].append(float(gap))
        scan_data["gap_x_N"].append(float(gap * n))

    dt = time.time() - t0

    # Find gap minimum
    gaps = np.array(scan_data["gap"])
    gs = np.array(scan_data["g"])
    imin = np.argmin(gaps)
    g_min = gs[imin]
    gap_min = gaps[imin]

    print(f"  Gap minimum: Δ={gap_min:.6f} at g={g_min:.3f} (gc={gc})")
    print(f"  Δ_min·N = {gap_min*n:.6f}")
    print(f"  Time: {dt:.1f}s")

    # Check E0/N for kink near gc
    E0s = np.array(scan_data["E0"]) / n
    dE = np.gradient(E0s, gs)
    d2E = np.gradient(dE, gs)
    i_gc = np.argmin(np.abs(gs - gc))
    print(f"  d²(E₀/N)/dg² at g≈{gs[i_gc]:.3f}: {d2E[i_gc]:.4f}")

    scan_data["g_min"] = float(g_min)
    scan_data["gap_min"] = float(gap_min)
    scan_data["gap_min_x_N"] = float(gap_min * n)
    scan_data["time_s"] = dt
    results["scan"][str(n)] = scan_data

# ============================================================
# ANALYSIS: Gap minimum scaling
# ============================================================
print(f"\n{'='*60}")
print("ANALYSIS: Gap minimum scaling")
print(f"{'='*60}")

ns_data = []
gap_mins = []
g_mins = []

print(f"\n{'n':>4} {'g_min':>8} {'Δ_min':>12} {'Δ_min·N':>10}")
print("-" * 36)
for n_str in sorted(results["scan"].keys(), key=int):
    data = results["scan"][n_str]
    n = int(n_str)
    ns_data.append(n)
    gap_mins.append(data["gap_min"])
    g_mins.append(data["g_min"])
    print(f"{n:>4} {data['g_min']:>8.3f} {data['gap_min']:>12.6f} {data['gap_min_x_N']:>10.6f}")

# Power-law fit: Δ_min ~ N^{-z}
if len(ns_data) >= 2:
    log_n = np.log(ns_data)
    log_gap = np.log(gap_mins)
    z, c = np.polyfit(log_n, log_gap, 1)
    print(f"\nPower-law fit: Δ_min ~ N^{{{z:.3f}}}")
    print(f"  z = {-z:.3f} (dynamical exponent)")
    if -z > 0.5:
        print(f"  → Polynomial gap closing: CONTINUOUS transition")
    else:
        print(f"  → Slow gap closing: may be exponential (first-order)")

    # Check if exponential fits better
    # exp fit: log(Δ) = -α·N + b
    alpha, b_exp = np.polyfit(ns_data, log_gap, 1)
    resid_pow = np.sum((log_gap - (z * log_n + c))**2)
    resid_exp = np.sum((log_gap - (alpha * np.array(ns_data) + b_exp))**2)
    print(f"\n  Power-law residual: {resid_pow:.6f}")
    print(f"  Exponential residual: {resid_exp:.6f}")
    if resid_exp < resid_pow * 0.5:
        print(f"  → Exponential fits BETTER: suggests first-order")
    elif resid_pow < resid_exp * 0.5:
        print(f"  → Power-law fits BETTER: suggests continuous")
    else:
        print(f"  → Fits comparable: inconclusive at these sizes")

# Pseudo-critical point drift
print(f"\nPseudo-critical point (g_min) drift:")
for i in range(len(g_mins)):
    diff = g_mins[i] - gc
    print(f"  n={ns_data[i]}: g_min={g_mins[i]:.3f}, drift from gc: {diff:+.3f}")

# Save
with open("results/sprint_066b_level_crossing_q10.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to results/sprint_066b_level_crossing_q10.json")

from db_utils import record
for n_str, data in results["scan"].items():
    n = int(n_str)
    record(sprint=66, model='hybrid', q=q, n=n, quantity='gap_min',
           value=data['gap_min'], method='gap_scan',
           notes=f'min gap over g=[0.5,0.9], g_min={data["g_min"]:.3f}')
    record(sprint=66, model='hybrid', q=q, n=n, quantity='g_pseudo_critical',
           value=data['g_min'], method='gap_minimum',
           notes=f'g where gap is minimized')
print("Recorded to DB.")
