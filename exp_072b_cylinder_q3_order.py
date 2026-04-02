#!/usr/bin/env python3
"""Sprint 072b: q=3 Ly=2 cylinder order parameter + cross-q comparison.

Compute <delta(s_i,s_j)> across transition for q=3 on Ly=2 cylinder.
Compare slope steepness with q=2 and q=5 (from Sprint 071c).
"""
import numpy as np
import json
import time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from gpu_utils import eigsh


def build_H_cylinder(Lx, Ly, q, g):
    """Build hybrid Hamiltonian on Lx x Ly cylinder."""
    n = Lx * Ly
    dim = q ** n
    potts_2site = np.zeros(q ** 2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q ** 2, q ** 2), format='csr')
    X = np.zeros((q, q))
    for s in range(q):
        X[(s + 1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)
    H = csr_matrix((dim, dim))
    bonds = set()
    for y in range(Ly):
        for x in range(Lx):
            site = y * Lx + x
            if x + 1 < Lx:
                nbr = y * Lx + (x + 1)
                bonds.add((min(site, nbr), max(site, nbr)))
            nbr_v = ((y + 1) % Ly) * Lx + x
            if nbr_v != site:
                bonds.add((min(site, nbr_v), max(site, nbr_v)))
    for (i, j) in bonds:
        if j == i + 1:
            left = q ** i
            right = q ** (n - j - 1)
            op = potts_op
            if left > 1:
                op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1:
                op = sp_kron(op, sp_eye(right), format='csr')
        else:
            indices = np.arange(dim)
            si = (indices // (q ** i)) % q
            sj = (indices // (q ** j)) % q
            diag_vals = (si == sj).astype(float)
            op = diags(diag_vals, 0, shape=(dim, dim), format='csr')
        H = H - op
    for i in range(n):
        left = q ** i
        right = q ** (n - i - 1)
        op = XpXd.copy()
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H


def compute_order_param(Lx, Ly, q, g):
    """Compute <delta(s_i, s_j)> for middle rung bond."""
    H = build_H_cylinder(Lx, Ly, q, g)
    dim = H.shape[0]
    n = Lx * Ly
    if dim < 500:
        evals, evecs = np.linalg.eigh(H.toarray())
        psi = evecs[:, 0]
    else:
        _, evecs = eigsh(H, k=1, which='SA')
        psi = evecs[:, 0]
    mid_x = Lx // 2
    i = 0 * Lx + mid_x
    j = 1 * Lx + mid_x
    indices = np.arange(dim)
    si = (indices // (q ** i)) % q
    sj = (indices // (q ** j)) % q
    return float(np.sum(psi ** 2 * (si == sj)))


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 072b: q=3 cylinder order parameter", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '072b_cylinder_q3_order',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': 3, 'Ly': 2,
}

gc_q3 = 0.565  # From 072a

# q=3 order parameter scan at Lx=5 (dim=59049)
Lx = 5
g_range = np.linspace(0.2, 1.0, 33)
print(f"\n--- q=3, Lx={Lx} order parameter scan ---", flush=True)
t0 = time.time()
op_data = []
for gi, g in enumerate(g_range):
    op = compute_order_param(Lx, 2, 3, g)
    op_data.append({'g': round(float(g), 4), 'delta_ij': float(op)})
    if (gi + 1) % 11 == 0:
        print(f"  {gi+1}/{len(g_range)} done", flush=True)
dt = time.time() - t0
print(f"  Done in {dt:.1f}s", flush=True)

results['q3_order_param'] = {'Lx': Lx, 'data': op_data}

# Compute slope
gs = np.array([d['g'] for d in op_data])
ops = np.array([d['delta_ij'] for d in op_data])
dop_dg = np.gradient(ops, gs)
max_slope = float(np.max(np.abs(dop_dg)))
max_slope_g = float(gs[np.argmax(np.abs(dop_dg))])
op_at_gc = float(np.interp(gc_q3, gs, ops))

print(f"\n  <delta> at g_c={gc_q3:.3f}: {op_at_gc:.4f} (random=0.333, ordered=1)", flush=True)
print(f"  Max |d<delta>/dg| = {max_slope:.3f} at g={max_slope_g:.3f}", flush=True)

results['q3_analysis'] = {
    'gc': gc_q3,
    'op_at_gc': op_at_gc,
    'max_slope': max_slope,
    'max_slope_g': max_slope_g,
}

# Also scan at Lx=6 for size comparison
print(f"\n--- q=3, Lx=6 order parameter (fewer points) ---", flush=True)
g_range_6 = np.linspace(0.3, 0.9, 13)
t0 = time.time()
op_data_6 = []
for g in g_range_6:
    op = compute_order_param(6, 2, 3, g)
    op_data_6.append({'g': round(float(g), 4), 'delta_ij': float(op)})
dt = time.time() - t0
print(f"  Done in {dt:.1f}s", flush=True)
results['q3_order_param_Lx6'] = {'Lx': 6, 'data': op_data_6}

gs6 = np.array([d['g'] for d in op_data_6])
ops6 = np.array([d['delta_ij'] for d in op_data_6])
dop6 = np.gradient(ops6, gs6)
max_slope_6 = float(np.max(np.abs(dop6)))
print(f"  Max |d<delta>/dg| = {max_slope_6:.3f} (Lx=5: {max_slope:.3f})", flush=True)

results['q3_slope_comparison'] = {
    'Lx5': max_slope, 'Lx6': max_slope_6,
    'ratio': max_slope_6 / max_slope if max_slope > 0 else None,
}

# Cross-q comparison
print(f"\n{'='*60}", flush=True)
print("Cross-q cylinder comparison (all Ly=2)", flush=True)
print(f"{'='*60}", flush=True)

# q=2 data from Sprint 071
gc_q2_cyl = 0.451
cyl_1d_q2 = gc_q2_cyl / 0.250
# q=3 from this sprint
gc_q3_cyl = gc_q3
cyl_1d_q3 = gc_q3_cyl / (1.0/3.0)
# q=5 data from Sprint 071
gc_q5_cyl = 0.714
cyl_1d_q5 = gc_q5_cyl / 0.441

# 2D values
gc_q2_2d = 0.771
gc_q3_2d = 1.267
gc_q5_2d = 1.588

print(f"\n{'q':>3} {'g_c(1D)':>8} {'g_c(cyl)':>9} {'g_c(2D)':>8} {'cyl/1D':>7} {'2D/cyl':>7} {'2D/1D':>6}", flush=True)
print("-" * 55, flush=True)
for q, gc1d, gc_cyl, gc2d in [(2, 0.250, gc_q2_cyl, gc_q2_2d),
                                (3, 1/3, gc_q3_cyl, gc_q3_2d),
                                (5, 0.441, gc_q5_cyl, gc_q5_2d)]:
    r_cyl = gc_cyl / gc1d
    r_2d_cyl = gc2d / gc_cyl
    r_2d = gc2d / gc1d
    print(f"{q:>3} {gc1d:>8.3f} {gc_cyl:>9.3f} {gc2d:>8.3f} {r_cyl:>7.3f} {r_2d_cyl:>7.3f} {r_2d:>6.2f}", flush=True)

results['cross_q_table'] = {
    'q2': {'gc_1d': 0.250, 'gc_cyl': gc_q2_cyl, 'gc_2d': gc_q2_2d,
           'cyl_1d': cyl_1d_q2, '2d_cyl': gc_q2_2d/gc_q2_cyl},
    'q3': {'gc_1d': 1/3, 'gc_cyl': gc_q3_cyl, 'gc_2d': gc_q3_2d,
           'cyl_1d': cyl_1d_q3, '2d_cyl': gc_q3_2d/gc_q3_cyl},
    'q5': {'gc_1d': 0.441, 'gc_cyl': gc_q5_cyl, 'gc_2d': gc_q5_2d,
           'cyl_1d': cyl_1d_q5, '2d_cyl': gc_q5_2d/gc_q5_cyl},
}

# Order param slope comparison
# Sprint 071c: q=2 max_slope=1.21 (Lx=6), q=5 max_slope=1.88 (Lx=4)
print(f"\nOrder parameter max slope:", flush=True)
print(f"  q=2: 1.21 (Lx=6, Sprint 071c)", flush=True)
print(f"  q=3: {max_slope:.2f} (Lx=5), {max_slope_6:.2f} (Lx=6)", flush=True)
print(f"  q=5: 1.88 (Lx=4, Sprint 071c)", flush=True)

print(f"\nCyl/1D ratio trend: {cyl_1d_q2:.3f} → {cyl_1d_q3:.3f} → {cyl_1d_q5:.3f}", flush=True)
print(f"  Pattern: q=2 highest (1.80), monotonically decreasing to q=5 (1.62)", flush=True)

# Save
with open("results/sprint_072b_cylinder_q3_order.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_072b_cylinder_q3_order.json", flush=True)

from db_utils import record
record(sprint=72, model='hybrid_cyl', q=3, n=2,
       quantity='op_max_slope', value=max_slope,
       method='order_param_Ly2', notes=f'Lx={Lx} cylinder')
