#!/usr/bin/env python3
"""Sprint 067a: Wide gap scan at q=5,7 looking for second critical point.

Scan energy gap across g=[0.01, 2.0] at multiple sizes (n=4,6,8 for q=5; n=4,6 for q=7).
Look for: (1) second gap minimum, (2) extended gapless region, (3) gap~1/N in any region.

Periodic BC for clean gap extraction.
"""
import numpy as np, json, time, sys
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh

def build_H_periodic(n, q, g):
    """Build H = -sum delta(s_i,s_{i+1}) - g*sum(X+X†) with periodic BC."""
    dim = q**n
    # Potts coupling: -delta(s_i, s_{i+1})
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

    # Nearest-neighbor Potts coupling (periodic)
    for i in range(n):
        j = (i + 1) % n
        if j == i + 1:
            # Consecutive sites
            left = q**i
            right = q**(n - i - 2)
            op = potts_op
            if left > 1:
                op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1:
                op = sp_kron(op, sp_eye(right), format='csr')
        else:
            # Boundary bond: sites n-1 and 0
            diag_boundary = np.zeros(dim)
            for idx in range(dim):
                s0 = idx % q
                sn = (idx // (q**(n-1))) % q
                diag_boundary[idx] = 1.0 if s0 == sn else 0.0
            op = diags(diag_boundary, 0, shape=(dim, dim), format='csr')
        H = H - op

    # Transverse field
    for i in range(n):
        left = q**i
        right = q**(n - i - 1)
        op = XpXd
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op

    return H

def scan_gaps(n, q, g_values):
    """Compute energy gap at each g value."""
    data = []
    for g in g_values:
        H = build_H_periodic(n, q, g)
        dim = q**n
        k = min(4, dim - 2)
        evals, _ = eigsh(H, k=k, which='SA')
        evals = np.sort(evals)
        gap = evals[1] - evals[0]
        gap2 = evals[2] - evals[0] if len(evals) > 2 else None
        data.append({
            "g": float(g), "E0": float(evals[0]), "E0_per_site": float(evals[0]/n),
            "gap1": float(gap), "gap1_x_N": float(gap * n),
            "gap2": float(gap2) if gap2 else None,
            "gap2_x_N": float(gap2 * n) if gap2 else None
        })
    return data

results = {"experiment": "067a_wide_gap_scan", "description": "Wide g scan looking for second critical point"}

# ============================================================
# q=5: n=4,6,8
# ============================================================
print("=" * 60)
print("q=5: Wide gap scan")
print("=" * 60)

g_values_q5 = np.concatenate([
    np.linspace(0.02, 0.20, 10),   # Deep ordered
    np.linspace(0.25, 0.60, 15),   # Around g_c=0.441
    np.linspace(0.65, 1.50, 18),   # Disordered side - looking for 2nd transition
    np.linspace(1.60, 3.00, 8),    # Far disordered
])
g_values_q5 = np.sort(np.unique(g_values_q5))

for n in [4, 6, 8]:
    dim = 5**n
    if dim > 500000 and n == 8:
        # q=5 n=8 = 390625, OK on CPU
        pass
    print(f"\nq=5, n={n} (dim={dim}): {len(g_values_q5)} points", flush=True)
    t0 = time.time()
    data = scan_gaps(n, 5, g_values_q5)
    dt = time.time() - t0
    print(f"  Done in {dt:.1f}s", flush=True)

    # Find gap minima
    gaps = [d["gap1"] for d in data]
    gs = [d["g"] for d in data]

    # Look for ALL local minima
    minima = []
    for i in range(1, len(gaps) - 1):
        if gaps[i] < gaps[i-1] and gaps[i] < gaps[i+1]:
            minima.append((gs[i], gaps[i], gaps[i]*n))

    print(f"  Local minima found: {len(minima)}")
    for gm, gapm, gapN in minima:
        print(f"    g={gm:.3f}: gap={gapm:.6f}, gap*N={gapN:.4f}")

    # Also report global minimum
    imin = np.argmin(gaps)
    print(f"  Global minimum: g={gs[imin]:.3f}, gap={gaps[imin]:.6f}, gap*N={gaps[imin]*n:.4f}")

    results[f"q5_n{n}"] = {"data": data, "time_s": dt, "minima": minima,
                            "global_min_g": float(gs[imin]), "global_min_gap": float(gaps[imin])}

# ============================================================
# q=7: n=4,6
# ============================================================
print("\n" + "=" * 60)
print("q=7: Wide gap scan")
print("=" * 60)

g_values_q7 = np.concatenate([
    np.linspace(0.02, 0.30, 8),
    np.linspace(0.35, 0.70, 15),   # Around g_c=0.535
    np.linspace(0.75, 1.50, 16),   # Disordered side
    np.linspace(1.60, 3.00, 8),
])
g_values_q7 = np.sort(np.unique(g_values_q7))

for n in [4, 6]:
    dim = 7**n
    print(f"\nq=7, n={n} (dim={dim}): {len(g_values_q7)} points", flush=True)
    t0 = time.time()
    data = scan_gaps(n, 7, g_values_q7)
    dt = time.time() - t0
    print(f"  Done in {dt:.1f}s", flush=True)

    gaps = [d["gap1"] for d in data]
    gs = [d["g"] for d in data]

    minima = []
    for i in range(1, len(gaps) - 1):
        if gaps[i] < gaps[i-1] and gaps[i] < gaps[i+1]:
            minima.append((gs[i], gaps[i], gaps[i]*n))

    print(f"  Local minima found: {len(minima)}")
    for gm, gapm, gapN in minima:
        print(f"    g={gm:.3f}: gap={gapm:.6f}, gap*N={gapN:.4f}")

    imin = np.argmin(gaps)
    print(f"  Global minimum: g={gs[imin]:.3f}, gap={gaps[imin]:.6f}, gap*N={gaps[imin]*n:.4f}")

    results[f"q7_n{n}"] = {"data": data, "time_s": dt, "minima": minima,
                            "global_min_g": float(gs[imin]), "global_min_gap": float(gaps[imin])}

# ============================================================
# ANALYSIS: Compare gap*N across sizes to detect gapless regions
# ============================================================
print("\n" + "=" * 60)
print("ANALYSIS: Gap*N comparison across sizes")
print("=" * 60)

for q_val, sizes in [(5, [4, 6, 8]), (7, [4, 6])]:
    print(f"\nq={q_val}:")
    # Find common g values
    ref_data = results[f"q{q_val}_n{sizes[0]}"]["data"]
    g_common = [d["g"] for d in ref_data]

    print(f"{'g':>6} ", end="")
    for n in sizes:
        print(f"{'gap*N(n='+str(n)+')':>14}", end="")
    if len(sizes) >= 2:
        print(f"{'ratio(last/first)':>18}", end="")
    print()

    for ig, g in enumerate(g_common):
        if ig % 3 != 0:
            continue  # Print every 3rd point
        vals = []
        for n in sizes:
            data = results[f"q{q_val}_n{n}"]["data"]
            vals.append(data[ig]["gap1_x_N"])
        print(f"{g:6.3f} ", end="")
        for v in vals:
            print(f"{v:14.6f}", end="")
        if len(vals) >= 2:
            ratio = vals[-1] / vals[0] if vals[0] > 0 else float('inf')
            print(f"{ratio:18.4f}", end="")
        print()

# Save
with open("results/sprint_067a_wide_gap_scan.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_067a_wide_gap_scan.json", flush=True)

# Record key findings to DB
from db_utils import record
for q_val in [5, 7]:
    for key in results:
        if key.startswith(f"q{q_val}_n"):
            n = int(key.split("n")[1])
            d = results[key]
            record(sprint=67, model='hybrid', q=q_val, n=n, quantity='n_gap_minima',
                   value=len(d['minima']), method='wide_scan',
                   notes=f'g_range=[0.02,3.0], {len(d["data"])} points')
            record(sprint=67, model='hybrid', q=q_val, n=n, quantity='global_gap_min_g',
                   value=d['global_min_g'], method='wide_scan')

print("Recorded to DB.", flush=True)
