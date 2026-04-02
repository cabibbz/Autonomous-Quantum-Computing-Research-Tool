#!/usr/bin/env python3
"""Sprint 067c: Clock model gap scan at q=5,7 — validate floating phase detection.

If the clock model shows two gap minima or an extended gapless region where the
hybrid model doesn't, it confirms the hybrid's single-transition is real.

Clock Hamiltonian: H = -sum cos(2pi(s_i-s_j)/q) - g*sum(X+X†)
Uses periodic BC for clean gap extraction.
"""
import numpy as np, json, time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh

def build_H_clock_periodic(n, q, g):
    """Clock model with cosine coupling + clock field, periodic BC."""
    dim = q**n
    omega = np.exp(2j * np.pi / q)

    # Clock coupling: -cos(2pi(s_i-s_j)/q) = -Re(Z_i Z_j†)
    # = -(1/2)(Z_i Z_j† + Z_i† Z_j)
    # Build 2-site coupling operator
    clock_2site = np.zeros((q**2, q**2))
    for si in range(q):
        for sj in range(q):
            idx = si * q + sj
            clock_2site[idx, idx] = np.cos(2 * np.pi * (si - sj) / q)
    clock_op = diags(clock_2site.flatten()[:q**2], 0, shape=(q**2, q**2), format='csr')
    # Reconstruct properly
    clock_diag = np.zeros(q**2)
    for si in range(q):
        for sj in range(q):
            clock_diag[si * q + sj] = np.cos(2 * np.pi * (si - sj) / q)
    clock_op = diags(clock_diag, 0, shape=(q**2, q**2), format='csr')

    # Transverse field: X + X†
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)

    H = csr_matrix((dim, dim))

    # Clock coupling (periodic)
    for i in range(n):
        j = (i + 1) % n
        if j == i + 1:
            left = q**i; right = q**(n - i - 2)
            op = clock_op
            if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        else:
            # Boundary: sites n-1 and 0
            diag_bd = np.zeros(dim)
            for idx in range(dim):
                s0 = idx % q
                sn = (idx // (q**(n-1))) % q
                diag_bd[idx] = np.cos(2 * np.pi * (sn - s0) / q)
            op = diags(diag_bd, 0, shape=(dim, dim), format='csr')
        H = H - op

    # Transverse field
    for i in range(n):
        left = q**i; right = q**(n - i - 1)
        op = XpXd
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op

    return H

def build_H_hybrid_periodic(n, q, g):
    """Hybrid model (Potts coupling + clock field), periodic BC."""
    dim = q**n
    potts_2site = np.zeros(q**2)
    for a in range(q): potts_2site[a*q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q**2, q**2), format='csr')
    X = np.zeros((q, q))
    for s in range(q): X[(s+1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)
    H = csr_matrix((dim, dim))
    for i in range(n):
        j = (i + 1) % n
        if j == i + 1:
            left = q**i; right = q**(n - i - 2)
            op = potts_op
            if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        else:
            diag_bd = np.zeros(dim)
            for idx in range(dim):
                s0 = idx % q; sn = (idx // (q**(n-1))) % q
                diag_bd[idx] = 1.0 if s0 == sn else 0.0
            op = diags(diag_bd, 0, shape=(dim, dim), format='csr')
        H = H - op
    for i in range(n):
        left = q**i; right = q**(n - i - 1)
        op = XpXd
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H

def scan_model(builder, n, q, g_values, label):
    """Scan gap for a model."""
    data = []
    for g in g_values:
        H = builder(n, q, g)
        dim = q**n
        k = min(4, dim - 2)
        evals, _ = eigsh(H, k=k, which='SA')
        evals = np.sort(evals)
        gap = evals[1] - evals[0]
        data.append({"g": float(g), "gap": float(gap), "gap_x_N": float(gap*n),
                     "E0": float(evals[0])})
    return data

# ============================================================
# q=5: Compare clock vs hybrid gap profiles
# ============================================================
q = 5
g_values = np.concatenate([
    np.linspace(0.05, 0.30, 6),
    np.linspace(0.35, 0.70, 15),
    np.linspace(0.75, 1.50, 10),
    np.linspace(1.60, 2.50, 5),
])
g_values = np.sort(np.unique(g_values))

results = {"experiment": "067c_clock_comparison", "q": q, "g_values": list(map(float, g_values))}

for n in [4, 6, 8]:
    dim = q**n
    if dim > 500000:
        continue  # q=5 n=8 = 390625, OK

    print(f"\n{'='*60}")
    print(f"q={q}, n={n} (dim={dim})")
    print(f"{'='*60}")

    print(f"  Clock model: {len(g_values)} points...", flush=True)
    t0 = time.time()
    clock_data = scan_model(build_H_clock_periodic, n, q, g_values, "clock")
    dt_c = time.time() - t0
    print(f"  Done in {dt_c:.1f}s", flush=True)

    print(f"  Hybrid model: {len(g_values)} points...", flush=True)
    t0 = time.time()
    hybrid_data = scan_model(build_H_hybrid_periodic, n, q, g_values, "hybrid")
    dt_h = time.time() - t0
    print(f"  Done in {dt_h:.1f}s", flush=True)

    # Find local gap minima for both models
    for model, data, name in [(clock_data, clock_data, "Clock"), (hybrid_data, hybrid_data, "Hybrid")]:
        gaps = [d["gap"] for d in data]
        gs = [d["g"] for d in data]
        minima = []
        for i in range(1, len(gaps)-1):
            if gaps[i] < gaps[i-1] and gaps[i] < gaps[i+1]:
                minima.append((gs[i], gaps[i]))
        print(f"\n  {name} local minima: {len(minima)}")
        for gm, gapm in minima:
            print(f"    g={gm:.3f}: gap={gapm:.6f}, gap*N={gapm*n:.4f}")
        if not minima:
            imin = np.argmin(gaps)
            print(f"    Global min: g={gs[imin]:.3f}, gap={gaps[imin]:.6f}")

    results[f"clock_n{n}"] = clock_data
    results[f"hybrid_n{n}"] = hybrid_data

# ============================================================
# Compare gap*N profiles at n=4,6 for both models
# ============================================================
print("\n" + "=" * 60)
print("Gap*N comparison: Clock vs Hybrid at q=5")
print("=" * 60)

print(f"\n{'g':>6} {'Clock(4)':>10} {'Clock(6)':>10} {'Hybrid(4)':>10} {'Hybrid(6)':>10} {'C_ratio':>8} {'H_ratio':>8}")
print("-" * 70)
for ig, g in enumerate(g_values):
    if ig % 2 != 0:
        continue
    c4 = results["clock_n4"][ig]["gap_x_N"]
    c6 = results["clock_n6"][ig]["gap_x_N"]
    h4 = results["hybrid_n4"][ig]["gap_x_N"]
    h6 = results["hybrid_n6"][ig]["gap_x_N"]
    c_rat = c6/c4 if c4 > 1e-8 else 0
    h_rat = h6/h4 if h4 > 1e-8 else 0
    flag = ""
    # Flag regions where gap*N ratio is near 1 (gapless) for clock but not hybrid
    if 0.8 < c_rat < 1.2 and not (0.8 < h_rat < 1.2):
        flag = " <-- CLOCK GAPLESS, HYBRID NOT"
    elif 0.8 < c_rat < 1.2 and 0.8 < h_rat < 1.2:
        flag = " <-- BOTH GAPLESS"
    print(f"{g:6.3f} {c4:10.4f} {c6:10.4f} {h4:10.4f} {h6:10.4f} {c_rat:8.3f} {h_rat:8.3f}{flag}")

# ============================================================
# Key question: Does clock show extended gapless region?
# ============================================================
print("\n" + "=" * 60)
print("Gapless region detection (gap*N ratio ∈ [0.7, 1.3])")
print("=" * 60)

for model_name, prefix in [("Clock", "clock"), ("Hybrid", "hybrid")]:
    gapless_region = []
    for ig, g in enumerate(g_values):
        gN4 = results[f"{prefix}_n4"][ig]["gap_x_N"]
        gN6 = results[f"{prefix}_n6"][ig]["gap_x_N"]
        if gN4 > 0.01:  # Not ordered phase
            ratio = gN6 / gN4
            if 0.7 < ratio < 1.3:
                gapless_region.append(g)
    if gapless_region:
        print(f"  {model_name}: gapless g ∈ [{min(gapless_region):.3f}, {max(gapless_region):.3f}]")
    else:
        print(f"  {model_name}: no gapless region detected")

# Save
with open("results/sprint_067c_clock_comparison.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to results/sprint_067c_clock_comparison.json", flush=True)

from db_utils import record
# Record clock gap minima
for n in [4, 6]:
    key = f"clock_n{n}"
    if key in results:
        gaps = [d["gap"] for d in results[key]]
        gs = [d["g"] for d in results[key]]
        imin = np.argmin(gaps)
        record(sprint=67, model='clock', q=q, n=n, quantity='gap_min_g',
               value=gs[imin], method='wide_scan')
        record(sprint=67, model='clock', q=q, n=n, quantity='gap_min',
               value=gaps[imin], method='wide_scan')

print("Recorded to DB.", flush=True)
