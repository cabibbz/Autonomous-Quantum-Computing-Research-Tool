#!/usr/bin/env python3
"""Sprint 053a: Energy gap data collapse for ν(q=3) at g_c=1/3.

At criticality, Δ·N = f((g - g_c)·N^{1/ν}).
With known g_c = 1/3, optimize ν to get best collapse of Δ·N curves.
Exact answer: ν = 5/6 ≈ 0.8333.

Use exact diag for N=4,6,8 (q=3: dims 81, 729, 6561).
"""
import numpy as np, json, time, warnings
warnings.filterwarnings('ignore')
from scipy.sparse.linalg import eigsh
from scipy.sparse import kron, eye, csr_matrix
from scipy.interpolate import interp1d

def potts_hamiltonian(q, n, g, J=1.0):
    """H = -J Σ δ(s_i,s_j) - g Σ (X + X†)"""
    dim = q**n
    X = np.zeros((q, q), dtype=complex)
    for a in range(q):
        X[(a+1) % q, a] = 1.0
    Xphc = X + X.conj().T

    projectors = [csr_matrix(np.diag([1.0 if b == a else 0.0 for b in range(q)])) for a in range(q)]
    H = csr_matrix((dim, dim), dtype=complex)
    I = csr_matrix(np.eye(q))

    for i in range(n - 1):
        for a in range(q):
            op = csr_matrix(np.eye(1))
            for j in range(n):
                if j == i:
                    op = kron(op, projectors[a])
                elif j == i + 1:
                    op = kron(op, projectors[a])
                else:
                    op = kron(op, I)
            H += -J * op

    Xphc_sp = csr_matrix(Xphc)
    for i in range(n):
        op = csr_matrix(np.eye(1))
        for j in range(n):
            if j == i:
                op = kron(op, Xphc_sp)
            else:
                op = kron(op, I)
        H += -g * op

    return H

def energy_gap(q, n, g, J=1.0):
    H = potts_hamiltonian(q, n, g, J)
    vals = eigsh(H.real if q == 2 else H, k=min(4, H.shape[0]-1), which='SA', return_eigenvectors=False)
    vals = np.sort(vals.real)
    return vals[1] - vals[0]

# === Timing test ===
print("=== Timing test ===", flush=True)
t0 = time.time()
gap = energy_gap(3, 8, 0.333)
print(f"q=3 n=8: gap={gap:.6f}, t={time.time()-t0:.1f}s", flush=True)

# === Compute gaps for N=4,6,8 near g_c=1/3 ===
g_c = 1.0 / 3.0
# Fine grid near g_c, wider grid farther away
g_inner = np.linspace(0.20, 0.48, 57)  # fine near g_c
g_values = np.sort(np.unique(np.round(g_inner, 4)))

results = {}
for n in [4, 6, 8]:
    dim = 3**n
    print(f"\n=== q=3 n={n} (dim={dim}) ===", flush=True)
    t0 = time.time()
    gaps = []
    for g in g_values:
        gap = energy_gap(3, n, g)
        gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * n)})
    dt = time.time() - t0
    print(f"  Done: {len(g_values)} points in {dt:.1f}s", flush=True)
    results[f'n{n}'] = gaps

# Save immediately
with open('results/sprint_053a_gap_collapse.json', 'w') as f:
    json.dump({'sprint': '053a', 'q': 3, 'g_c': g_c, 'sizes': [4, 6, 8],
               'data': results}, f, indent=2)
print("Results saved.", flush=True)

# === Data collapse optimization ===
print("\n=== Data collapse: optimize ν ===", flush=True)

def collapse_quality(nu, g_c, sizes, data):
    """Measure collapse quality by interpolation overlap."""
    # For each size, compute x = (g - g_c) * N^{1/ν}, y = Δ·N
    curves = {}
    for n in sizes:
        pts = data[f'n{n}']
        x = np.array([(p['g'] - g_c) * n**(1.0/nu) for p in pts])
        y = np.array([p['gap_x_n'] for p in pts])
        curves[n] = (x, y)

    # Compute overlap: for each pair of sizes, interpolate and compare
    total_err = 0.0
    n_pts = 0
    for i, n1 in enumerate(sizes):
        for n2 in sizes[i+1:]:
            x1, y1 = curves[n1]
            x2, y2 = curves[n2]
            # Find common x range
            x_min = max(x1.min(), x2.min())
            x_max = min(x1.max(), x2.max())
            if x_min >= x_max:
                continue
            x_common = np.linspace(x_min, x_max, 50)
            f1 = interp1d(x1, y1, kind='cubic', bounds_error=False, fill_value=np.nan)
            f2 = interp1d(x2, y2, kind='cubic', bounds_error=False, fill_value=np.nan)
            y1i = f1(x_common)
            y2i = f2(x_common)
            mask = np.isfinite(y1i) & np.isfinite(y2i)
            if mask.sum() < 5:
                continue
            # Normalized residual
            y_mean = (np.abs(y1i[mask]) + np.abs(y2i[mask])) / 2
            y_mean = np.maximum(y_mean, 1e-10)
            err = np.mean(((y1i[mask] - y2i[mask]) / y_mean)**2)
            total_err += err
            n_pts += 1

    return total_err / max(n_pts, 1)

sizes = [4, 6, 8]
nu_values = np.linspace(0.3, 2.5, 221)
quality = []
for nu in nu_values:
    q_val = collapse_quality(nu, g_c, sizes, results)
    quality.append(q_val)

quality = np.array(quality)
best_idx = np.argmin(quality)
best_nu = nu_values[best_idx]
best_q = quality[best_idx]

print(f"\nBest ν = {best_nu:.4f} (quality = {best_q:.6f})", flush=True)
print(f"Exact ν = {5/6:.4f}", flush=True)
print(f"Error: {abs(best_nu - 5/6) / (5/6) * 100:.1f}%", flush=True)

# Show top 5 candidates
sorted_idx = np.argsort(quality)[:5]
print("\nTop 5 ν values:")
for idx in sorted_idx:
    print(f"  ν = {nu_values[idx]:.4f}, quality = {quality[idx]:.6f}", flush=True)

# Also try with g_c as free parameter
print("\n=== Joint (ν, g_c) optimization ===", flush=True)
best_joint = (1.0, g_c, 1e10)
for nu in np.linspace(0.4, 2.0, 81):
    for gc_test in np.linspace(0.28, 0.40, 61):
        q_val = collapse_quality(nu, gc_test, sizes, results)
        if q_val < best_joint[2]:
            best_joint = (nu, gc_test, q_val)

print(f"Best joint: ν = {best_joint[0]:.4f}, g_c = {best_joint[1]:.4f} (quality = {best_joint[2]:.6f})", flush=True)
print(f"Exact: ν = {5/6:.4f}, g_c = {1/3:.4f}", flush=True)

# Compute Δ·N crossing points
print("\n=== Δ·N crossing points ===", flush=True)
for n1, n2 in [(4, 6), (4, 8), (6, 8)]:
    d1 = results[f'n{n1}']
    d2 = results[f'n{n2}']
    for j in range(len(g_values) - 1):
        diff_j = d1[j]['gap_x_n'] - d2[j]['gap_x_n']
        diff_j1 = d1[j+1]['gap_x_n'] - d2[j+1]['gap_x_n']
        if diff_j * diff_j1 < 0:
            g_cross = g_values[j] - diff_j * (g_values[j+1] - g_values[j]) / (diff_j1 - diff_j)
            print(f"  n={n1},{n2}: g_cross = {g_cross:.5f}", flush=True)

# Slope ratio for independent ν check
# dΔ·N/dg at g_c should scale as N^{1/ν}
print("\n=== Slope ratio method ===", flush=True)
slopes = {}
for n in sizes:
    pts = results[f'n{n}']
    g_arr = np.array([p['g'] for p in pts])
    y_arr = np.array([p['gap_x_n'] for p in pts])
    # Numerical derivative at g_c
    f = interp1d(g_arr, y_arr, kind='cubic')
    dg = 0.005
    slope = (f(g_c + dg) - f(g_c - dg)) / (2 * dg)
    slopes[n] = slope
    print(f"  n={n}: d(Δ·N)/dg at g_c = {slope:.4f}", flush=True)

# ν from slope ratios
for n1, n2 in [(4, 6), (4, 8), (6, 8)]:
    ratio = slopes[n2] / slopes[n1]
    nu_est = np.log(n2 / n1) / np.log(ratio)
    print(f"  n={n1},{n2}: slope ratio = {ratio:.4f}, ν = {nu_est:.4f}", flush=True)

# Save final results
with open('results/sprint_053a_gap_collapse.json', 'w') as f:
    json.dump({
        'sprint': '053a', 'q': 3, 'g_c_exact': g_c, 'sizes': sizes,
        'best_nu_fixed_gc': float(best_nu),
        'best_joint': {'nu': best_joint[0], 'g_c': best_joint[1], 'quality': best_joint[2]},
        'slopes': {str(k): float(v) for k, v in slopes.items()},
        'data': results
    }, f, indent=2)

print("\nDone!", flush=True)
