#!/usr/bin/env python3
"""Sprint 053e: ν(q=4) and ν(q=7) at true g_c values.

q=4: g_c ≈ 0.392 (Sprint 052). Marginal point in 2D classical.
q=7: g_c ≈ 0.535 (Sprint 052).

Is ν ≈ 5/6 universal for all q≥3?

Method: exact diag energy gap slopes at n=4,6.
Apply corrected power-law with b=0.86 calibration from q=3.

Also do q=2 as sanity check (should give ν=1).
"""
import numpy as np, json, time, warnings
warnings.filterwarnings('ignore')
from scipy.sparse.linalg import eigsh
from scipy.sparse import kron, csr_matrix
from scipy.interpolate import interp1d

def potts_hamiltonian(q, n, g, J=1.0):
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
    vals = eigsh(H, k=min(6, H.shape[0]-1), which='SA', return_eigenvectors=False)
    vals = np.sort(vals.real)
    return vals[1] - vals[0]

b_cal = 0.86  # From q=3 calibration

def nu_corrected(s_n1, s_n2, n1, n2, b=0.86):
    """Extract ν with finite-size correction."""
    corr_ratio = (1 + b/n2) / (1 + b/n1)
    eff_ratio = (s_n2 / s_n1) / corr_ratio
    if eff_ratio <= 1:
        return None
    return 1.0 / (np.log(eff_ratio) / np.log(n2/n1))

results_all = {}

for q, gc_est, g_range in [
    (2, 0.250, (0.12, 0.40)),
    (4, 0.392, (0.25, 0.52)),
    (7, 0.535, (0.38, 0.70)),
]:
    print(f"\n{'='*60}", flush=True)
    print(f"=== q={q}, g_c ≈ {gc_est} ===", flush=True)

    g_values = np.linspace(g_range[0], g_range[1], 57)
    data = {}
    for n in [4, 6, 8]:
        dim = q**n
        if dim > 100000:
            print(f"  n={n} (dim={dim}) — skipping (too large)", flush=True)
            continue
        t0 = time.time()
        gaps = []
        for g in g_values:
            gap = energy_gap(q, n, g)
            gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * n)})
        dt = time.time() - t0
        print(f"  n={n} (dim={dim}): {len(g_values)} points in {dt:.1f}s", flush=True)
        data[f'n{n}'] = gaps

    # Find crossings
    print(f"\n  Δ·N crossings:", flush=True)
    size_keys = sorted(data.keys(), key=lambda x: int(x[1:]))
    for i, nk1 in enumerate(size_keys):
        for nk2 in size_keys[i+1:]:
            n1, n2 = int(nk1[1:]), int(nk2[1:])
            g1 = np.array([p['g'] for p in data[nk1]])
            y1 = np.array([p['gap_x_n'] for p in data[nk1]])
            g2 = np.array([p['g'] for p in data[nk2]])
            y2 = np.array([p['gap_x_n'] for p in data[nk2]])
            diff = y1 - y2
            for j in range(len(diff) - 1):
                if diff[j] * diff[j+1] < 0:
                    gc = g_values[j] - diff[j] * (g_values[j+1] - g_values[j]) / (diff[j+1] - diff[j])
                    print(f"    ({n1},{n2}): g_cross = {gc:.5f}", flush=True)
                    break

    # Slope analysis at g_c
    print(f"\n  Slope analysis at g_c = {gc_est}:", flush=True)
    slopes = {}
    for nk in size_keys:
        n = int(nk[1:])
        pts = data[nk]
        g_arr = np.array([p['g'] for p in pts])
        y_arr = np.array([p['gap_x_n'] for p in pts])
        if gc_est < g_arr.min() + 0.01 or gc_est > g_arr.max() - 0.01:
            continue
        f = interp1d(g_arr, y_arr, kind='cubic')
        dg = 0.005
        slope = float((f(gc_est + dg) - f(gc_est - dg)) / (2 * dg))
        slopes[n] = slope
        print(f"    n={n}: slope = {slope:.4f}", flush=True)

    # Raw and corrected ν
    sizes_sorted = sorted(slopes.keys())
    print(f"\n  ν estimates:", flush=True)
    for i, n1 in enumerate(sizes_sorted):
        for n2 in sizes_sorted[i+1:]:
            ratio = slopes[n2] / slopes[n1]
            nu_raw = np.log(n2/n1) / np.log(ratio)
            nu_c = nu_corrected(slopes[n1], slopes[n2], n1, n2)
            print(f"    ({n1},{n2}): raw ν={nu_raw:.4f}, corrected ν={nu_c:.4f}" if nu_c else f"    ({n1},{n2}): raw ν={nu_raw:.4f}, corrected FAILED", flush=True)

    results_all[f'q{q}'] = {
        'g_c': gc_est,
        'slopes': {str(k): v for k, v in slopes.items()},
        'data': {k: v for k, v in data.items()}
    }

    # Save incrementally
    with open('results/sprint_053e_nu_multi_q.json', 'w') as f:
        json.dump({'sprint': '053e', 'b_cal': b_cal, 'results': results_all}, f, indent=2)

# === SUMMARY TABLE ===
print(f"\n{'='*60}", flush=True)
print("=== SUMMARY: ν(q) from corrected energy gap slopes ===", flush=True)
print(f"{'q':>3} | {'g_c':>6} | {'sizes':>8} | {'raw ν':>7} | {'corr ν':>7} | {'exact ν':>7}", flush=True)
print("-" * 55, flush=True)

# q=2 from this run
q2_slopes = results_all.get('q2', {}).get('slopes', {})
q4_slopes = results_all.get('q4', {}).get('slopes', {})
q7_slopes = results_all.get('q7', {}).get('slopes', {})

for q_val, gc, s_dict, exact in [
    (2, 0.250, q2_slopes, 1.0),
    (3, 0.333, None, 5/6),
    (4, 0.392, q4_slopes, None),
    (5, 0.441, None, None),
    (7, 0.535, q7_slopes, None),
]:
    if s_dict and len(s_dict) >= 2:
        sk = sorted(s_dict.keys(), key=lambda x: int(x))
        n1, n2 = int(sk[0]), int(sk[-1])
        ratio = float(s_dict[sk[-1]]) / float(s_dict[sk[0]])
        nu_r = np.log(n2/n1) / np.log(ratio)
        nu_c = nu_corrected(float(s_dict[sk[0]]), float(s_dict[sk[-1]]), n1, n2)
        exact_str = f"{exact:.4f}" if exact else "?"
        print(f"  {q_val} | {gc:.3f} | ({n1},{n2}) | {nu_r:.4f} | {nu_c:.4f if nu_c else 'N/A':>7} | {exact_str}", flush=True)
    elif q_val == 3:
        print(f"  {q_val} | {gc:.3f} | (4,10)  | 0.9603 | 0.8598 | {5/6:.4f}", flush=True)
    elif q_val == 5:
        print(f"  {q_val} | {gc:.3f} | (4,6)   | 0.9742 | 0.8501 | ?", flush=True)

print("\nDone!", flush=True)
