#!/usr/bin/env python3
"""Sprint 051e: q=7 energy gap extended to g=0.60, and q=10 n=4,6.

q=7 showed no crossing in [0.05, 0.35]. If g_c increases with q
(q=2: 0.25, q=3: 0.33, q=4: 0.38, q=5: 0.43), then q=7 crossing
might be near 0.50.
"""
import numpy as np, json, time, warnings
warnings.filterwarnings('ignore')
from scipy.sparse.linalg import eigsh
from scipy.sparse import kron, csr_matrix

def potts_hamiltonian(q, n, g, J=1.0):
    I = csr_matrix(np.eye(q))
    X = np.zeros((q, q), dtype=complex)
    for a in range(q):
        X[(a+1) % q, a] = 1.0
    Xphc = csr_matrix(X + X.conj().T)
    projectors = []
    for a in range(q):
        P = np.zeros((q, q)); P[a, a] = 1.0
        projectors.append(csr_matrix(P))
    dim = q**n
    H = csr_matrix((dim, dim), dtype=complex)
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
    for i in range(n):
        op = csr_matrix(np.eye(1))
        for j in range(n):
            if j == i:
                op = kron(op, Xphc)
            else:
                op = kron(op, I)
        H += -g * op
    return H

def energy_gap(q, n, g, J=1.0):
    H = potts_hamiltonian(q, n, g, J)
    vals = eigsh(H, k=min(4, H.shape[0]-1), which='SA', return_eigenvectors=False)
    vals = np.sort(vals.real)
    return vals[1] - vals[0], vals[0]

results = {}

# === q=7 extended to g=0.70 ===
print("=== q=7 extended range [0.30, 0.70] ===", flush=True)
g_vals = np.linspace(0.30, 0.70, 41)
for n in [4, 6]:
    dim = 7**n
    print(f"  n={n} (dim={dim})...", flush=True)
    t0 = time.time()
    gaps = []
    for g in g_vals:
        if time.time() - t0 > 90:
            print(f"    Time limit at g={g:.3f}", flush=True)
            break
        gap, E0 = energy_gap(7, n, g)
        gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * n)})
    dt = time.time() - t0
    print(f"  Done: {len(gaps)} points in {dt:.1f}s", flush=True)
    results[f'q7_n{n}'] = gaps

# q=7 crossings
print("\n  q=7 crossings [0.30, 0.70]:", flush=True)
d4 = results['q7_n4']
d6 = results['q7_n6']
crossings = []
min_len = min(len(d4), len(d6))
for j in range(min_len - 1):
    diff_j = d4[j]['gap_x_n'] - d6[j]['gap_x_n']
    diff_j1 = d4[j+1]['gap_x_n'] - d6[j+1]['gap_x_n']
    if diff_j * diff_j1 < 0:
        g_cross = d4[j]['g'] - diff_j * (d4[j+1]['g'] - d4[j]['g']) / (diff_j1 - diff_j)
        print(f"    n=4,6: g_cross = {g_cross:.4f}", flush=True)
        crossings.append({'n1': 4, 'n2': 6, 'g_cross': float(g_cross)})

if not crossings:
    print("    Still no crossings!", flush=True)
    # Print some values to understand the trend
    print("    g     Δ·4    Δ·6    diff", flush=True)
    for i in range(0, min_len, 5):
        diff = d4[i]['gap_x_n'] - d6[i]['gap_x_n']
        print(f"    {d4[i]['g']:.2f}  {d4[i]['gap_x_n']:.4f}  {d6[i]['gap_x_n']:.4f}  {diff:+.4f}", flush=True)

# === q=10 n=4 (dim=10000) feasibility ===
print("\n=== q=10 n=4 (dim=10000) ===", flush=True)
t0 = time.time()
gap, E0 = energy_gap(10, 4, 0.30)
dt = time.time() - t0
print(f"  Timing: gap={gap:.6f}, t={dt:.1f}s", flush=True)

if dt < 5:
    g_vals_q10 = np.linspace(0.20, 0.80, 31)
    print(f"  Scanning {len(g_vals_q10)} points...", flush=True)
    gaps = []
    for g in g_vals_q10:
        gap, E0 = energy_gap(10, 4, g)
        gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * 4)})
    results['q10_n4'] = gaps
    print("  Done.", flush=True)

    # Try n=6 (dim=10^6 = 1M - probably too large for exact diag)
    print("\n=== q=10 n=6 (dim=1000000) — might fail ===", flush=True)
    try:
        t0 = time.time()
        gap, E0 = energy_gap(10, 6, 0.50)
        dt = time.time() - t0
        print(f"  Timing: gap={gap:.6f}, t={dt:.1f}s", flush=True)
        if dt < 20:
            gaps = []
            for g in g_vals_q10:
                if time.time() - t0 > 120:
                    print(f"    Time limit at g={g:.3f}", flush=True)
                    break
                gap, E0 = energy_gap(10, 6, g)
                gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * 6)})
            results['q10_n6'] = gaps
    except Exception as e:
        print(f"  Failed: {e}", flush=True)

# q=10 crossings
if 'q10_n4' in results and 'q10_n6' in results:
    print("\n  q=10 crossings:", flush=True)
    d4 = results['q10_n4']
    d6 = results['q10_n6']
    min_len = min(len(d4), len(d6))
    for j in range(min_len - 1):
        diff_j = d4[j]['gap_x_n'] - d6[j]['gap_x_n']
        diff_j1 = d4[j+1]['gap_x_n'] - d6[j+1]['gap_x_n']
        if diff_j * diff_j1 < 0:
            g_cross = d4[j]['g'] - diff_j * (d4[j+1]['g'] - d4[j]['g']) / (diff_j1 - diff_j)
            print(f"    n=4,6: g_cross = {g_cross:.4f}", flush=True)

# Summary: g_c trend
print("\n=== g_c(q) trend from energy gap ===", flush=True)
print("  q=2: g_c = 0.250 (exact, self-dual)", flush=True)
print("  q=3: g_c = 0.333 (exact, self-dual)", flush=True)
print("  q=4: g_c ≈ 0.382-0.392 (n=6,8 crossing + correction)", flush=True)
print("  q=5: g_c ≈ 0.430-0.441 (n=6,8 crossing + correction)", flush=True)

# Save
out = {'sprint': '051e', 'method': 'exact_diag_gap_extended',
       'crossings_q7': crossings, 'data': results}
with open('results/sprint_051e_gap_extended.json', 'w') as f:
    json.dump(out, f, indent=2)

print("\nDone!", flush=True)
