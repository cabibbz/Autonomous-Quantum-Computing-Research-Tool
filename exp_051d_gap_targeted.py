#!/usr/bin/env python3
"""Sprint 051d: Targeted gap scans.

1. q=5 n=6,8 near g=0.42 (narrow range, fewer points)
2. q=7 n=4,6: check if crossing is below 0.10 or if gap is just very small
3. q=10 n=4: feasibility check (10^4 = 10000)
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

# === q=5 n=8 targeted near crossing ===
print("=== q=5 n=8 targeted (g=0.35-0.50) ===", flush=True)
g_vals_target = np.linspace(0.35, 0.50, 16)
for n in [6, 8]:
    dim = 5**n
    print(f"  n={n} (dim={dim})...", flush=True)
    t0 = time.time()
    gaps = []
    for g in g_vals_target:
        if time.time() - t0 > 120:
            print(f"    Time limit at g={g:.3f}", flush=True)
            break
        gap, E0 = energy_gap(5, n, g)
        gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * n)})
        print(f"    g={g:.3f}: Δ·n={gap*n:.4f}", flush=True)
    results[f'q5_n{n}'] = gaps

# q=5 n=6,8 crossing
print("\n  q=5 n=6,8 crossings:", flush=True)
d6 = results.get('q5_n6', [])
d8 = results.get('q5_n8', [])
crossings = []
if d6 and d8:
    min_len = min(len(d6), len(d8))
    for j in range(min_len - 1):
        diff_j = d6[j]['gap_x_n'] - d8[j]['gap_x_n']
        diff_j1 = d6[j+1]['gap_x_n'] - d8[j+1]['gap_x_n']
        if diff_j * diff_j1 < 0:
            g_cross = d6[j]['g'] - diff_j * (d6[j+1]['g'] - d6[j]['g']) / (diff_j1 - diff_j)
            print(f"    n=6,8: g_cross = {g_cross:.4f}", flush=True)
            crossings.append({'n1': 6, 'n2': 8, 'g_cross': float(g_cross)})

# === q=7 extended range ===
print("\n=== q=7 gap check: very small g range ===", flush=True)
# Check g near 0 to see ordering
for g in [0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]:
    gap4, _ = energy_gap(7, 4, g)
    gap6, _ = energy_gap(7, 6, g)
    sign = "n4 > n6" if gap4*4 > gap6*6 else "n6 > n4"
    print(f"  g={g:.2f}: Δ·4={gap4*4:.4f}, Δ·6={gap6*6:.4f} ({sign})", flush=True)

# Dense scan around expected crossing
print("\n=== q=7 dense scan ===", flush=True)
g_vals_q7 = np.linspace(0.05, 0.35, 61)
for n in [4, 6]:
    dim = 7**n
    print(f"  n={n} (dim={dim})...", flush=True)
    t0 = time.time()
    gaps = []
    for g in g_vals_q7:
        if time.time() - t0 > 90:
            print(f"    Time limit at g={g:.3f}", flush=True)
            break
        gap, E0 = energy_gap(7, n, g)
        gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * n)})
    print(f"  Done: {len(gaps)} points in {time.time()-t0:.1f}s", flush=True)
    results[f'q7_n{n}'] = gaps

# q=7 crossings
print("\n  q=7 crossings:", flush=True)
d4 = results.get('q7_n4', [])
d6 = results.get('q7_n6', [])
crossings_q7 = []
if d4 and d6:
    min_len = min(len(d4), len(d6))
    for j in range(min_len - 1):
        diff_j = d4[j]['gap_x_n'] - d6[j]['gap_x_n']
        diff_j1 = d4[j+1]['gap_x_n'] - d6[j+1]['gap_x_n']
        if diff_j * diff_j1 < 0:
            g_cross = d4[j]['g'] - diff_j * (d4[j+1]['g'] - d4[j]['g']) / (diff_j1 - diff_j)
            print(f"    n=4,6: g_cross = {g_cross:.4f}", flush=True)
            crossings_q7.append({'n1': 4, 'n2': 6, 'g_cross': float(g_cross)})

if not crossings_q7:
    print("    No crossings found!", flush=True)
    # Check if curves are always ordered one way
    if d4 and d6:
        print(f"    At g=0.05: n4={d4[0]['gap_x_n']:.4f}, n6={d6[0]['gap_x_n']:.4f}", flush=True)
        print(f"    At g=0.35: n4={d4[-1]['gap_x_n']:.4f}, n6={d6[-1]['gap_x_n']:.4f}", flush=True)

# Save
out = {'sprint': '051d', 'method': 'exact_diag_gap_targeted',
       'crossings_q5': crossings, 'crossings_q7': crossings_q7,
       'data': results}
with open('results/sprint_051d_gap_targeted.json', 'w') as f:
    json.dump(out, f, indent=2)

print("\nDone!", flush=True)
