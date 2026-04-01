#!/usr/bin/env python3
"""Sprint 051c: Energy gap for q=5 Potts to find g_c.

q=5: dim = 5^n. n=4: 625, n=6: 15625, n=8: 390625 (might be too large).
Also try q=7 and q=10 at smaller sizes if feasible.
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

all_results = {}

# === q=5 ===
print("=== q=5 Timing ===", flush=True)
for n in [4, 6, 8]:
    dim = 5**n
    t0 = time.time()
    gap, E0 = energy_gap(5, n, 0.30)
    dt = time.time() - t0
    print(f"  n={n} (dim={dim}): gap={gap:.6f}, t={dt:.1f}s", flush=True)
    if dt > 20:
        print(f"  n={n} too slow for full scan, stopping here", flush=True)
        max_n_q5 = n - 2 if n > 4 else n
        break
else:
    max_n_q5 = 8

sizes_q5 = [n for n in [4, 6, 8] if n <= max_n_q5]
if 8 not in sizes_q5 and 5**8 <= 400000:
    # Try n=8 with fewer g points
    sizes_q5_sparse = [8]
else:
    sizes_q5_sparse = []

g_vals = np.linspace(0.15, 0.55, 41)
g_vals_sparse = np.linspace(0.20, 0.45, 26)  # Fewer points for large n

for n in sizes_q5:
    dim = 5**n
    print(f"\n=== q=5 n={n} (dim={dim}) ===", flush=True)
    t0 = time.time()
    gaps = []
    for g in g_vals:
        if time.time() - t0 > 120:
            print(f"  Time limit at g={g:.3f}", flush=True)
            break
        gap, E0 = energy_gap(5, n, g)
        gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * n), 'E0': float(E0)})
    dt = time.time() - t0
    print(f"  Done: {len(gaps)} points in {dt:.1f}s", flush=True)
    all_results[f'q5_n{n}'] = gaps

for n in sizes_q5_sparse:
    dim = 5**n
    print(f"\n=== q=5 n={n} SPARSE (dim={dim}) ===", flush=True)
    t0 = time.time()
    gaps = []
    for g in g_vals_sparse:
        if time.time() - t0 > 120:
            print(f"  Time limit at g={g:.3f}", flush=True)
            break
        gap, E0 = energy_gap(5, n, g)
        gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * n), 'E0': float(E0)})
    dt = time.time() - t0
    print(f"  Done: {len(gaps)} points in {dt:.1f}s", flush=True)
    all_results[f'q5_n{n}'] = gaps

# q=5 crossings
print("\n=== q=5 Δ·N crossings ===", flush=True)
q5_keys = sorted([k for k in all_results if k.startswith('q5_')], key=lambda x: int(x.split('n')[1]))
crossings_q5 = []
for i, k1 in enumerate(q5_keys):
    for k2 in q5_keys[i+1:]:
        d1 = all_results[k1]; d2 = all_results[k2]
        # Match g values
        g1_set = {round(d['g'], 4) for d in d1}
        g2_set = {round(d['g'], 4) for d in d2}
        common_g = sorted(g1_set & g2_set)
        if len(common_g) < 2: continue
        g1_map = {round(d['g'], 4): d['gap_x_n'] for d in d1}
        g2_map = {round(d['g'], 4): d['gap_x_n'] for d in d2}
        for j in range(len(common_g) - 1):
            ga, gb = common_g[j], common_g[j+1]
            diff_a = g1_map[ga] - g2_map[ga]
            diff_b = g1_map[gb] - g2_map[gb]
            if diff_a * diff_b < 0:
                g_cross = ga - diff_a * (gb - ga) / (diff_b - diff_a)
                n1 = int(k1.split('n')[1]); n2 = int(k2.split('n')[1])
                print(f"  n={n1},{n2}: g_cross = {g_cross:.4f}", flush=True)
                crossings_q5.append({'n1': n1, 'n2': n2, 'g_cross': float(g_cross)})

# === q=7 at n=4,6 ===
print("\n=== q=7 Timing ===", flush=True)
for n in [4, 6]:
    dim = 7**n
    t0 = time.time()
    gap, E0 = energy_gap(7, n, 0.20)
    dt = time.time() - t0
    print(f"  n={n} (dim={dim}): gap={gap:.6f}, t={dt:.1f}s", flush=True)

g_vals_q7 = np.linspace(0.10, 0.40, 31)
for n in [4, 6]:
    dim = 7**n
    if dim > 200000:
        print(f"  q=7 n={n} too large (dim={dim}), skipping", flush=True)
        continue
    print(f"\n=== q=7 n={n} (dim={dim}) ===", flush=True)
    t0 = time.time()
    gaps = []
    for g in g_vals_q7:
        if time.time() - t0 > 120:
            print(f"  Time limit at g={g:.3f}", flush=True)
            break
        gap, E0 = energy_gap(7, n, g)
        gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * n), 'E0': float(E0)})
    dt = time.time() - t0
    print(f"  Done: {len(gaps)} points in {dt:.1f}s", flush=True)
    all_results[f'q7_n{n}'] = gaps

# q=7 crossings
print("\n=== q=7 Δ·N crossings ===", flush=True)
q7_keys = sorted([k for k in all_results if k.startswith('q7_')], key=lambda x: int(x.split('n')[1]))
crossings_q7 = []
if len(q7_keys) >= 2:
    k1, k2 = q7_keys[0], q7_keys[1]
    d1 = all_results[k1]; d2 = all_results[k2]
    for j in range(min(len(d1), len(d2)) - 1):
        diff_j = d1[j]['gap_x_n'] - d2[j]['gap_x_n']
        diff_j1 = d1[j+1]['gap_x_n'] - d2[j+1]['gap_x_n']
        if diff_j * diff_j1 < 0:
            g_cross = d1[j]['g'] - diff_j * (d1[j+1]['g'] - d1[j]['g']) / (diff_j1 - diff_j)
            n1 = int(k1.split('n')[1]); n2 = int(k2.split('n')[1])
            print(f"  n={n1},{n2}: g_cross = {g_cross:.4f}", flush=True)
            crossings_q7.append({'n1': n1, 'n2': n2, 'g_cross': float(g_cross)})

# Save all
out = {'sprint': '051c', 'method': 'exact_diag_gap',
       'crossings_q5': crossings_q5, 'crossings_q7': crossings_q7,
       'data': all_results}
with open('results/sprint_051c_gap_q5_q7.json', 'w') as f:
    json.dump(out, f, indent=2)

print("\n=== Summary ===", flush=True)
if crossings_q5:
    best = crossings_q5[-1]
    corr = best['g_cross'] / (1 - 0.025)
    print(f"q=5: best crossing g={best['g_cross']:.4f} (n={best['n1']},{best['n2']}), corrected ≈ {corr:.4f}")
if crossings_q7:
    best = crossings_q7[-1]
    corr = best['g_cross'] / (1 - 0.025)
    print(f"q=7: best crossing g={best['g_cross']:.4f} (n={best['n1']},{best['n2']}), corrected ≈ {corr:.4f}")

print("\nDone!", flush=True)
