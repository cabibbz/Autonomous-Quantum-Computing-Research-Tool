#!/usr/bin/env python3
"""Sprint 052b: Energy gap crossing for q=10 to verify g_c prediction.

Prediction from scaling law: g_c(10) ≈ 0.58-0.62.
Strategy:
  1. Quick coarse scan at n=4 (0.1s/pt) over g=[0.3, 0.8]
  2. Quick coarse scan at n=6 (23s/pt) — only ~8-10 points affordable
  3. Find crossing of Δ·N curves for (n=4, n=6)
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

# === n=4 coarse scan (fast: 0.1s/pt) ===
print("=== q=10 n=4 (dim=10000): coarse scan g=[0.30, 0.80] ===", flush=True)
g_vals_n4 = np.linspace(0.30, 0.80, 26)
data_n4 = []
t0 = time.time()
for g in g_vals_n4:
    gap, E0 = energy_gap(10, 4, g)
    data_n4.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * 4), 'E0': float(E0)})
    print(f"  g={g:.3f}: Δ·4={gap*4:.4f}", flush=True)
print(f"  n=4 done in {time.time()-t0:.1f}s\n", flush=True)
results['q10_n4'] = data_n4

# Identify approximate crossing region from n=4 curve shape
# The crossing with n=6 should be where Δ·N changes curvature
# First, run n=6 at a few coarse points to bracket the crossing
print("=== q=10 n=6 (dim=1000000): targeted scan ===", flush=True)
# Start with wide spacing, then refine
g_coarse = [0.35, 0.45, 0.55, 0.60, 0.65, 0.70, 0.75]
data_n6 = []
t0 = time.time()
for g in g_coarse:
    elapsed = time.time() - t0
    if elapsed > 200:
        print(f"  Time limit at g={g:.2f} ({elapsed:.0f}s)", flush=True)
        break
    gap, E0 = energy_gap(10, 6, g)
    data_n6.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * 6), 'E0': float(E0)})
    print(f"  g={g:.3f}: Δ·6={gap*6:.4f} ({time.time()-t0:.1f}s)", flush=True)
print(f"  n=6 done in {time.time()-t0:.1f}s\n", flush=True)
results['q10_n6'] = data_n6

# Find crossings between n=4 and n=6
print("=== Finding n=4,6 crossings ===", flush=True)
# Interpolate n=4 data at n=6 g values
from scipy.interpolate import interp1d

g_n4_arr = np.array([d['g'] for d in data_n4])
gxn_n4_arr = np.array([d['gap_x_n'] for d in data_n4])
f_n4 = interp1d(g_n4_arr, gxn_n4_arr, kind='cubic', fill_value='extrapolate')

crossings = []
for i in range(len(data_n6) - 1):
    g1, g2 = data_n6[i]['g'], data_n6[i+1]['g']
    gxn6_1, gxn6_2 = data_n6[i]['gap_x_n'], data_n6[i+1]['gap_x_n']
    gxn4_1, gxn4_2 = float(f_n4(g1)), float(f_n4(g2))
    diff1 = gxn4_1 - gxn6_1
    diff2 = gxn4_2 - gxn6_2
    print(f"  g=[{g1:.2f},{g2:.2f}]: Δ·4-Δ·6 = [{diff1:.4f}, {diff2:.4f}]", flush=True)
    if diff1 * diff2 < 0:
        g_cross = g1 - diff1 * (g2 - g1) / (diff2 - diff1)
        print(f"  *** CROSSING at g = {g_cross:.4f} ***", flush=True)
        crossings.append({'n1': 4, 'n2': 6, 'g_cross': float(g_cross)})

if not crossings:
    print("  No crossings found in scanned range!", flush=True)
    # Report the ordering
    for d6 in data_n6:
        g = d6['g']
        gxn4 = float(f_n4(g))
        gxn6 = d6['gap_x_n']
        sign = "n4 > n6" if gxn4 > gxn6 else "n6 > n4"
        print(f"  g={g:.2f}: Δ·4={gxn4:.4f}, Δ·6={gxn6:.4f} ({sign})", flush=True)

# Summary
print(f"\n=== Summary ===")
if crossings:
    raw_gc = crossings[0]['g_cross']
    corrected_gc = raw_gc * 1.025  # ~2.5% FSS correction from validated q=2,3
    print(f"Raw crossing: g_c = {raw_gc:.4f}")
    print(f"Corrected (×1.025): g_c ≈ {corrected_gc:.4f}")
    print(f"Prediction (power_offset): 0.6164")
    print(f"Prediction (logarithmic): 0.5808")
    print(f"Error vs power_offset: {abs(corrected_gc - 0.6164)/0.6164*100:.1f}%")
    print(f"Error vs logarithmic: {abs(corrected_gc - 0.5808)/0.5808*100:.1f}%")
else:
    print("No crossing found — prediction may be outside scanned range")

results['crossings'] = crossings
results['prediction'] = {'power_offset': 0.6164, 'logarithmic': 0.5808}

with open('results/sprint_052b_gap_q10.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved to results/sprint_052b_gap_q10.json")
