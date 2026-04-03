#!/usr/bin/env python3
"""Sprint 051b: Energy gap for q=4 Potts to find g_c.

q=4: dim = 4^n. n=4: 256, n=6: 4096, n=8: 65536.
Self-duality broken, g_c unknown. Sprint 050 pseudo-critical: g≈0.34 (n=8).
Scan g in [0.15, 0.55] to find Δ·N crossing.
"""
import numpy as np, json, time, warnings
warnings.filterwarnings('ignore')
from scipy.sparse.linalg import eigsh
from scipy.sparse import kron, csr_matrix

def potts_hamiltonian(q, n, g, J=1.0):
    """Build H = -J Σ δ(s_i,s_j) - g Σ (X_i + X_i†)."""
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

    # NN coupling: -J Σ_a P_a⊗P_a
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

    # Field: -g Σ (X + X†)
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
    vals = eigsh(H.real if q == 2 else H, k=min(4, H.shape[0]-1), which='SA', return_eigenvectors=False)
    vals = np.sort(vals.real)
    return vals[1] - vals[0], vals[0]

# Timing test
print("=== Timing: q=4, n=8 (dim=65536) ===", flush=True)
t0 = time.time()
gap, E0 = energy_gap(4, 8, 0.34)
dt = time.time() - t0
print(f"  gap={gap:.6f}, E0={E0:.6f}, t={dt:.1f}s", flush=True)

if dt > 15:
    print("  n=8 too slow, using n=4,6 only", flush=True)
    sizes = [4, 6]
else:
    sizes = [4, 6, 8]

results = {}
g_vals = np.linspace(0.15, 0.55, 41)

for n in sizes:
    dim = 4**n
    print(f"\n=== q=4 n={n} (dim={dim}) ===", flush=True)
    t0 = time.time()
    gaps = []
    for g in g_vals:
        if time.time() - t0 > 180:
            print(f"  Time limit at g={g:.3f}", flush=True)
            break
        gap, E0 = energy_gap(4, n, g)
        gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * n), 'E0': float(E0)})
    dt = time.time() - t0
    print(f"  Done: {len(gaps)} points in {dt:.1f}s", flush=True)
    results[f'n{n}'] = gaps

# Find crossings
print("\n=== q=4 Δ·N crossings ===", flush=True)
size_keys = sorted(results.keys(), key=lambda x: int(x[1:]))
crossings = []
for i, k1 in enumerate(size_keys):
    for k2 in size_keys[i+1:]:
        d1 = results[k1]
        d2 = results[k2]
        min_len = min(len(d1), len(d2))
        for j in range(min_len - 1):
            diff_j = d1[j]['gap_x_n'] - d2[j]['gap_x_n']
            diff_j1 = d1[j+1]['gap_x_n'] - d2[j+1]['gap_x_n']
            if diff_j * diff_j1 < 0:
                g_cross = d1[j]['g'] - diff_j * (d1[j+1]['g'] - d1[j]['g']) / (diff_j1 - diff_j)
                n1 = int(k1[1:]); n2 = int(k2[1:])
                print(f"  n={n1},{n2}: g_cross = {g_cross:.4f}", flush=True)
                crossings.append({'n1': n1, 'n2': n2, 'g_cross': float(g_cross)})

# Save
out = {'sprint': '051b', 'q': 4, 'method': 'exact_diag_gap',
       'crossings': crossings, 'data': results}
with open('results/sprint_051b_gap_q4.json', 'w') as f:
    json.dump(out, f, indent=2)

# Extrapolate g_c
if len(crossings) >= 2:
    # Use largest-size crossing as best estimate
    best = crossings[-1]
    print(f"\n*** Best estimate: g_c(q=4) ≈ {best['g_cross']:.4f} (n={best['n1']},{best['n2']})")
    # With ~2.5% finite-size correction (from q=3 validation):
    corrected = best['g_cross'] / (1 - 0.025)
    print(f"    Corrected (~+2.5%): g_c(q=4) ≈ {corrected:.4f}")

print("\nDone!", flush=True)
