#!/usr/bin/env python3
"""Sprint 069c: 2D entropy vs coupling — critical peak detection

Scan entanglement entropy S(g) for strip partition (w=1) across coupling g.
The entropy should peak at g_c and be suppressed in both ordered and disordered phases.

q=2 L=4 (validation, 16 sites, dim=65536): scan g=[0.2, 1.5], g_c=0.771
q=5 L=2 (fast, 4 sites, dim=625): scan g=[0.2, 4.0], g_c=1.588
q=5 L=3 (main, 9 sites, dim=2M): scan g=[0.5, 3.0], ~15 points

Also: compare S at g_c between L=2 and L=3 for q=5 to test if S grows with L
(evidence for area-law scaling at criticality).
"""
import numpy as np
import json
import time
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from scipy.sparse.linalg import eigsh

try:
    import cupy as cp
    from cupyx.scipy.sparse import csr_matrix as cp_csr
    from cupyx.scipy.sparse.linalg import eigsh as cp_eigsh
    HAS_GPU = True
    print("GPU available")
except ImportError:
    HAS_GPU = False
    print("No GPU, using CPU")

GPU_THRESHOLD = 50000


def build_H_2d(Lx, Ly, q, g):
    n = Lx * Ly
    dim = q**n
    potts_2site = np.zeros(q**2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q**2, q**2), format='csr')
    X = np.zeros((q, q))
    for s in range(q):
        X[(s + 1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)
    H = csr_matrix((dim, dim))
    bonds = []
    for y in range(Ly):
        for x in range(Lx):
            site = y * Lx + x
            nbr_h = y * Lx + (x + 1) % Lx
            if nbr_h != site:
                bonds.append((min(site, nbr_h), max(site, nbr_h)))
            nbr_v = ((y + 1) % Ly) * Lx + x
            if nbr_v != site:
                bonds.append((min(site, nbr_v), max(site, nbr_v)))
    bonds = list(set(bonds))
    for (i, j) in bonds:
        if j == i + 1:
            left = q**i
            right = q**(n - i - 2)
            op = potts_op
            if left > 1:
                op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1:
                op = sp_kron(op, sp_eye(right), format='csr')
        else:
            diag_vals = np.zeros(dim)
            for idx in range(dim):
                si = (idx // (q**i)) % q
                sj = (idx // (q**j)) % q
                diag_vals[idx] = 1.0 if si == sj else 0.0
            op = diags(diag_vals, 0, shape=(dim, dim), format='csr')
        H = H - op
    for i in range(n):
        left = q**i
        right = q**(n - i - 1)
        op = XpXd.copy()
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H


def get_ground_state(H, dim):
    if dim <= 2000:
        evals, evecs = np.linalg.eigh(H.toarray())
        return evals[0], evecs[:, 0]
    if HAS_GPU and dim > GPU_THRESHOLD:
        try:
            H_gpu = cp_csr(H)
            evals_gpu, evecs_gpu = cp_eigsh(H_gpu, k=1, which='SA')
            return float(cp.asnumpy(evals_gpu)[0]), cp.asnumpy(evecs_gpu[:, 0])
        except:
            pass
    evals, evecs = eigsh(H, k=1, which='SA')
    return evals[0], evecs[:, 0]


def entanglement_entropy(psi, q, n, subsystem_A):
    nA = len(subsystem_A)
    subsystem_B = [i for i in range(n) if i not in subsystem_A]
    dimA = q**nA
    dimB = q**(n - nA)
    psi_matrix = np.zeros((dimA, dimB))
    for idx in range(q**n):
        site_vals = []
        tmp = idx
        for i in range(n):
            site_vals.append(tmp % q)
            tmp //= q
        a_idx = sum(site_vals[site] * (q**ki) for ki, site in enumerate(subsystem_A))
        b_idx = sum(site_vals[site] * (q**ki) for ki, site in enumerate(subsystem_B))
        psi_matrix[a_idx, b_idx] += psi[idx]
    s = np.linalg.svd(psi_matrix, compute_uv=False)
    s = s[s > 1e-15]
    probs = s**2
    probs /= probs.sum()
    return float(-np.sum(probs * np.log(probs)))


def strip_partition(Lx, Ly, w):
    return [y * Lx + x for y in range(Ly) for x in range(w)]


results = {
    'experiment': '069c_2d_entropy_scan',
    'description': 'Entropy vs g scan in 2D — critical peak detection',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

# === q=2 L=4: validation scan ===
print("=" * 60)
print("q=2, L=4 (dim=65536): entropy scan g=[0.2, 1.5]")
print("=" * 60)

q, L = 2, 4
n = L*L
A_sites = strip_partition(L, L, 1)  # w=1 strip
g_vals = np.linspace(0.2, 1.5, 25)

q2_scan = []
for gi, g in enumerate(g_vals):
    t0 = time.time()
    H = build_H_2d(L, L, q, g)
    E0, psi = get_ground_state(H, q**n)
    S = entanglement_entropy(psi, q, n, A_sites)
    dt = time.time() - t0
    q2_scan.append({'g': round(float(g), 5), 'S': round(S, 8), 'E0_per_site': round(float(E0/n), 8)})
    if (gi+1) % 5 == 0:
        print(f"  {gi+1}/{len(g_vals)}: g={g:.3f}, S={S:.6f} ({dt:.1f}s)")

# Find peak
S_vals = [d['S'] for d in q2_scan]
peak_idx = np.argmax(S_vals)
print(f"\n  Peak: g={q2_scan[peak_idx]['g']:.3f}, S={q2_scan[peak_idx]['S']:.6f}")
print(f"  Expected g_c=0.771")
results['q2_L4_scan'] = q2_scan
results['q2_L4_peak'] = q2_scan[peak_idx]

# === q=5 L=2: fast scan ===
print(f"\n{'='*60}")
print("q=5, L=2 (dim=625): entropy scan g=[0.2, 4.0]")
print("=" * 60)

q, L = 5, 2
n = L*L
A_sites = strip_partition(L, L, 1)
g_vals = np.linspace(0.2, 4.0, 30)

q5_L2_scan = []
for gi, g in enumerate(g_vals):
    H = build_H_2d(L, L, q, g)
    E0, psi = get_ground_state(H, q**n)
    S = entanglement_entropy(psi, q, n, A_sites)
    q5_L2_scan.append({'g': round(float(g), 5), 'S': round(S, 8)})

S_vals = [d['S'] for d in q5_L2_scan]
peak_idx = np.argmax(S_vals)
print(f"  Peak: g={q5_L2_scan[peak_idx]['g']:.3f}, S={q5_L2_scan[peak_idx]['S']:.6f}")
print(f"  Expected g_c=1.588")
results['q5_L2_scan'] = q5_L2_scan
results['q5_L2_peak'] = q5_L2_scan[peak_idx]

# === q=5 L=3: main scan (dim=2M, ~15s per point) ===
print(f"\n{'='*60}")
print("q=5, L=3 (dim=1953125): entropy scan g=[0.5, 3.0]")
print("=" * 60)

q, L = 5, 3
n = L*L
A_sites = strip_partition(L, L, 1)
g_vals = np.linspace(0.5, 3.0, 15)

q5_L3_scan = []
t_total = time.time()
for gi, g in enumerate(g_vals):
    t0 = time.time()
    H = build_H_2d(L, L, q, g)
    E0, psi = get_ground_state(H, q**n)
    S = entanglement_entropy(psi, q, n, A_sites)
    dt = time.time() - t0
    q5_L3_scan.append({'g': round(float(g), 5), 'S': round(S, 8), 'E0_per_site': round(float(E0/n), 8)})
    print(f"  {gi+1}/{len(g_vals)}: g={g:.3f}, S={S:.6f} ({dt:.1f}s)")

dt_total = time.time() - t_total
S_vals = [d['S'] for d in q5_L3_scan]
peak_idx = np.argmax(S_vals)
print(f"\n  Peak: g={q5_L3_scan[peak_idx]['g']:.3f}, S={q5_L3_scan[peak_idx]['S']:.6f}")
print(f"  Expected g_c=1.588")
print(f"  Total: {dt_total:.0f}s")
results['q5_L3_scan'] = q5_L3_scan
results['q5_L3_peak'] = q5_L3_scan[peak_idx]

# === Analysis ===
print(f"\n{'='*60}")
print("ANALYSIS")
print("=" * 60)

# 1. Entropy peak location vs known g_c
print("\n1. Entropy peak location vs known g_c(2D):")
print(f"   q=2, L=4: peak g={results['q2_L4_peak']['g']:.3f}, known g_c=0.771")
print(f"   q=5, L=2: peak g={results['q5_L2_peak']['g']:.3f}, known g_c=1.588")
print(f"   q=5, L=3: peak g={results['q5_L3_peak']['g']:.3f}, known g_c=1.588")

# 2. Peak entropy growth: L=2 vs L=3 for q=5
S_peak_L2 = results['q5_L2_peak']['S']
S_peak_L3 = results['q5_L3_peak']['S']
print(f"\n2. Peak entropy growth with L (q=5):")
print(f"   L=2: S_peak={S_peak_L2:.6f} (boundary=4)")
print(f"   L=3: S_peak={S_peak_L3:.6f} (boundary=6)")
print(f"   S_peak(L=3)/S_peak(L=2) = {S_peak_L3/S_peak_L2:.3f}")
print(f"   boundary ratio = 6/4 = 1.500")
print(f"   Area law prediction: S ratio = boundary ratio = 1.5")

# 3. Sharpness: entropy at g_c vs g=2*g_c
print(f"\n3. Entropy profile sharpness (q=5 L=3):")
for d in q5_L3_scan:
    bar = '#' * int(d['S'] * 80)
    print(f"   g={d['g']:5.2f}: S={d['S']:.4f} {bar}")

# Save
with open("results/sprint_069c_2d_entropy_scan.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_069c_2d_entropy_scan.json")
