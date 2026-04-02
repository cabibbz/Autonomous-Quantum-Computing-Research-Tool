#!/usr/bin/env python3
"""Sprint 070b: 2D fidelity susceptibility & dE/dg discontinuity test.

Two probes:
1. Fidelity susceptibility χ_F = Σ_{n>0} |<n|V|0>|²/(E_n-E_0)² where V=∂H/∂g
   - Computed from eigenvectors (full diag for small dim, sparse for large)
2. dE₀/dg jump: latent heat test (first-order → step in dE/dg)

Focus on: q=2 L=3,4 (known continuous) and q=5 L=2,3 (test case)
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
    print("No GPU")


def build_H_2d(Lx, Ly, q, g):
    """Build 2D hybrid Hamiltonian on periodic square lattice."""
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
            if nbr_h != site: bonds.append((min(site, nbr_h), max(site, nbr_h)))
            nbr_v = ((y + 1) % Ly) * Lx + x
            if nbr_v != site: bonds.append((min(site, nbr_v), max(site, nbr_v)))
    bonds = list(set(bonds))
    for (i, j) in bonds:
        if j == i + 1:
            left = q**i; right = q**(n - i - 2)
            op = potts_op
            if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
            if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        else:
            diag_vals = np.zeros(dim)
            for idx in range(dim):
                si = (idx // (q**i)) % q
                sj = (idx // (q**j)) % q
                diag_vals[idx] = 1.0 if si == sj else 0.0
            op = diags(diag_vals, 0, shape=(dim, dim), format='csr')
        H = H - op
    for i in range(n):
        left = q**i; right = q**(n - i - 1)
        op = XpXd.copy()
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H


def build_V(Lx, Ly, q):
    """V = ∂H/∂g = -Σ_i (X_i + X_i†)"""
    n = Lx * Ly
    dim = q**n
    X = np.zeros((q, q))
    for s in range(q): X[(s + 1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)
    V = csr_matrix((dim, dim))
    for i in range(n):
        left = q**i; right = q**(n - i - 1)
        op = XpXd.copy()
        if left > 1: op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1: op = sp_kron(op, sp_eye(right), format='csr')
        V = V - op
    return V


def fidelity_and_energy(L, q, g, k_states=20):
    """Compute χ_F, E0, gap at given g."""
    n = L * L
    dim = q**n
    H = build_H_2d(L, L, q, g)

    if dim <= 5000:
        H_dense = H.toarray()
        evals, evecs = np.linalg.eigh(H_dense)
        V = build_V(L, L, q)
        V_dense = V.toarray()
        psi0 = evecs[:, 0]
        V_psi0 = V_dense @ psi0
        chi_F = 0.0
        for ni in range(1, min(dim, 200)):
            me = np.dot(evecs[:, ni], V_psi0)
            dE = evals[ni] - evals[0]
            if dE > 1e-12:
                chi_F += me**2 / dE**2
        return evals[0], evals[1] - evals[0], chi_F
    else:
        k = min(k_states, dim - 1)
        if HAS_GPU and dim > 100000:
            try:
                H_gpu = cp_csr(H)
                evals_gpu, evecs_gpu = cp_eigsh(H_gpu, k=k, which='SA')
                evals = cp.asnumpy(evals_gpu)
                evecs = cp.asnumpy(evecs_gpu)
            except:
                evals, evecs = eigsh(H, k=k, which='SA')
        else:
            evals, evecs = eigsh(H, k=k, which='SA')
        order = np.argsort(evals)
        evals = evals[order]
        evecs = evecs[:, order]
        V = build_V(L, L, q)
        psi0 = evecs[:, 0]
        V_psi0 = V @ psi0
        chi_F = 0.0
        for ni in range(1, k):
            me = np.dot(evecs[:, ni], V_psi0)
            dE = evals[ni] - evals[0]
            if dE > 1e-12:
                chi_F += me**2 / dE**2
        return evals[0], evals[1] - evals[0], chi_F


# ---- MAIN ----
results = {
    'experiment': '070b_2d_fidelity',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
}

gc_2d = {2: 0.771, 3: 1.267, 5: 1.588}

# Focus: q=2 L=3,4 and q=5 L=2,3 with dense scan near g_c
systems = [
    (2, 3, 40, 0.3, 10),    # q, L, npts, half_width, k_states
    (2, 4, 25, 0.2, 10),    # q=2 L=4 dim=65k — fewer points
    (5, 2, 40, 0.5, 10),    # q=5 L=2 dim=625
    (5, 3, 12, 0.4, 10),    # q=5 L=3 dim=2M — very few points
]

for q, L, npts, hw, kst in systems:
    gc = gc_2d[q]
    n = L * L
    dim = q**n
    g_lo = max(0.05, gc - hw)
    g_hi = gc + hw
    g_values = np.linspace(g_lo, g_hi, npts)

    print(f"\nq={q}, L={L} (dim={dim}, {npts} pts, g=[{g_lo:.2f},{g_hi:.2f}])")
    E0_arr, gap_arr, chiF_arr = [], [], []

    t0 = time.time()
    for gi, g in enumerate(g_values):
        t1 = time.time()
        E0, gap, chi_F = fidelity_and_energy(L, q, g, k_states=kst)
        E0_arr.append(E0)
        gap_arr.append(gap)
        chiF_arr.append(chi_F)
        dt1 = time.time() - t1
        if (gi + 1) % 5 == 0 or dt1 > 5:
            print(f"  {gi+1}/{npts} g={g:.3f} χ_F={chi_F:.2f} ({dt1:.1f}s)")
    dt = time.time() - t0
    print(f"  Total: {dt:.1f}s")

    E0_arr = np.array(E0_arr)
    chiF_arr = np.array(chiF_arr)
    dg = g_values[1] - g_values[0]
    dEdg = np.gradient(E0_arr / n, dg)

    peak_idx = np.argmax(chiF_arr / n)
    peak_g = g_values[peak_idx]
    peak_chiF = chiF_arr[peak_idx] / n

    print(f"  χ_F/N peak: {peak_chiF:.4f} at g={peak_g:.3f}")

    key = f'q{q}_L{L}'
    results[key] = {
        'q': q, 'L': L, 'n': n, 'dim': dim,
        'g_values': g_values.tolist(),
        'E0': E0_arr.tolist(),
        'gap': [float(x) for x in gap_arr],
        'chi_F': chiF_arr.tolist(),
        'chi_F_per_site': (chiF_arr / n).tolist(),
        'dEdg_per_site': dEdg.tolist(),
        'peak_g': round(float(peak_g), 4),
        'peak_chiF_per_site': round(float(peak_chiF), 6),
        'time_s': round(dt, 1)
    }

# === ANALYSIS ===
print("\n" + "=" * 60)
print("ANALYSIS")
print("=" * 60)

# q=2 L=3 vs L=4
print("\nq=2 χ_F/N scaling (L=3 vs L=4):")
p3 = results['q2_L3']['peak_chiF_per_site']
p4 = results['q2_L4']['peak_chiF_per_site']
slope_q2 = np.log(p4 / p3) / np.log(4 / 3)
print(f"  L=3: χ_F/N = {p3:.4f}")
print(f"  L=4: χ_F/N = {p4:.4f}")
print(f"  Scaling: L^{slope_q2:.2f}")
print(f"  (continuous ν=1: L^0; first-order: L^2)")
results['q2_scaling_3_4'] = round(float(slope_q2), 3)

# q=5 L=2 vs L=3
print("\nq=5 χ_F/N scaling (L=2 vs L=3):")
p2 = results['q5_L2']['peak_chiF_per_site']
p3 = results['q5_L3']['peak_chiF_per_site']
slope_q5 = np.log(p3 / p2) / np.log(3 / 2)
print(f"  L=2: χ_F/N = {p2:.4f}")
print(f"  L=3: χ_F/N = {p3:.4f}")
print(f"  Scaling: L^{slope_q5:.2f}")
print(f"  (continuous ν~0.85: L^{2/0.85-2:.2f}; first-order: L^2)")
results['q5_scaling_2_3'] = round(float(slope_q5), 3)

# dE/dg comparison at g_c
print("\ndE₀/dg at g_c (latent heat test):")
for q_val in [2, 5]:
    gc = gc_2d[q_val]
    for L in ([3, 4] if q_val == 2 else [2, 3]):
        key = f'q{q_val}_L{L}'
        g_arr = np.array(results[key]['g_values'])
        dEdg = np.array(results[key]['dEdg_per_site'])
        val = np.interp(gc, g_arr, dEdg)
        print(f"  q={q_val} L={L}: dE₀/dg(g_c) = {val:.6f}")

# Save
with open("results/sprint_070b_fidelity.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_070b_fidelity.json")

from db_utils import record
record(sprint=70, model='hybrid_2d', q=2, n=0,
       quantity='chiF_scaling_L3L4', value=results['q2_scaling_3_4'],
       method='fidelity_suscept')
record(sprint=70, model='hybrid_2d', q=5, n=0,
       quantity='chiF_scaling_L2L3', value=results['q5_scaling_2_3'],
       method='fidelity_suscept')
