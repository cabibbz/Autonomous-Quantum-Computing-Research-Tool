#!/usr/bin/env python3
"""Sprint 070c: Ground state overlap F(g,g+δg) across 2D transition.

F = |<ψ₀(g)|ψ₀(g+δg)>| is the ground state fidelity.
- At first-order: F drops sharply at g_c (level crossing)
- At continuous: F stays close to 1, smooth dip

χ_F = -2 ln F / (N δg²) includes ALL states implicitly.

Also: dE₀/dg from Hellmann-Feynman: dE/dg = <ψ₀|∂H/∂g|ψ₀> = <ψ₀|V|ψ₀>
This uses the wavefunction directly, avoiding numerical differentiation noise.

Systems: q=2 L=3,4; q=3 L=2,3; q=5 L=2,3
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
    """Build 2D hybrid Hamiltonian."""
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


def get_ground_state(L, q, g):
    """Return (E0, psi0) for L×L lattice."""
    n = L * L
    dim = q**n
    H = build_H_2d(L, L, q, g)
    if dim <= 5000:
        evals, evecs = np.linalg.eigh(H.toarray())
        return evals[0], evecs[:, 0]
    if HAS_GPU and dim > 100000:
        try:
            H_gpu = cp_csr(H)
            evals_gpu, evecs_gpu = cp_eigsh(H_gpu, k=2, which='SA')
            evals = cp.asnumpy(evals_gpu)
            evecs = cp.asnumpy(evecs_gpu)
            idx = np.argmin(evals)
            return float(evals[idx]), evecs[:, idx]
        except:
            pass
    evals, evecs = eigsh(H, k=2, which='SA')
    idx = np.argmin(evals)
    return float(evals[idx]), evecs[:, idx]


# ---- MAIN ----
results = {
    'experiment': '070c_2d_overlap',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
}

gc_2d = {2: 0.771, 3: 1.267, 5: 1.588}

# Systems
systems = [
    (2, 3, 50, 0.4),   # q, L, npts, half_width
    (2, 4, 30, 0.25),
    (3, 2, 50, 0.5),
    (3, 3, 30, 0.4),
    (5, 2, 50, 0.6),
    (5, 3, 15, 0.5),   # dim=2M, 15 pts
]

for q, L, npts, hw in systems:
    gc = gc_2d[q]
    n = L * L
    dim = q**n
    g_lo = max(0.05, gc - hw)
    g_hi = gc + hw
    g_values = np.linspace(g_lo, g_hi, npts)
    dg = g_values[1] - g_values[0]

    print(f"\nq={q}, L={L} (dim={dim}, {npts} pts, dg={dg:.4f})")

    # Get ground state at each g
    E0_arr = []
    psi_arr = []
    dEdg_HF = []  # Hellmann-Feynman dE/dg
    V = build_V(L, L, q)

    t0 = time.time()
    for gi, g in enumerate(g_values):
        t1 = time.time()
        E0, psi0 = get_ground_state(L, q, g)
        E0_arr.append(E0)
        psi_arr.append(psi0)
        # Hellmann-Feynman: dE/dg = <ψ₀|V|ψ₀>
        dEdg = float(psi0 @ (V @ psi0))
        dEdg_HF.append(dEdg / n)  # per site
        dt1 = time.time() - t1
        if (gi + 1) % 10 == 0 or dt1 > 10:
            print(f"  {gi+1}/{npts} g={g:.3f} E0/N={E0/n:.6f} ({dt1:.1f}s)")
    dt = time.time() - t0
    print(f"  Total: {dt:.1f}s")

    # Overlaps F(g, g+dg) = |<ψ₀(g)|ψ₀(g+dg)>|
    F_arr = []
    chiF_arr = []
    for i in range(len(psi_arr) - 1):
        overlap = abs(np.dot(psi_arr[i], psi_arr[i + 1]))
        F_arr.append(overlap)
        if overlap > 1e-15:
            chi = -2.0 * np.log(overlap) / (n * dg**2)
        else:
            chi = float('inf')
        chiF_arr.append(chi)
    g_mid = [(g_values[i] + g_values[i + 1]) / 2 for i in range(len(F_arr))]

    # d²E/dg² via Hellmann-Feynman (smoother than numerical diff)
    d2Edg2_HF = np.gradient(np.array(dEdg_HF), dg)

    # Peak fidelity susceptibility
    chiF_arr_np = np.array(chiF_arr)
    valid = np.isfinite(chiF_arr_np)
    if valid.any():
        peak_idx = np.argmax(chiF_arr_np[valid])
        # Map back to original indices
        valid_indices = np.where(valid)[0]
        peak_real_idx = valid_indices[peak_idx]
        peak_g = g_mid[peak_real_idx]
        peak_chiF = chiF_arr_np[peak_real_idx]
    else:
        peak_g, peak_chiF = None, None

    # Min overlap
    min_F = min(F_arr)
    min_F_g = g_mid[F_arr.index(min_F)]

    print(f"  Min overlap: F={min_F:.6f} at g={min_F_g:.3f}")
    print(f"  χ_F/N peak: {peak_chiF:.4f} at g={peak_g:.3f}" if peak_chiF else "  No valid χ_F")
    print(f"  dE/dg at g_c: {np.interp(gc, g_values, dEdg_HF):.6f}")

    # d²E/dg² peak
    d2_peak_idx = np.argmax(np.abs(d2Edg2_HF))
    d2_peak_g = g_values[d2_peak_idx]
    d2_peak_val = d2Edg2_HF[d2_peak_idx]
    print(f"  d²E/dg² peak: {d2_peak_val:.4f} at g={d2_peak_g:.3f}")

    key = f'q{q}_L{L}'
    results[key] = {
        'q': q, 'L': L, 'n': n, 'dim': dim,
        'g_values': [round(float(x), 6) for x in g_values],
        'E0_per_site': [round(float(x) / n, 8) for x in E0_arr],
        'dEdg_HF_per_site': [round(float(x), 8) for x in dEdg_HF],
        'd2Edg2_HF': [round(float(x), 6) for x in d2Edg2_HF],
        'g_mid': [round(float(x), 6) for x in g_mid],
        'overlap_F': [round(float(x), 8) for x in F_arr],
        'chiF_per_site': [round(float(x), 6) if np.isfinite(x) else None for x in chiF_arr],
        'min_F': round(float(min_F), 8),
        'min_F_g': round(float(min_F_g), 4),
        'peak_chiF_g': round(float(peak_g), 4) if peak_g else None,
        'peak_chiF': round(float(peak_chiF), 6) if peak_chiF else None,
        'd2_peak_g': round(float(d2_peak_g), 4),
        'd2_peak_val': round(float(d2_peak_val), 6),
        'time_s': round(dt, 1)
    }

    # Free wavefunctions to save memory
    del psi_arr

# === ANALYSIS ===
print("\n" + "=" * 60)
print("ANALYSIS: Overlap-derived fidelity susceptibility")
print("=" * 60)

# χ_F/N peak scaling
print("\nχ_F/N peak scaling:")
for q_val in [2, 3, 5]:
    Ls = []
    peaks = []
    for L in [2, 3, 4]:
        key = f'q{q_val}_L{L}'
        if key in results and results[key].get('peak_chiF') is not None:
            Ls.append(L)
            peaks.append(results[key]['peak_chiF'])
    if len(Ls) >= 2:
        log_L = np.log(np.array(Ls, dtype=float))
        log_p = np.log(np.array(peaks))
        if len(Ls) == 2:
            slope = (log_p[-1] - log_p[-2]) / (log_L[-1] - log_L[-2])
        else:
            slope = np.polyfit(log_L, log_p, 1)[0]
        print(f"  q={q_val}: {dict(zip(Ls, [f'{p:.4f}' for p in peaks]))}")
        print(f"    Scaling: L^{slope:.2f}")
        results[f'q{q_val}_chiF_slope'] = round(float(slope), 3)

# d²E/dg² peak scaling
print("\nd²E/dg² peak scaling:")
for q_val in [2, 3, 5]:
    Ls = []
    peaks = []
    for L in [2, 3, 4]:
        key = f'q{q_val}_L{L}'
        if key in results:
            Ls.append(L)
            peaks.append(abs(results[key]['d2_peak_val']))
    if len(Ls) >= 2:
        log_L = np.log(np.array(Ls, dtype=float))
        log_p = np.log(np.array(peaks))
        if len(Ls) == 2:
            slope = (log_p[-1] - log_p[-2]) / (log_L[-1] - log_L[-2])
        else:
            slope = np.polyfit(log_L, log_p, 1)[0]
        print(f"  q={q_val}: {dict(zip(Ls, [f'{p:.4f}' for p in peaks]))}")
        print(f"    Scaling: L^{slope:.2f}")
        results[f'q{q_val}_d2E_slope'] = round(float(slope), 3)

# Min overlap comparison
print("\nMinimum overlap (first-order → F→0 with L):")
for q_val in [2, 3, 5]:
    for L in [2, 3, 4]:
        key = f'q{q_val}_L{L}'
        if key in results:
            F = results[key]['min_F']
            g = results[key]['min_F_g']
            print(f"  q={q_val} L={L}: F_min={F:.6f} at g={g:.3f}")

# dE/dg at g_c (latent heat test)
print("\ndE₀/dg per site at g_c (latent heat = discontinuity):")
for q_val in [2, 3, 5]:
    gc = gc_2d[q_val]
    for L in [2, 3, 4]:
        key = f'q{q_val}_L{L}'
        if key in results:
            g_arr = results[key]['g_values']
            dE_arr = results[key]['dEdg_HF_per_site']
            val = np.interp(gc, g_arr, dE_arr)
            print(f"  q={q_val} L={L}: dE/dg(g_c)={val:.6f}")

# Save
with open("results/sprint_070c_overlap.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_070c_overlap.json")

from db_utils import record
for q_val in [2, 3, 5]:
    if f'q{q_val}_chiF_slope' in results:
        record(sprint=70, model='hybrid_2d', q=q_val, n=0,
               quantity='chiF_overlap_slope', value=results[f'q{q_val}_chiF_slope'],
               method='ground_state_overlap')
    if f'q{q_val}_d2E_slope' in results:
        record(sprint=70, model='hybrid_2d', q=q_val, n=0,
               quantity='d2E_HF_slope', value=results[f'q{q_val}_d2E_slope'],
               method='hellmann_feynman')
