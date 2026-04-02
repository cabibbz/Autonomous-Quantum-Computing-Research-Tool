#!/usr/bin/env python3
"""Sprint 070a: 2D energy second derivative d²E₀/dg² near critical points.

Probe transition nature: continuous → d²E/dg² diverges as L^{α/ν}
                         first-order → d²E/dg² has delta-function peak ~ L^d

Systems: q=2 L=2,3,4; q=3 L=2,3; q=5 L=2,3
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

GPU_THRESHOLD = 100000


def build_H_2d(Lx, Ly, q, g):
    """Build 2D hybrid Hamiltonian on Lx x Ly periodic square lattice."""
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


def get_E0(Lx, Ly, q, g):
    """Get ground state energy."""
    H = build_H_2d(Lx, Ly, q, g)
    dim = H.shape[0]
    if dim < 500:
        evals = np.linalg.eigvalsh(H.toarray())
        return evals[0]
    if HAS_GPU and dim > GPU_THRESHOLD:
        try:
            H_gpu = cp_csr(H)
            evals_gpu, _ = cp_eigsh(H_gpu, k=2, which='SA')
            return float(cp.asnumpy(evals_gpu).min())
        except:
            pass
    evals, _ = eigsh(H, k=2, which='SA')
    return float(evals.min())


def energy_scan(L, q, g_values):
    """Compute E0 at each g."""
    n = L * L
    dim = q**n
    energies = []
    for gi, g in enumerate(g_values):
        t0 = time.time()
        E0 = get_E0(L, L, q, g)
        dt = time.time() - t0
        energies.append(E0)
        if (gi + 1) % 20 == 0:
            print(f"      {gi+1}/{len(g_values)} done ({dt:.1f}s/pt)")
    return np.array(energies)


def numerical_derivatives(g_values, energies, n_sites):
    """Compute dE/dg and d²E/dg² per site via finite differences."""
    dg = g_values[1] - g_values[0]
    # Central differences
    dEdg = np.gradient(energies / n_sites, dg)
    d2Edg2 = np.gradient(dEdg, dg)
    return dEdg, d2Edg2


# ---- MAIN ----
results = {
    'experiment': '070a_2d_energy_deriv',
    'description': '2D energy 2nd derivative near g_c',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
}

gc_2d = {2: 0.771, 3: 1.267, 5: 1.588}

# Systems to scan
systems = [
    (2, [2, 3, 4], 0.35),  # q, L_list, half-width of scan around g_c
    (3, [2, 3], 0.4),
    (5, [2, 3], 0.5),
]

for q, L_list, hw in systems:
    gc = gc_2d[q]
    print(f"\n{'='*60}")
    print(f"q={q}, g_c={gc}")
    print(f"{'='*60}")

    q_results = {}

    for L in L_list:
        n = L * L
        dim = q**n
        # More points for smaller systems
        if dim < 1000:
            npts = 80
        elif dim < 100000:
            npts = 60
        else:
            npts = 30  # L=3 q=5 is ~2M dim

        g_lo = max(0.05, gc - hw)
        g_hi = gc + hw
        g_values = np.linspace(g_lo, g_hi, npts)

        print(f"\n  L={L} (dim={dim}, {npts} pts, g=[{g_lo:.2f}, {g_hi:.2f}])")
        t0 = time.time()
        energies = energy_scan(L, q, g_values)
        dt = time.time() - t0
        print(f"    Total: {dt:.1f}s ({dt/npts:.1f}s/pt)")

        dEdg, d2Edg2 = numerical_derivatives(g_values, energies, n)

        # Find peak of |d²E/dg²|
        peak_idx = np.argmax(np.abs(d2Edg2))
        peak_g = g_values[peak_idx]
        peak_val = d2Edg2[peak_idx]

        print(f"    d²E/dg² peak: {peak_val:.4f} at g={peak_g:.3f}")
        print(f"    dE/dg at g_c: {np.interp(gc, g_values, dEdg):.6f}")

        q_results[L] = {
            'g_values': g_values.tolist(),
            'E0': energies.tolist(),
            'dEdg': dEdg.tolist(),
            'd2Edg2': d2Edg2.tolist(),
            'peak_g': round(float(peak_g), 4),
            'peak_d2E': round(float(peak_val), 6),
            'n_sites': n,
            'dim': dim,
            'time_s': round(dt, 1),
            'npts': npts
        }

    # Peak scaling
    print(f"\n  Peak |d²E/dg²| scaling (q={q}):")
    for L in L_list:
        p = q_results[L]['peak_d2E']
        pg = q_results[L]['peak_g']
        print(f"    L={L}: peak={p:.4f} at g={pg:.3f}")

    if len(L_list) >= 2:
        peaks = [abs(q_results[L]['peak_d2E']) for L in L_list]
        Ls = np.array(L_list, dtype=float)
        # Log-log slope
        if len(L_list) >= 2:
            log_L = np.log(Ls)
            log_p = np.log(np.abs(peaks))
            if len(L_list) == 2:
                slope = (log_p[1] - log_p[0]) / (log_L[1] - log_L[0])
            else:
                slope = np.polyfit(log_L, log_p, 1)[0]
            print(f"    Scaling exponent: {slope:.2f}")
            print(f"    (continuous q=2: expect ~0 (log), first-order: expect ~2)")
            q_results['scaling_exponent'] = round(float(slope), 3)

    results[f'q{q}'] = {str(L): v for L, v in q_results.items() if isinstance(L, int)}
    if 'scaling_exponent' in q_results:
        results[f'q{q}']['scaling_exponent'] = q_results['scaling_exponent']

# Save
with open("results/sprint_070a_energy_deriv.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_070a_energy_deriv.json")

# Record key quantities
from db_utils import record
for q in [2, 3, 5]:
    qr = results.get(f'q{q}', {})
    if 'scaling_exponent' in qr:
        record(sprint=70, model='hybrid_2d', q=q, n=0,
               quantity='d2E_scaling', value=qr['scaling_exponent'],
               method='energy_deriv_peak_vs_L')
