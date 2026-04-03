#!/usr/bin/env python3
"""Sprint 053b: Energy gap at larger sizes for ν(q=3).

Push to n=10 (exact diag, dim=59049) and n=12,16 via DMRG excited states.
The n=4,6,8 slope ratio gave ν≈0.97 instead of exact 5/6.
Need larger sizes to see convergence toward 5/6.

Strategy:
1. Exact diag n=10 for q=3 (time it first)
2. DMRG ground + first excited state for n=12,16
3. Recompute slope ratios with larger sizes
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
    vals = eigsh(H, k=min(4, H.shape[0]-1), which='SA', return_eigenvectors=False)
    vals = np.sort(vals.real)
    return vals[1] - vals[0]

# === Timing: n=10 ===
print("=== Timing: q=3, n=10 (dim=59049) ===", flush=True)
t0 = time.time()
gap = energy_gap(3, 10, 0.333)
dt = time.time() - t0
print(f"  gap={gap:.6f}, t={dt:.1f}s", flush=True)

if dt > 20:
    print("  n=10 too slow for full scan, reducing g points", flush=True)
    g_values = np.linspace(0.25, 0.42, 18)
else:
    g_values = np.linspace(0.25, 0.42, 35)

# === Exact diag n=10 ===
print(f"\n=== q=3 n=10: {len(g_values)} g points ===", flush=True)
results_n10 = []
t0 = time.time()
for g in g_values:
    gap = energy_gap(3, 10, g)
    results_n10.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * 10)})
    if time.time() - t0 > 120:
        print(f"  Time limit at g={g:.3f}", flush=True)
        break
print(f"  Done: {len(results_n10)} points in {time.time()-t0:.1f}s", flush=True)

# Load previous n=4,6,8 data
try:
    with open('results/sprint_053a_gap_collapse.json') as f:
        prev = json.load(f)
    results = prev['data']
except:
    print("  Warning: could not load previous data", flush=True)
    results = {}

results['n10'] = results_n10

# Save
with open('results/sprint_053b_gap_larger.json', 'w') as f:
    json.dump({'sprint': '053b', 'q': 3, 'g_c': 1/3, 'data': results}, f, indent=2)
print("Saved.", flush=True)

# === DMRG for n=12, 16: energy gap via orthogonal excited state ===
print("\n=== DMRG excited state method ===", flush=True)
from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

class PottsSite(Site):
    def __init__(self, q):
        leg = npc.LegCharge.from_trivial(q)
        Site.__init__(self, leg, [str(a) for a in range(q)], sort_charge=False)
        for a in range(q):
            P = np.zeros((q, q), dtype=complex); P[a, a] = 1.0
            self.add_op(f'P{a}', P)
        X = np.zeros((q, q), dtype=complex)
        for a in range(q):
            X[(a + 1) % q, a] = 1.0
        self.add_op('X', X, hc='Xhc')
        self.add_op('Xhc', X.conj().T, hc='X')
        self.add_op('Xphc', X + X.conj().T, hc='Xphc')

class PottsChain(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, model_params):
        return PottsSite(model_params.get('q', 3))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 1.0)
        q = model_params.get('q', 3)
        for a in range(q):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'Xphc')

def dmrg_gap(n, q, g, chi_max=60):
    """Get E0 and E1 via DMRG with orthogonal_to."""
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})

    # Ground state
    np.random.seed(42 + n + int(g * 1000))
    init = [np.random.randint(q) for _ in range(n)]
    psi0 = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng0 = dmrg.TwoSiteDMRGEngine(psi0, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 40,
    })
    E0, _ = eng0.run()

    # First excited state orthogonal to ground state
    np.random.seed(123 + n + int(g * 1000))
    init1 = [(np.random.randint(q) + 1) % q for _ in range(n)]
    psi1 = MPS.from_product_state(model.lat.mps_sites(), init1, bc='finite')
    eng1 = dmrg.TwoSiteDMRGEngine(psi1, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 60,
        'orthogonal_to': [psi0],
    })
    E1, _ = eng1.run()

    return float(E0), float(E1), float(E1 - E0)

# Test DMRG gap at n=8 against exact diag
print("\n--- Validating DMRG gap vs exact diag at n=8 ---", flush=True)
t0 = time.time()
E0_dmrg, E1_dmrg, gap_dmrg = dmrg_gap(8, 3, 1/3, chi_max=60)
gap_exact = energy_gap(3, 8, 1/3)
print(f"  DMRG: E0={E0_dmrg:.8f}, E1={E1_dmrg:.8f}, gap={gap_dmrg:.6f}", flush=True)
print(f"  Exact: gap={gap_exact:.6f}", flush=True)
print(f"  Diff: {abs(gap_dmrg - gap_exact):.2e} ({abs(gap_dmrg - gap_exact)/gap_exact*100:.2f}%)", flush=True)
print(f"  Time: {time.time()-t0:.1f}s", flush=True)

# DMRG n=12 scan
g_dmrg = np.linspace(0.28, 0.40, 13)
for n in [12]:
    print(f"\n=== DMRG n={n}: {len(g_dmrg)} g points ===", flush=True)
    results_dmrg = []
    t_start = time.time()
    for g in g_dmrg:
        if time.time() - t_start > 150:
            print(f"  Time limit at g={g:.3f}", flush=True)
            break
        t0 = time.time()
        E0, E1, gap = dmrg_gap(n, 3, g, chi_max=60)
        dt = time.time() - t0
        print(f"  g={g:.3f}: gap={gap:.6f}, Δ·N={gap*n:.4f}, t={dt:.1f}s", flush=True)
        results_dmrg.append({'g': float(g), 'E0': E0, 'E1': E1, 'gap': gap, 'gap_x_n': float(gap * n)})

    results[f'n{n}'] = results_dmrg

    # Save incrementally
    with open('results/sprint_053b_gap_larger.json', 'w') as f:
        json.dump({'sprint': '053b', 'q': 3, 'g_c': 1/3, 'data': results}, f, indent=2)

# === Recompute slope ratios with all sizes ===
print("\n=== Slope ratios at g_c=1/3 ===", flush=True)
g_c = 1/3
slopes = {}
for nk in ['n4', 'n6', 'n8', 'n10', 'n12']:
    if nk not in results or len(results[nk]) < 5:
        continue
    n = int(nk[1:])
    pts = results[nk]
    g_arr = np.array([p['g'] for p in pts])
    y_arr = np.array([p['gap_x_n'] for p in pts])
    if g_c < g_arr.min() or g_c > g_arr.max():
        continue
    f = interp1d(g_arr, y_arr, kind='cubic')
    dg = 0.005
    if g_c - dg >= g_arr.min() and g_c + dg <= g_arr.max():
        slope = float((f(g_c + dg) - f(g_c - dg)) / (2 * dg))
        slopes[n] = slope
        print(f"  n={n}: d(Δ·N)/dg = {slope:.4f}", flush=True)

sizes_sorted = sorted(slopes.keys())
print("\n  Slope ratios → ν estimates:", flush=True)
for i, n1 in enumerate(sizes_sorted):
    for n2 in sizes_sorted[i+1:]:
        ratio = slopes[n2] / slopes[n1]
        if ratio > 0:
            nu_est = np.log(n2 / n1) / np.log(ratio)
            print(f"  n={n1},{n2}: ratio={ratio:.4f}, ν={nu_est:.4f}", flush=True)

# === Crossing points ===
print("\n=== Δ·N crossing points ===", flush=True)
for n1k, n2k in [('n4','n6'), ('n4','n8'), ('n6','n8'), ('n6','n10'), ('n8','n10'), ('n8','n12'), ('n10','n12')]:
    if n1k not in results or n2k not in results:
        continue
    d1, d2 = results[n1k], results[n2k]
    g1 = np.array([p['g'] for p in d1])
    g2 = np.array([p['g'] for p in d2])
    y1 = np.array([p['gap_x_n'] for p in d1])
    y2 = np.array([p['gap_x_n'] for p in d2])
    # Interpolate to common grid
    g_min = max(g1.min(), g2.min())
    g_max = min(g1.max(), g2.max())
    g_common = np.linspace(g_min, g_max, 200)
    f1 = interp1d(g1, y1, kind='cubic')
    f2 = interp1d(g2, y2, kind='cubic')
    diff = f1(g_common) - f2(g_common)
    for j in range(len(diff) - 1):
        if diff[j] * diff[j+1] < 0:
            gc = g_common[j] - diff[j] * (g_common[j+1] - g_common[j]) / (diff[j+1] - diff[j])
            n1, n2 = int(n1k[1:]), int(n2k[1:])
            print(f"  {n1k},{n2k}: g_cross = {gc:.5f}", flush=True)
            break

# Save final
with open('results/sprint_053b_gap_larger.json', 'w') as f:
    json.dump({
        'sprint': '053b', 'q': 3, 'g_c': 1/3,
        'slopes': {str(k): v for k, v in slopes.items()},
        'data': results
    }, f, indent=2)

print("\nDone!", flush=True)
