#!/usr/bin/env python3
"""Sprint 054d: Push c(q=3) to larger n with higher chi.

Previous: c pairwise 0.934 → 0.907 → 0.884 at chi=80, n=8-24.
Need n=32,48,64 to see convergence toward c=4/5=0.800.
Also: exact diag cross-check at n=8 to verify DMRG.
"""
import numpy as np, json, time, warnings
warnings.filterwarnings('ignore')
from scipy.sparse.linalg import eigsh
from scipy.sparse import kron, csr_matrix

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


def run_dmrg(q, n, g, chi_max=120):
    t0 = time.time()
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n + int(g * 1000))
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 50,
    })
    E0, _ = eng.run()
    S_half = float(psi.entanglement_entropy()[n // 2 - 1])
    chi_actual = max(psi.chi)
    return float(E0), S_half, chi_actual, time.time() - t0


# === Exact diag cross-check at n=8 ===
print("=== Exact diag cross-check: q=3, n=8, g=1/3 ===", flush=True)
q, n, g = 3, 8, 1.0/3
dim = q**n  # 6561
I_q = csr_matrix(np.eye(q))
X = np.zeros((q, q), dtype=complex)
for a in range(q):
    X[(a+1) % q, a] = 1.0
Xphc = csr_matrix(X + X.conj().T)
projs = [csr_matrix(np.diag([1.0 if b == a else 0.0 for b in range(q)])) for a in range(q)]
H = csr_matrix((dim, dim), dtype=complex)
for i in range(n-1):
    for a in range(q):
        op = csr_matrix(np.eye(1))
        for j in range(n):
            if j == i or j == i+1:
                op = kron(op, projs[a])
            else:
                op = kron(op, I_q)
        H += -op
for i in range(n):
    op = csr_matrix(np.eye(1))
    for j in range(n):
        if j == i:
            op = kron(op, Xphc)
        else:
            op = kron(op, I_q)
    H += -g * op

vals, vecs = eigsh(H, k=1, which='SA')
psi_gs = vecs[:, 0]
# Compute half-chain entropy
rho = np.outer(psi_gs, psi_gs.conj())
d_A = q**(n//2)
d_B = q**(n - n//2)
rho_A = rho.reshape(d_A, d_B, d_A, d_B).trace(axis1=1, axis2=3)
eigvals_rho = np.linalg.eigvalsh(rho_A)
eigvals_rho = eigvals_rho[eigvals_rho > 1e-14]
S_exact = -np.sum(eigvals_rho * np.log(eigvals_rho))
print(f"  Exact: E={vals[0]:.8f}, S_half={S_exact:.6f}")

# DMRG at same n=8
E_dmrg, S_dmrg, chi, dt = run_dmrg(3, 8, 1.0/3, chi_max=120)
print(f"  DMRG:  E={E_dmrg:.8f}, S_half={S_dmrg:.6f}, chi={chi}")
print(f"  Diff:  ΔE={abs(E_dmrg-vals[0]):.2e}, ΔS={abs(S_dmrg-S_exact):.6f}")

# === Large sizes with higher chi ===
results = []
sizes = [32, 48, 64]
chi_max = 120

total_start = time.time()
for n in sizes:
    if time.time() - total_start > 200:
        print(f"Time limit reached at n={n}", flush=True)
        break
    print(f"\n=== q=3, n={n}, g_c=0.333333, chi_max={chi_max} ===", flush=True)
    E, S, chi, dt = run_dmrg(3, n, 1.0/3, chi_max)
    print(f"  E={E:.8f}, S_half={S:.6f}, chi={chi}, t={dt:.1f}s", flush=True)
    results.append({'n': n, 'E': float(E), 'S_half': float(S), 'chi_actual': chi, 'time': dt})

    with open('results/sprint_054d_c_q3_larger.json', 'w') as f:
        json.dump({'sprint': '054d', 'q': 3, 'g_c': 1.0/3, 'chi_max': chi_max,
                   'exact_diag_n8': {'E': float(vals[0]), 'S_half': float(S_exact)},
                   'dmrg_n8': {'E': E_dmrg, 'S_half': S_dmrg},
                   'large_n_data': results}, f, indent=2)

# === Combine with chi=80 data and refit ===
# Previous data from 054a
prev = [
    {'n': 8, 'S_half': 0.516117},
    {'n': 12, 'S_half': 0.579230},
    {'n': 16, 'S_half': 0.622697},
    {'n': 24, 'S_half': 0.682416},
]
# Use chi=120 for large n
all_data = prev + results
ns = np.array([d['n'] for d in all_data])
Ss = np.array([d['S_half'] for d in all_data])
ln_n = np.log(ns)

A = np.vstack([ln_n, np.ones(len(ln_n))]).T
slope, intercept = np.linalg.lstsq(A, Ss, rcond=None)[0]
c_full = 6 * slope

print(f"\n=== Combined central charge (n={ns[0]}-{ns[-1]}) ===")
print(f"Full fit: c = {c_full:.4f} (CFT: 0.800)")
c_pairs = []
for i in range(len(ns) - 1):
    c_pair = 6 * (Ss[i+1] - Ss[i]) / (np.log(ns[i+1]) - np.log(ns[i]))
    c_pairs.append(c_pair)
    print(f"  n={ns[i]},{ns[i+1]}: c = {c_pair:.4f}")

# Save final
with open('results/sprint_054d_c_q3_larger.json', 'w') as f:
    json.dump({'sprint': '054d', 'q': 3, 'g_c': 1.0/3, 'chi_max': chi_max,
               'exact_diag_n8': {'E': float(vals[0]), 'S_half': float(S_exact)},
               'dmrg_n8': {'E': E_dmrg, 'S_half': S_dmrg},
               'large_n_data': results,
               'combined_analysis': {
                   'c_full_fit': float(c_full),
                   'c_pairwise': [float(c) for c in c_pairs],
                   'sizes': [int(n) for n in ns],
                   'S_halfs': [float(s) for s in Ss]}}, f, indent=2)

print("\nDone.", flush=True)
