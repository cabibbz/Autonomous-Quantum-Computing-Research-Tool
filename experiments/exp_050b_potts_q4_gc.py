#!/usr/bin/env python3
"""Sprint 050b: Find true q=4 Potts critical point via entropy FSS.

q=4 Potts (2D classical) is marginal — exactly at boundary of 1st/2nd order.
1D quantum: expected c=1 if second-order (from self-dual Potts CFT).
Use exact diag for n=4,6,8 (dim=4^n), DMRG for n=12,16.
"""
import numpy as np, json, time, warnings
warnings.filterwarnings('ignore')

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


def exact_diag_potts(q, n, g, J=1.0):
    """Full exact diag for q-state Potts."""
    dim = q**n
    H = np.zeros((dim, dim))
    for i in range(n - 1):
        for s in range(dim):
            digits = [(s // q**k) % q for k in range(n)]
            if digits[i] == digits[i + 1]:
                H[s, s] += -J
    for i in range(n):
        for s in range(dim):
            digits = [(s // q**k) % q for k in range(n)]
            new_d = digits.copy()
            new_d[i] = (digits[i] + 1) % q
            s2 = sum(new_d[k] * q**k for k in range(n))
            H[s, s2] += -g
            new_d[i] = (digits[i] - 1) % q
            s2 = sum(new_d[k] * q**k for k in range(n))
            H[s, s2] += -g

    from scipy.linalg import eigh
    evals, evecs = eigh(H)
    gs_vec = evecs[:, 0]
    gs_mat = gs_vec.reshape(q**(n // 2), q**(n // 2))
    rho = gs_mat @ gs_mat.T
    rho_evals = np.linalg.eigvalsh(rho)
    rho_evals = rho_evals[rho_evals > 1e-15]
    S = -np.sum(rho_evals * np.log(rho_evals))
    return evals[0], S


def run_dmrg(q, n, g, chi_max=80):
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
    S = float(psi.entanglement_entropy()[n // 2 - 1])
    chi_used = max(psi.chi)
    return float(E0), S, chi_used, time.time() - t0


# g scan — similar range as q=3 but shifted (q=4 g_c may differ)
g_values = [0.10, 0.15, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30,
            0.35, 0.40, 0.50, 0.70, 1.00]

results = {}

# === Exact diag: n=4,6,8 ===
for n in [4, 6, 8]:
    if 4**n > 100000:  # 4^8 = 65536, OK
        continue
    print(f"\n=== q=4 Potts n={n} (exact diag, dim={4**n}) ===", flush=True)
    results[n] = []
    for g in g_values:
        t0 = time.time()
        E, S = exact_diag_potts(4, n, float(g))
        dt = time.time() - t0
        print(f"  g={g:.2f}: S={S:.6f}, E={E:.4f}, t={dt:.1f}s", flush=True)
        results[n].append({'g': float(g), 'E': float(E), 'S_half': float(S),
                          'method': 'exact', 'time': dt})
        if dt > 60:
            print("  Too slow, stopping this size", flush=True)
            break

    with open('results/sprint_050b_q4_gc.json', 'w') as f:
        json.dump({'sprint': '050b', 'q': 4, 'data': {str(k): v for k, v in results.items()}}, f, indent=2)

# === DMRG: n=12,16 ===
# Time test first
print(f"\n=== Timing test: q=4 n=12, g=0.24, chi=80 ===", flush=True)
E, S, chi, dt = run_dmrg(4, 12, 0.24, chi_max=80)
print(f"  S={S:.6f}, chi={chi}, t={dt:.1f}s", flush=True)

chi_use = 80 if dt < 40 else 60 if dt < 80 else 40

for n in [12, 16]:
    chi = chi_use if n == 12 else min(chi_use, 60)
    print(f"\n=== q=4 Potts n={n}, chi={chi} (DMRG) ===", flush=True)
    results[n] = []
    t_start = time.time()
    for g in g_values:
        if time.time() - t_start > 200:
            print(f"  Time limit reached", flush=True)
            break
        E, S, chi_used, dt = run_dmrg(4, n, g, chi_max=chi)
        print(f"  g={g:.2f}: S={S:.6f}, chi={chi_used}, t={dt:.1f}s", flush=True)
        results[n].append({'g': float(g), 'E': float(E), 'S_half': float(S),
                          'chi_used': chi_used, 'method': 'DMRG', 'time': dt})

    with open('results/sprint_050b_q4_gc.json', 'w') as f:
        json.dump({'sprint': '050b', 'q': 4, 'data': {str(k): v for k, v in results.items()}}, f, indent=2)

# === Analysis ===
print("\n=== Central charge c(g) for q=4 Potts ===", flush=True)
all_data = {}
for n, data_list in results.items():
    for d in data_list:
        g_r = round(d['g'], 2)
        if g_r not in all_data:
            all_data[g_r] = {}
        all_data[g_r][n] = d['S_half']

print(f"\n{'g':>6s} | {'c':>6s} | sizes", flush=True)
print("-" * 50, flush=True)
c_values = []
for g in sorted(all_data.keys()):
    ns = sorted(all_data[g].keys())
    if len(ns) < 3:
        continue
    ns_a = np.array(ns, dtype=float)
    Ss_a = np.array([all_data[g][n] for n in ns])
    coeffs = np.polyfit(np.log(ns_a), Ss_a, 1)
    c = 6 * coeffs[0]
    print(f"  {g:.2f} | {c:>6.3f} | {ns}", flush=True)
    c_values.append((g, c))

if c_values:
    best_g, best_c = min(c_values, key=lambda x: abs(x[1] - 1.000))
    print(f"\n*** Best match to c=1: g={best_g:.2f}, c={best_c:.3f} ***", flush=True)

print("\nDone!", flush=True)
