#!/usr/bin/env python3
"""Sprint 054b: Central charge c(q=4) at true g_c=0.392.

Extract c from entropy scaling S(n/2) = (c/6)ln(n) + const.
CFT prediction: c = 1 (Ashkin-Teller/free boson, marginal).
Use DMRG (PottsChain) at g_c=0.392 for n=8,12,16,24,32.
q=4 is marginal — logarithmic corrections expected.
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
    S_half = float(psi.entanglement_entropy()[n // 2 - 1])
    chi_actual = max(psi.chi)
    return float(E0), S_half, chi_actual, time.time() - t0


q = 4
g_c = 0.392  # from energy gap Sprint 051
results = []
sizes = [8, 12, 16, 24, 32, 48]
chi_max = 80

total_start = time.time()
for n in sizes:
    if time.time() - total_start > 240:
        print(f"Time limit reached, stopping at n={n}", flush=True)
        break
    print(f"=== q={q}, n={n}, g_c={g_c}, chi_max={chi_max} ===", flush=True)
    E, S, chi, dt = run_dmrg(q, n, g_c, chi_max)
    print(f"  E={E:.8f}, S_half={S:.6f}, chi={chi}, t={dt:.1f}s", flush=True)
    results.append({'n': n, 'E': float(E), 'S_half': float(S), 'chi_actual': chi, 'time': dt})

    with open('results/sprint_054b_central_charge_q4.json', 'w') as f:
        json.dump({'sprint': '054b', 'q': q, 'g_c': g_c, 'chi_max': chi_max,
                   'data': results}, f, indent=2)

# === Extract central charge ===
if len(results) >= 3:
    ns = np.array([r['n'] for r in results])
    Ss = np.array([r['S_half'] for r in results])
    ln_n = np.log(ns)

    A = np.vstack([ln_n, np.ones(len(ln_n))]).T
    slope, intercept = np.linalg.lstsq(A, Ss, rcond=None)[0]
    c_full = 6 * slope

    print(f"\n=== Central charge extraction ===")
    print(f"Full fit: c = {c_full:.4f} (CFT: 1.000)")
    print(f"Pairwise estimates:")
    c_pairs = []
    for i in range(len(ns) - 1):
        c_pair = 6 * (Ss[i+1] - Ss[i]) / (np.log(ns[i+1]) - np.log(ns[i]))
        c_pairs.append(c_pair)
        print(f"  n={ns[i]},{ns[i+1]}: c = {c_pair:.4f}")

    with open('results/sprint_054b_central_charge_q4.json', 'w') as f:
        json.dump({'sprint': '054b', 'q': q, 'g_c': g_c, 'chi_max': chi_max,
                   'data': results, 'analysis': {
                       'c_full_fit': float(c_full),
                       'c_pairwise': [float(c) for c in c_pairs],
                       'c_exact': 1.000}}, f, indent=2)

print("\nDone.", flush=True)
