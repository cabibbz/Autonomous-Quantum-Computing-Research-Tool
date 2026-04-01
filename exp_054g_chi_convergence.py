#!/usr/bin/env python3
"""Sprint 054g: Chi convergence test for entropy at q=3 critical point.

Check if chi=60/80 is sufficient for entropy at g_c. If S changes with chi,
our c estimates are chi-limited, not size-limited.
Test at n=16 (fixed) with chi=20,40,60,80,100,120.
Also test q=4 n=12 to see chi sensitivity at larger d.
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


def run_dmrg(q, n, g, chi_max):
    t0 = time.time()
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n + int(g * 1000) + chi_max)
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 40,
    })
    E0, _ = eng.run()
    S_half = float(psi.entanglement_entropy()[n // 2 - 1])
    chi_actual = max(psi.chi)
    return float(E0), S_half, chi_actual, time.time() - t0


results = {}

# === q=3, n=16 chi convergence ===
print("=== q=3, n=16, g_c=1/3: chi convergence ===", flush=True)
q3_data = []
for chi in [20, 40, 60, 80]:
    E, S, chi_act, dt = run_dmrg(3, 16, 1.0/3, chi)
    print(f"  chi={chi}: S={S:.6f}, E={E:.8f}, chi_act={chi_act}, t={dt:.1f}s", flush=True)
    q3_data.append({'chi_max': chi, 'S_half': S, 'E': E, 'chi_actual': chi_act, 'time': dt})
    results['q3_n16'] = q3_data
    with open('results/sprint_054g_chi_convergence.json', 'w') as f:
        json.dump({'sprint': '054g', 'results': results}, f, indent=2)

# === q=4, n=12 chi convergence ===
print("\n=== q=4, n=12, g_c=0.392: chi convergence ===", flush=True)
q4_data = []
for chi in [20, 40, 60]:
    E, S, chi_act, dt = run_dmrg(4, 12, 0.392, chi)
    print(f"  chi={chi}: S={S:.6f}, E={E:.8f}, chi_act={chi_act}, t={dt:.1f}s", flush=True)
    q4_data.append({'chi_max': chi, 'S_half': S, 'E': E, 'chi_actual': chi_act, 'time': dt})
    results['q4_n12'] = q4_data
    with open('results/sprint_054g_chi_convergence.json', 'w') as f:
        json.dump({'sprint': '054g', 'results': results}, f, indent=2)

print("\n=== Summary ===")
print("q=3 n=16:", [(d['chi_max'], f"S={d['S_half']:.6f}") for d in q3_data])
print("q=4 n=12:", [(d['chi_max'], f"S={d['S_half']:.6f}") for d in q4_data])
print("Done.", flush=True)
