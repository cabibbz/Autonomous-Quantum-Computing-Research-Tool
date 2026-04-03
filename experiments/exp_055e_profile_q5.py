#!/usr/bin/env python3
"""Sprint 055e: Central charge from entropy profile — Potts q=5 at g_c=0.441.
NO CFT prediction exists (2D classical q=5 is first-order).
Sprint 054 found c≈1.34. Profile method should confirm or deny c>1.
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
        J = model_params.get('J', 1.0); g = model_params.get('g', 1.0); q = model_params.get('q', 3)
        for a in range(q):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'Xphc')


def run_and_extract(q, n, g, chi_max=20):
    t0 = time.time()
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n + int(g * 1000))
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 30,
    })
    E0, _ = eng.run()
    S_profile = [float(s) for s in psi.entanglement_entropy()]
    chi_actual = max(psi.chi)
    dt = time.time() - t0
    return float(E0), S_profile, chi_actual, dt


def extract_c(S_profile, n):
    ls = np.arange(1, n)
    l_chord = (2 * n / np.pi) * np.sin(np.pi * ls / n)
    ln_chord = np.log(l_chord)
    S = np.array(S_profile)
    quarter = n // 4
    mask = (ls >= quarter) & (ls <= 3 * quarter)
    A = np.vstack([ln_chord[mask], np.ones(mask.sum())]).T
    slope, _ = np.linalg.lstsq(A, S[mask], rcond=None)[0]
    return float(6 * slope)


data = []
for n in [16, 24]:
    E0, S_prof, chi_act, dt = run_and_extract(5, n, 0.441, 20)
    c = extract_c(S_prof, n)
    print(f"n={n:3d}: c={c:.4f}, S_half={S_prof[n//2-1]:.6f}, chi={chi_act}, time={dt:.1f}s", flush=True)
    data.append({
        'n': n, 'E0': E0, 'S_half': S_prof[n//2-1],
        'c_profile': c, 'chi_actual': chi_act, 'time': dt,
        'S_profile': S_prof,
    })
    results = {'experiment': '055e', 'model': 'Potts', 'q': 5, 'g_c': 0.441,
               'c_cft': 'none (2D is first-order)', 'data': data}
    with open('results/exp_055e.json', 'w') as f:
        json.dump(results, f, indent=2)

print(f"\nc values: {' -> '.join(f'{d['c_profile']:.4f}' for d in data)}")
print("Sprint 054 FSS: c(raw)=1.335")
print("Saved results/exp_055e.json")
