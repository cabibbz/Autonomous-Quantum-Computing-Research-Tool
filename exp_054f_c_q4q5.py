#!/usr/bin/env python3
"""Sprint 054f: Central charge c(q=4,5) with reduced chi.

q=3 already converging (c_pairwise: 0.934→0.907→0.884 toward 0.800).
Focus on q=4 and q=5. Use chi=40 for speed.
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


def run_dmrg(q, n, g, chi_max=40):
    t0 = time.time()
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n + int(g * 1000))
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 30,
    })
    E0, _ = eng.run()
    S_half = float(psi.entanglement_entropy()[n // 2 - 1])
    chi_actual = max(psi.chi)
    return float(E0), S_half, chi_actual, time.time() - t0


def extract_c(ns, Ss, label=""):
    ln_n = np.log(ns)
    A = np.vstack([ln_n, np.ones(len(ln_n))]).T
    slope, intercept = np.linalg.lstsq(A, Ss, rcond=None)[0]
    c_full = 6 * slope
    print(f"\n  {label} Full fit: c = {c_full:.4f}")
    c_pairs = []
    for i in range(len(ns) - 1):
        c_pair = 6 * (Ss[i+1] - Ss[i]) / (np.log(ns[i+1]) - np.log(ns[i]))
        c_pairs.append(float(c_pair))
        print(f"    n={ns[i]},{ns[i+1]}: c = {c_pair:.4f}")
    return float(c_full), c_pairs


all_results = {}

# === q=4 ===
print("=" * 60)
print("q=4, g_c=0.392, c_pred=1.0 (Ashkin-Teller)")
print("=" * 60, flush=True)
data4 = []
for n in [8, 12, 16, 24]:
    E, S, chi, dt = run_dmrg(4, n, 0.392, chi_max=40)
    print(f"  n={n}: S={S:.6f}, E={E:.8f}, chi={chi}, t={dt:.1f}s", flush=True)
    data4.append({'n': n, 'E': E, 'S_half': S, 'chi_actual': chi, 'time': dt})
    all_results['4'] = {'g_c': 0.392, 'data': data4}
    with open('results/sprint_054f_c_q4q5.json', 'w') as f:
        json.dump({'sprint': '054f', 'results': all_results}, f, indent=2)

ns4 = np.array([d['n'] for d in data4])
Ss4 = np.array([d['S_half'] for d in data4])
c4, c4_pairs = extract_c(ns4, Ss4, "q=4")
all_results['4']['c_full_fit'] = c4
all_results['4']['c_pairwise'] = c4_pairs

# === q=5 ===
print("\n" + "=" * 60)
print("q=5, g_c=0.441, NO CFT prediction")
print("=" * 60, flush=True)
data5 = []
for n in [8, 12, 16, 24]:
    E, S, chi, dt = run_dmrg(5, n, 0.441, chi_max=40)
    print(f"  n={n}: S={S:.6f}, E={E:.8f}, chi={chi}, t={dt:.1f}s", flush=True)
    data5.append({'n': n, 'E': E, 'S_half': S, 'chi_actual': chi, 'time': dt})
    all_results['5'] = {'g_c': 0.441, 'data': data5}
    with open('results/sprint_054f_c_q4q5.json', 'w') as f:
        json.dump({'sprint': '054f', 'results': all_results}, f, indent=2)

ns5 = np.array([d['n'] for d in data5])
Ss5 = np.array([d['S_half'] for d in data5])
c5, c5_pairs = extract_c(ns5, Ss5, "q=5")
all_results['5']['c_full_fit'] = c5
all_results['5']['c_pairwise'] = c5_pairs

# Final save
with open('results/sprint_054f_c_q4q5.json', 'w') as f:
    json.dump({'sprint': '054f', 'results': all_results}, f, indent=2)

print(f"\n{'='*60}")
print("SUMMARY")
print(f"  q=4: c = {c4:.4f} (CFT: 1.000)")
print(f"  q=5: c = {c5:.4f} (no prediction)")
print(f"{'='*60}")
print("Done.", flush=True)
