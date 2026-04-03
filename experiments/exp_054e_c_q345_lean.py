#!/usr/bin/env python3
"""Sprint 054e: Central charge c(q=3,4,5) — lean version.

Use chi=60, n=8,12,16,24,32. Run sequentially to avoid CPU contention.
q=3 cross-checks against chi=80 data.
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


def run_dmrg(q, n, g, chi_max=60):
    t0 = time.time()
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n + int(g * 1000))
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
configs = [
    (3, 1.0/3, 0.800, "3-state Potts CFT"),
    (4, 0.392, 1.000, "Ashkin-Teller/free boson"),
    (5, 0.441, None, "NO 2D CFT prediction"),
]

grand_start = time.time()
for q, g_c, c_pred, label in configs:
    if time.time() - grand_start > 240:
        print(f"\nGlobal time limit, stopping before q={q}", flush=True)
        break

    print(f"\n{'='*60}")
    print(f"q={q}, g_c={g_c:.4f}, c_pred={c_pred} ({label})")
    print(f"{'='*60}", flush=True)

    data = []
    sizes = [8, 12, 16, 24, 32]
    q_start = time.time()
    for n in sizes:
        if time.time() - q_start > 70:  # 70s per q
            print(f"  Time limit for q={q}, stopping at n={n}", flush=True)
            break
        t0 = time.time()
        E, S, chi, dt = run_dmrg(q, n, g_c, chi_max=60)
        print(f"  n={n}: S={S:.6f}, E={E:.8f}, chi={chi}, t={dt:.1f}s", flush=True)
        data.append({'n': n, 'E': E, 'S_half': S, 'chi_actual': chi, 'time': dt})

        # Save incrementally
        all_results[str(q)] = {'g_c': g_c, 'c_prediction': c_pred, 'label': label,
                               'chi_max': 60, 'data': data}
        with open('results/sprint_054e_central_charges.json', 'w') as f:
            json.dump({'sprint': '054e', 'results': all_results}, f, indent=2)

    if len(data) >= 3:
        ns = np.array([d['n'] for d in data])
        Ss = np.array([d['S_half'] for d in data])
        c_full, c_pairs = extract_c(ns, Ss, f"q={q}")
        all_results[str(q)]['c_full_fit'] = c_full
        all_results[str(q)]['c_pairwise'] = c_pairs

# Final save
with open('results/sprint_054e_central_charges.json', 'w') as f:
    json.dump({'sprint': '054e', 'results': all_results}, f, indent=2)

# Summary
print(f"\n{'='*60}")
print("SUMMARY")
print(f"{'='*60}")
for q_str, res in all_results.items():
    c_fit = res.get('c_full_fit', '?')
    c_pred = res.get('c_prediction', '?')
    print(f"  q={q_str}: c_fit={c_fit:.4f} vs c_pred={c_pred}")

print("\nDone.", flush=True)
