#!/usr/bin/env python3
"""Sprint 049c: Continue entropy FSS — q=2 (n=32,48,64), q=3 (n=12,16,24), q=4 (n=8,12,16).

Loads existing results from sprint_049b_entropy_fss.json and appends new data.
"""
import numpy as np
import json, time, sys
import warnings
warnings.filterwarnings('ignore')

from tenpy.models.tf_ising import TFIChain
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


def run_entropy(q, n, g, chi_max=40):
    t0 = time.time()
    if q == 2:
        model = TFIChain({'L': n, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    else:
        model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    init = [i % max(q, 2) for i in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 30,
    })
    E0, _ = eng.run()
    dt = time.time() - t0
    S_half = float(psi.entanglement_entropy()[n//2 - 1])
    return float(E0), S_half, dt


# Load existing results
with open('results/sprint_049b_entropy_fss.json') as f:
    all_results = json.load(f)

def save():
    with open('results/sprint_049b_entropy_fss.json', 'w') as f:
        json.dump(all_results, f, indent=2)

# === Timing test ===
print("=== Timing: q=2 n=48 g=1.0 chi=40 ===", flush=True)
_, S, dt = run_entropy(2, 48, 1.0, 40)
print(f"  S={S:.4f}, t={dt:.1f}s", flush=True)

print("=== Timing: q=3 n=24 g=1.0 chi=40 ===", flush=True)
_, S, dt = run_entropy(3, 24, 1.0, 40)
print(f"  S={S:.4f}, t={dt:.1f}s", flush=True)

print("=== Timing: q=4 n=16 g=0.89 chi=40 ===", flush=True)
_, S, dt = run_entropy(4, 16, 0.89, 40)
print(f"  S={S:.4f}, t={dt:.1f}s", flush=True)

# ================================================================
# q=2 TFIM continuation: n=32, 48, 64
# ================================================================
g_values_q2 = [0.70, 0.80, 0.85, 0.90, 0.93, 0.95, 0.97, 0.99,
               1.00, 1.01, 1.03, 1.05, 1.07, 1.10, 1.15, 1.20, 1.30, 1.50]

for n in [32, 48, 64]:
    key = str(n)
    if key in all_results.get('q2', {}):
        print(f"\nq=2 n={n}: already done, skipping", flush=True)
        continue
    print(f"\n=== q=2 n={n} ===", flush=True)
    n_data = []
    t_start = time.time()
    for g in g_values_q2:
        if time.time() - t_start > 240:
            print(f"  Time limit for n={n}", flush=True)
            break
        E0, S, dt = run_entropy(2, n, g, chi_max=40)
        print(f"  g={g:.2f}: S={S:.4f}, t={dt:.1f}s", flush=True)
        n_data.append({'g': g, 'E0': E0, 'S_half': S, 'time': dt})
    if 'q2' not in all_results:
        all_results['q2'] = {}
    all_results['q2'][key] = n_data
    save()

# ================================================================
# q=3 Potts: n=12, 16, 24
# ================================================================
g_values_q3 = [0.70, 0.80, 0.85, 0.90, 0.93, 0.95, 0.97, 0.99,
               1.00, 1.01, 1.03, 1.05, 1.10, 1.20, 1.40]

for n in [12, 16, 24]:
    key = str(n)
    if key in all_results.get('q3', {}) and len(all_results['q3'][key]) >= 10:
        print(f"\nq=3 n={n}: already done, skipping", flush=True)
        continue
    print(f"\n=== q=3 n={n} ===", flush=True)
    n_data = []
    t_start = time.time()
    for g in g_values_q3:
        if time.time() - t_start > 240:
            print(f"  Time limit for n={n}", flush=True)
            break
        E0, S, dt = run_entropy(3, n, g, chi_max=40)
        print(f"  g={g:.2f}: S={S:.4f}, t={dt:.1f}s", flush=True)
        n_data.append({'g': g, 'E0': E0, 'S_half': S, 'time': dt})
    if 'q3' not in all_results:
        all_results['q3'] = {}
    all_results['q3'][key] = n_data
    save()

# ================================================================
# q=4 Potts: n=8, 12, 16
# ================================================================
g_values_q4 = [0.50, 0.60, 0.70, 0.75, 0.80, 0.83, 0.86, 0.89,
               0.92, 0.95, 1.00, 1.10, 1.20]

for n in [8, 12, 16]:
    key = str(n)
    if key in all_results.get('q4', {}) and len(all_results['q4'][key]) >= 10:
        print(f"\nq=4 n={n}: already done, skipping", flush=True)
        continue
    print(f"\n=== q=4 n={n} ===", flush=True)
    n_data = []
    t_start = time.time()
    for g in g_values_q4:
        if time.time() - t_start > 240:
            print(f"  Time limit for n={n}", flush=True)
            break
        E0, S, dt = run_entropy(4, n, g, chi_max=40)
        print(f"  g={g:.2f}: S={S:.4f}, t={dt:.1f}s", flush=True)
        n_data.append({'g': g, 'E0': E0, 'S_half': S, 'time': dt})
    if 'q4' not in all_results:
        all_results['q4'] = {}
    all_results['q4'][key] = n_data
    save()

print("\n=== All entropy FSS data saved ===", flush=True)
