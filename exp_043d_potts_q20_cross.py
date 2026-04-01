#!/usr/bin/env python3
"""Sprint 043d: q=20 Potts — verify model + entropy + crossing test.
Minimal version: chi=10, n=8 only, few g points.
"""
import numpy as np
import numpy.linalg as la
import json, time, sys

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

class PottsModel(CouplingMPOModel):
    def init_sites(self, model_params):
        return PottsSite(model_params.get('q', 20))
    def init_terms(self, model_params):
        J, g, q = model_params.get('J', 1.0), model_params.get('g', 1.0), model_params.get('q', 20)
        for a in range(q):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'Xphc')

class PottsChain(PottsModel, NearestNeighborModel):
    pass

t_start = time.time()

# ── Verify: q=10 vs q=20 have DIFFERENT energies ──
print("=== Verification: q=10 vs q=20 give different E0 ===")
for q_test in [10, 20]:
    model = PottsChain({'L': 8, 'q': q_test, 'J': 1.0, 'g': 0.15, 'bc_MPS': 'finite'})
    d = model.lat.site(0).dim
    print(f"  q={q_test}: site dim = {d}")
    psi = MPS.from_product_state(model.lat.mps_sites(), [0]*8, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-5,
        'trunc_params': {'chi_max': 10, 'svd_min': 1e-8},
        'max_sweeps': 6,
    })
    E0, _ = eng.run()
    S_half = psi.entanglement_entropy()[3]  # entropy at bond 3-4 (half chain)
    chi = max(psi.chi)
    print(f"  q={q_test}: E0={E0:.6f}, S_half={S_half:.6f}, chi_max_used={chi}")
    sys.stdout.flush()

print(f"  (If E0 values differ, model is correct)")

# ── q=20 entropy sweep (cheap — no MI needed) ──
print("\n=== q=20 half-chain entropy sweep, n=8 ===")
q = 20
n = 8
g_values = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
results = []

for g in g_values:
    t0 = time.time()
    try:
        model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
        psi = MPS.from_product_state(model.lat.mps_sites(), [0]*n, bc='finite')
        eng = dmrg.TwoSiteDMRGEngine(psi, model, {
            'mixer': True, 'max_E_err': 1e-5,
            'trunc_params': {'chi_max': 10, 'svd_min': 1e-8},
            'max_sweeps': 6,
        })
        E0, _ = eng.run()
        S_half = psi.entanglement_entropy()[n//2 - 1]
        chi = max(psi.chi)
        dt = time.time() - t0
        print(f"  g={g:.2f}: S_half={S_half:.4f}, E0={E0:.6f}, chi={chi}, time={dt:.1f}s")
        results.append({'g': g, 'S_half': float(S_half), 'E0': E0, 'chi': chi, 'time': dt})
    except Exception as e:
        dt = time.time() - t0
        print(f"  g={g:.2f}: FAILED ({type(e).__name__}), time={dt:.1f}s")
        results.append({'g': g, 'status': 'failed', 'time': dt})
    sys.stdout.flush()
    with open('results/sprint_043d_q20_entropy.json', 'w') as f:
        json.dump({'sprint': '043d', 'q': q, 'n': n, 'chi_max': 10,
                   'data': results, 'total_time': time.time()-t_start}, f, indent=2)

print(f"\nTotal: {time.time()-t_start:.1f}s")

# Derivative analysis
ok = [r for r in results if 'S_half' in r]
if len(ok) >= 2:
    print("\n=== ENTROPY DERIVATIVE ===")
    for i in range(len(ok)-1):
        dsdg = (ok[i+1]['S_half'] - ok[i]['S_half']) / (ok[i+1]['g'] - ok[i]['g'])
        print(f"  g={ok[i]['g']:.2f}-{ok[i+1]['g']:.2f}: dS/dg = {dsdg:.2f}")
