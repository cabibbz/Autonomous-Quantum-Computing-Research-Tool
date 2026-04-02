#!/usr/bin/env python3
"""Sprint 081b: DMRG gap×N(q=6) at n=10,12 g_c=1/6 (S_q Potts).

Sprint 080 found gap×N INCREASES with q even as walking breaks.
Test: does q=6 gap×N stay stable or increase from n=8→12?
Exact diag gives gap×N at n=6,7,8 (Sprint 080). Extend with DMRG.

Uses orthogonal_to for excited state (same pattern as Sprint 078b).
"""
import numpy as np
import json, time

q = 6
gc = 1.0 / q

results = {
    'experiment': '081b_sq_dmrg_gap_q6',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'gc': gc,
    'data': [],
}

def save():
    with open("results/sprint_081b_sq_dmrg_gap_q6.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

class SqPottsSite(Site):
    def __init__(self, q):
        leg = npc.LegCharge.from_trivial(q)
        Site.__init__(self, leg, [str(a) for a in range(q)], sort_charge=False)
        for a in range(q):
            P = np.zeros((q, q), dtype=complex); P[a, a] = 1.0
            self.add_op(f'P{a}', P)
        Sq_field = np.ones((q, q), dtype=complex) - np.eye(q, dtype=complex)
        self.add_op('SqField', Sq_field, hc='SqField')

class SqPottsChain(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, model_params):
        return SqPottsSite(model_params.get('q', 6))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', gc)
        q_val = model_params.get('q', 6)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')

def dmrg_gap(n, q_val, g, chi_max=40):
    model = SqPottsChain({'L': n, 'q': q_val, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})

    # Ground state
    np.random.seed(42 + n)
    init = [np.random.randint(q_val) for _ in range(n)]
    psi0 = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng0 = dmrg.TwoSiteDMRGEngine(psi0, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 30,
    })
    E0, _ = eng0.run()

    # Excited state
    np.random.seed(123 + n)
    init1 = [(np.random.randint(q_val) + 1) % q_val for _ in range(n)]
    psi1 = MPS.from_product_state(model.lat.mps_sites(), init1, bc='finite')
    eng1 = dmrg.TwoSiteDMRGEngine(psi1, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 40,
    }, orthogonal_to=[psi0])
    E1, _ = eng1.run()

    return float(E0), float(E1)

# ---- MAIN ----
print(f"Sprint 081b: S_q Potts DMRG gap×N at q={q}, g_c=1/{q}={gc:.6f}")
print("=" * 60, flush=True)

# Reference: exact diag from Sprint 080
print("\nExact diag reference (Sprint 080):")
print("  n=6: gap×N ≈ 2.01 (from c_eff compilation)")
print("  n=8: gap×N from exact diag")

runs = [(10, 40), (12, 30)]

for n, chi in runs:
    print(f"\n--- n={n}, chi_max={chi} ---", flush=True)
    t0 = time.time()
    E0, E1 = dmrg_gap(n, q, gc, chi_max=chi)
    dt = time.time() - t0

    gap = E1 - E0
    gapN = gap * n
    print(f"  {dt:.1f}s, E0={E0:.8f}, E1={E1:.8f}", flush=True)
    print(f"  gap = {gap:.6f}, gap×N = {gapN:.4f}", flush=True)

    entry = {
        'n': n, 'chi_max': chi, 'E0': E0, 'E1': E1,
        'gap': gap, 'gapN': gapN, 'time_s': dt,
    }
    results['data'].append(entry)
    save()
    print(f"  Saved.", flush=True)

    if dt > 250:
        print(f"  Approaching limit ({dt:.0f}s), stopping.", flush=True)
        break

# Summary
print(f"\n{'='*60}", flush=True)
print("SUMMARY: gap×N vs n for q=6 S_q Potts", flush=True)
print(f"{'='*60}", flush=True)
for d in results['data']:
    print(f"  n={d['n']:3d}: gap={d['gap']:.6f}, gap×N={d['gapN']:.4f} ({d['time_s']:.0f}s)", flush=True)

if len(results['data']) >= 2:
    gN = [d['gapN'] for d in results['data']]
    drift = 100 * (gN[-1] - gN[0]) / gN[0]
    print(f"\ngap×N drift: {gN[0]:.4f} → {gN[-1]:.4f} ({drift:+.1f}%)", flush=True)
    results['gapN_drift_pct'] = drift

save()
print(f"\nFinal save to results/sprint_081b_sq_dmrg_gap_q6.json", flush=True)

from db_utils import record
for d in results['data']:
    record(sprint=81, model='sq_potts', q=q, n=d['n'],
           quantity='gapN_dmrg', value=d['gapN'],
           method=f'dmrg_chi{d["chi_max"]}',
           notes=f'gc=1/{q}, open BC, orthogonal_to')
    record(sprint=81, model='sq_potts', q=q, n=d['n'],
           quantity='gap', value=d['gap'],
           method=f'dmrg_chi{d["chi_max"]}',
           notes=f'gc=1/{q}')
print("Recorded to DB.", flush=True)
