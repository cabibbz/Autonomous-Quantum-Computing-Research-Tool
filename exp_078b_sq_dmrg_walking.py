#!/usr/bin/env python3
"""Sprint 078b: DMRG for S_q Potts at q=5 — find walking length.

S_q Potts at q>4 is first-order (Baxter), but looks CFT-like at n≤8.
DMRG at exact g_c = 1/5 = 0.2 for n=8-24 to find walking breakdown.
Open BC, S_q field = Σ_{k=1}^{q-1} X^k.

Fix: orthogonal_to must be keyword arg, not in config dict.
"""
import numpy as np
import json, time

q = 5
gc = 0.2  # 1/q exact

results = {
    'experiment': '078b_sq_dmrg_walking',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'gc': gc,
    'gap_data': [],
}

def save():
    with open("results/sprint_078b_sq_dmrg_walking.json", "w") as f:
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
        return SqPottsSite(model_params.get('q', 5))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 0.2)
        q_val = model_params.get('q', 5)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')

def dmrg_gap(n, q_val, g, chi_max=60):
    model = SqPottsChain({'L': n, 'q': q_val, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})

    # Ground state
    np.random.seed(42 + n)
    init = [np.random.randint(q_val) for _ in range(n)]
    psi0 = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng0 = dmrg.TwoSiteDMRGEngine(psi0, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 40,
    })
    E0, _ = eng0.run()
    S_profile = [float(psi0.entanglement_entropy()[i]) for i in range(n - 1)]

    # Excited state — orthogonal_to is keyword-only arg
    np.random.seed(123 + n)
    init1 = [(np.random.randint(q_val) + 1) % q_val for _ in range(n)]
    psi1 = MPS.from_product_state(model.lat.mps_sites(), init1, bc='finite')
    eng1 = dmrg.TwoSiteDMRGEngine(psi1, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 50,
    }, orthogonal_to=[psi0])
    E1, _ = eng1.run()
    return float(E0), float(E1), S_profile

# ---- MAIN ----
print("=" * 60, flush=True)
print(f"Sprint 078b: S_q Potts DMRG walking length at q={q}, g_c={gc}", flush=True)
print("=" * 60, flush=True)

sizes = [8, 12, 16, 20, 24]
chi_map = {8: 40, 12: 60, 16: 80, 20: 80, 24: 100}

for n in sizes:
    chi = chi_map[n]
    print(f"\n--- n={n}, chi={chi} ---", flush=True)
    t0 = time.time()
    E0, E1, S = dmrg_gap(n, q, gc, chi_max=chi)
    dt = time.time() - t0
    gap = E1 - E0
    gapN = gap * n
    S_mid = S[n // 2 - 1]
    print(f"  {dt:.1f}s, E0={E0:.8f}, E1={E1:.8f}", flush=True)
    print(f"  gap = {gap:.8f}, gap×N = {gapN:.4f}, S_mid = {S_mid:.4f}", flush=True)

    entry = {'n': n, 'chi': chi, 'E0': E0, 'E1': E1,
             'gap': gap, 'gapN': gapN, 'S_mid': S_mid,
             'S_profile': S, 'time_s': dt}
    results['gap_data'].append(entry)
    save()
    print(f"  Saved.", flush=True)

    if dt > 250:
        print(f"  Approaching time limit, stopping.", flush=True)
        break

# Summary
print(f"\n{'='*60}", flush=True)
print("SUMMARY: gap×N vs system size", flush=True)
print(f"{'='*60}", flush=True)
print(f"{'n':>4} {'chi':>4} {'gap':>12} {'gap×N':>10} {'S_mid':>8} {'time':>8}", flush=True)
for d in results['gap_data']:
    print(f"{d['n']:4d} {d['chi']:4d} {d['gap']:12.6f} {d['gapN']:10.4f} {d['S_mid']:8.4f} {d['time_s']:8.1f}s", flush=True)

# Walking analysis
gN_vals = [d['gapN'] for d in results['gap_data']]
if len(gN_vals) >= 2:
    gN_first = gN_vals[0]
    for d in results['gap_data']:
        if d['gapN'] < 0.5 * gN_first and gN_first > 0.01:
            print(f"\n  Walking breakdown at n={d['n']}: gap×N = {d['gapN']:.4f}", flush=True)
            results['walking_breakdown_n'] = d['n']
            break
    else:
        print(f"\n  No walking breakdown up to n={results['gap_data'][-1]['n']}", flush=True)
        results['walking_breakdown_n'] = None

    # Entropy scaling: S_mid ~ (c/6)ln(n) for open BC
    if len(results['gap_data']) >= 3:
        ns = np.array([d['n'] for d in results['gap_data']])
        Ss = np.array([d['S_mid'] for d in results['gap_data']])
        ln_n = np.log(ns)
        # Linear fit: S = (c/6)ln(n) + const
        A = np.vstack([ln_n, np.ones_like(ln_n)]).T
        slope, intercept = np.linalg.lstsq(A, Ss, rcond=None)[0]
        c_eff = 6 * slope
        print(f"\n  Entropy scaling: S_mid = {slope:.4f}·ln(n) + {intercept:.4f}", flush=True)
        print(f"  c_eff = {c_eff:.4f}", flush=True)
        results['c_eff'] = c_eff

save()
print("\nFinal save to results/sprint_078b_sq_dmrg_walking.json", flush=True)

from db_utils import record
for d in results['gap_data']:
    record(sprint=78, model='sq_potts', q=q, n=d['n'],
           quantity='gapN_dmrg', value=d['gapN'],
           method=f'dmrg_chi{d["chi"]}', notes=f'gc=1/{q}, open BC, orthogonal_to fixed')
    record(sprint=78, model='sq_potts', q=q, n=d['n'],
           quantity='S_mid', value=d['S_mid'],
           method=f'dmrg_chi{d["chi"]}', notes=f'gc=1/{q}, midchain entropy')
if results.get('c_eff'):
    record(sprint=78, model='sq_potts', q=q, n=0,
           quantity='c_eff', value=results['c_eff'],
           method='dmrg_entropy_fit', notes=f'from S_mid ~ (c/6)ln(n)')
print("Recorded to DB.", flush=True)
