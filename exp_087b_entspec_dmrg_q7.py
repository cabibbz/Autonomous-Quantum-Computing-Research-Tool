#!/usr/bin/env python3
"""Sprint 087b: Entanglement spectrum scaling via DMRG at q=7, n=6,8,10,12,16.

Track: λ_max, (q-1)-fold multiplet fraction, tail weight, entanglement gap.
Compare degradation rate to q=5 from 087a.
q=7 is in walking-broken regime — expect faster tail growth.
"""
import numpy as np
import json, time

q = 7
gc = 1.0 / q

results = {
    'experiment': '087b_entspec_dmrg_q7',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'gc': gc,
    'data': [],
}

def save():
    with open("results/sprint_087b_entspec_dmrg_q7.json", "w") as f:
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
        return SqPottsSite(model_params.get('q', q))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', gc)
        q_val = model_params.get('q', q)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')

def analyze_spectrum(spec, q_val):
    """Analyze entanglement spectrum: λ_max, multiplet, tail."""
    spec = np.sort(spec)[::-1]
    spec = spec[spec > 1e-30]

    lam_max = spec[0]
    multiplet = spec[1:q_val] if len(spec) >= q_val else spec[1:]
    tail = spec[q_val:] if len(spec) > q_val else np.array([])

    S_total = -np.sum(spec * np.log(np.clip(spec, 1e-30, None)))
    s_i = -spec * np.log(np.clip(spec, 1e-30, None))
    S_lev0 = s_i[0] if len(s_i) > 0 else 0
    S_lev1 = np.sum(s_i[1:q_val]) if len(s_i) >= q_val else np.sum(s_i[1:])
    S_tail = np.sum(s_i[q_val:]) if len(s_i) > q_val else 0

    w_max = lam_max
    w_mult = np.sum(multiplet)
    w_tail = np.sum(tail)

    xi = -np.log(np.clip(spec, 1e-30, None))
    ent_gap = xi[1] - xi[0] if len(xi) > 1 else 0

    if len(multiplet) > 0:
        mult_spread = (np.max(multiplet) - np.min(multiplet)) / np.mean(multiplet)
    else:
        mult_spread = 0

    return {
        'lam_max': float(lam_max),
        'multiplet_mean': float(np.mean(multiplet)) if len(multiplet) > 0 else 0,
        'multiplet_spread': float(mult_spread),
        'w_max': float(w_max),
        'w_mult': float(w_mult),
        'w_tail': float(w_tail),
        'S_total': float(S_total),
        'S_lev0_frac': float(S_lev0 / S_total) if S_total > 0 else 0,
        'S_lev1_frac': float(S_lev1 / S_total) if S_total > 0 else 0,
        'S_tail_frac': float(S_tail / S_total) if S_total > 0 else 0,
        'ent_gap': float(ent_gap),
        'schmidt_rank': len(spec),
        'n_tail': len(tail),
    }

# Complex CFT Re(c) for q=7
alpha_cft = np.arccosh(np.sqrt(q) / 2)
Re_c = 1 + 6 * alpha_cft**2 / (np.pi**2 + alpha_cft**2)
results['Re_c'] = float(Re_c)

print(f"Sprint 087b: DMRG entanglement spectrum scaling for S_q Potts q={q}")
print(f"g_c = 1/{q} = {gc:.6f}, Re(c) = {Re_c:.4f}")
print("=" * 70, flush=True)

# q=7, d=7 per site. Sprint 086b: n=8 chi=50 took 172s. Be conservative.
# n=6: fast. n=8: ~170s. n=10: chi=40 ~200s. n=12: chi=30 ~200s. n=16: risky.
runs = [(6, 50), (8, 56), (10, 50), (12, 40)]

total_time = 0
for n, chi_max in runs:
    print(f"\n{'='*60}")
    print(f"n={n}, chi_max={chi_max}", flush=True)
    t0 = time.time()

    model = SqPottsChain({'L': n, 'q': q, 'J': 1.0, 'g': gc, 'bc_MPS': 'finite'})
    np.random.seed(42 + n)
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 30,
    })
    E0, _ = eng.run()
    dt_dmrg = time.time() - t0
    total_time += dt_dmrg
    chi_used = max(psi.chi)

    print(f"  DMRG: {dt_dmrg:.1f}s, E0={E0:.8f}, chi_used={chi_used}", flush=True)

    mid_bond = n // 2 - 1
    sv_mid = psi.get_SL(mid_bond)
    spec_mid = sv_mid**2

    analysis = analyze_spectrum(spec_mid, q)

    print(f"  Schmidt rank: {analysis['schmidt_rank']}, ent_gap: {analysis['ent_gap']:.4f}")
    print(f"  λ_max: {analysis['lam_max']:.6f}, w_mult: {analysis['w_mult']:.6f}, w_tail: {analysis['w_tail']:.6f}")
    print(f"  S = {analysis['S_total']:.6f}")
    print(f"  Entropy fracs: lev0={analysis['S_lev0_frac']:.4f}, lev1={analysis['S_lev1_frac']:.4f}, tail={analysis['S_tail_frac']:.4f}")

    entry = {
        'n': n, 'chi_max': chi_max, 'chi_used': chi_used,
        'E0': float(E0), 'time_dmrg': dt_dmrg,
        'spec_mid_top20': [float(x) for x in spec_mid[:20]],
        **analysis,
    }
    results['data'].append(entry)
    save()
    print(f"  Saved. (cumulative: {total_time:.0f}s)", flush=True)

    # Don't attempt next size if this one was too slow
    if total_time > 500:
        print(f"  Stopping: cumulative time {total_time:.0f}s.", flush=True)
        break

# Summary table
print(f"\n{'='*70}")
print(f"SUMMARY: Entanglement spectrum scaling q={q}")
print(f"{'='*70}")
print(f"{'n':>4} {'chi':>4} {'λ_max':>8} {'w_mult':>8} {'w_tail':>10} {'%S(l0)':>8} {'%S(l1)':>8} {'%S(tail)':>8} {'Δξ':>8} {'S':>8}")
for d in results['data']:
    print(f"{d['n']:4d} {d['chi_used']:4d} {d['lam_max']:8.5f} {d['w_mult']:8.5f} {d['w_tail']:10.6f} "
          f"{d['S_lev0_frac']:8.4f} {d['S_lev1_frac']:8.4f} {d['S_tail_frac']:8.4f} {d['ent_gap']:8.4f} {d['S_total']:8.4f}")

save()
print(f"\nSaved to results/sprint_087b_entspec_dmrg_q7.json")

from db_utils import record
for d in results['data']:
    record(sprint=87, model='sq_potts', q=q, n=d['n'],
           quantity='lam_max', value=d['lam_max'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='entanglement spectrum largest eigenvalue')
    record(sprint=87, model='sq_potts', q=q, n=d['n'],
           quantity='w_tail', value=d['w_tail'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='tail weight (levels >= q)')
    record(sprint=87, model='sq_potts', q=q, n=d['n'],
           quantity='S_tail_frac', value=d['S_tail_frac'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='tail entropy fraction')
    record(sprint=87, model='sq_potts', q=q, n=d['n'],
           quantity='ent_gap', value=d['ent_gap'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='entanglement gap')
print("Recorded to DB.")
