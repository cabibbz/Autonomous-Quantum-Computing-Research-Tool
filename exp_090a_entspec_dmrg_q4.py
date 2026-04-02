#!/usr/bin/env python3
"""Sprint 090a: Entanglement spectrum scaling via DMRG at q=4, n=8,12,16,20,24.

q=4 is the real-to-complex CFT boundary: c=1 exactly, (q-1)=3-fold multiplet.
Key prediction: M/[(q-1)/q] ≈ 1.0 at q=4 (democracy index crossover).

Uses open-BC finite DMRG. Schmidt spectrum at midchain bond.
"""
import numpy as np
import json, time

q = 4
gc = 1.0 / q

results = {
    'experiment': '090a_entspec_dmrg_q4',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'gc': gc,
    'data': [],
}

def save():
    with open("results/sprint_090a_entspec_dmrg_q4.json", "w") as f:
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
    """Analyze entanglement spectrum: lambda_max, multiplet, tail."""
    spec = np.sort(spec)[::-1]  # descending
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

# q=4 is the boundary case: c = 1 exactly (real, Ashkin-Teller universality)
Re_c = 1.0
results['Re_c'] = Re_c

print(f"Sprint 090a: DMRG entanglement spectrum scaling for S_q Potts q={q}")
print(f"g_c = 1/{q} = {gc:.6f}, c = {Re_c} (exact)")
print("=" * 70, flush=True)

# q=4, d=4 per site. Intermediate speed.
# Start with a timing test at n=8
runs = [(8, 50), (12, 70), (16, 90), (20, 110), (24, 130)]

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
        'max_sweeps': 40,
    })
    E0, _ = eng.run()
    dt_dmrg = time.time() - t0
    chi_used = max(psi.chi)

    print(f"  DMRG: {dt_dmrg:.1f}s, E0={E0:.8f}, chi_used={chi_used}", flush=True)

    # Check timing — abort if too slow
    if dt_dmrg > 250 and n < 24:
        print(f"  WARNING: {dt_dmrg:.0f}s — may need to skip larger sizes")

    # Get Schmidt spectrum at midchain bond
    mid_bond = n // 2 - 1
    sv_mid = psi.get_SL(mid_bond)
    spec_mid = sv_mid**2

    analysis = analyze_spectrum(spec_mid, q)

    # Compute M and democracy index
    M = analysis['S_lev1_frac'] / (analysis['S_lev0_frac'] + analysis['S_lev1_frac']) \
        if (analysis['S_lev0_frac'] + analysis['S_lev1_frac']) > 0 else 0
    qm1_q = (q - 1) / q
    democracy = M / qm1_q if qm1_q > 0 else 0

    print(f"  Schmidt rank: {analysis['schmidt_rank']}, ent_gap: {analysis['ent_gap']:.4f}")
    print(f"  lam_max: {analysis['lam_max']:.6f}, w_mult: {analysis['w_mult']:.6f}, w_tail: {analysis['w_tail']:.6f}")
    print(f"  S = {analysis['S_total']:.6f}")
    print(f"  Entropy fracs: lev0={analysis['S_lev0_frac']:.4f}, lev1={analysis['S_lev1_frac']:.4f}, tail={analysis['S_tail_frac']:.4f}")
    print(f"  M = {M:.4f}, (q-1)/q = {qm1_q:.4f}, M/[(q-1)/q] = {democracy:.4f}")
    print(f"  Multiplet spread: {analysis['multiplet_spread']:.6f}")

    entry = {
        'n': n, 'chi_max': chi_max, 'chi_used': chi_used,
        'E0': float(E0), 'time_dmrg': dt_dmrg,
        'M': float(M), 'democracy': float(democracy),
        'spec_mid_top20': [float(x) for x in spec_mid[:20]],
        **analysis,
    }
    results['data'].append(entry)
    save()
    print(f"  Saved.", flush=True)

# Summary table
print(f"\n{'='*70}")
print(f"SUMMARY: Entanglement spectrum scaling q={q}")
print(f"{'='*70}")
print(f"{'n':>4} {'chi':>4} {'lam_max':>8} {'w_mult':>8} {'w_tail':>10} {'%S(l0)':>8} {'%S(l1)':>8} {'%S(tail)':>8} {'M':>8} {'M/ref':>8} {'Dxi':>8} {'S':>8} {'t(s)':>6}")
for d in results['data']:
    print(f"{d['n']:4d} {d['chi_used']:4d} {d['lam_max']:8.5f} {d['w_mult']:8.5f} {d['w_tail']:10.6f} "
          f"{d['S_lev0_frac']:8.4f} {d['S_lev1_frac']:8.4f} {d['S_tail_frac']:8.4f} {d['M']:8.4f} {d['democracy']:8.4f} "
          f"{d['ent_gap']:8.4f} {d['S_total']:8.4f} {d['time_dmrg']:6.1f}")

# c_eff from size pairs
print(f"\nc_eff from size pairs (open BC: S = c/6 * ln(n) + const):")
dd = results['data']
for i in range(1, len(dd)):
    n1, S1 = dd[i-1]['n'], dd[i-1]['S_total']
    n2, S2 = dd[i]['n'], dd[i]['S_total']
    c_eff = 6 * (S2 - S1) / (np.log(n2) - np.log(n1))
    ratio = c_eff / Re_c
    print(f"  ({n1},{n2}): c_eff = {c_eff:.4f}, c_eff/Re(c) = {ratio:.4f}")

save()
print(f"\nSaved to results/sprint_090a_entspec_dmrg_q4.json")

from db_utils import record
for d in results['data']:
    record(sprint=90, model='sq_potts', q=q, n=d['n'],
           quantity='lam_max', value=d['lam_max'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='entanglement spectrum largest eigenvalue')
    record(sprint=90, model='sq_potts', q=q, n=d['n'],
           quantity='w_tail', value=d['w_tail'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='tail weight (levels >= q)')
    record(sprint=90, model='sq_potts', q=q, n=d['n'],
           quantity='S_tail_frac', value=d['S_tail_frac'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='tail entropy fraction')
    record(sprint=90, model='sq_potts', q=q, n=d['n'],
           quantity='ent_gap', value=d['ent_gap'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='entanglement gap')
    record(sprint=90, model='sq_potts', q=q, n=d['n'],
           quantity='M_dominance', value=d['M'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='multiplet dominance S_lev1/(S_lev0+S_lev1)')
    record(sprint=90, model='sq_potts', q=q, n=d['n'],
           quantity='democracy_index', value=d['democracy'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='M/[(q-1)/q] democracy index')
print("Recorded to DB.")
