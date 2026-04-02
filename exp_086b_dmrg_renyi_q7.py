#!/usr/bin/env python3
"""Sprint 086b: DMRG Rényi entropies at q=7, n=6,8,10,12 (open BC, S_q Potts).

Extract Schmidt spectrum from DMRG MPS at midchain bond, compute S_α for
multiple α. Focus on midchain entropy for size-pair c_α extraction.

Key question: Does α=3 remain the most accurate Rényi index for walking-broken q=7
at larger system sizes?
"""
import numpy as np
import json, time

q = 7
gc = 1.0 / q

results = {
    'experiment': '086b_dmrg_renyi_q7',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'gc': gc,
    'data': [],
}

def save():
    with open("results/sprint_086b_dmrg_renyi_q7.json", "w") as f:
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

def renyi_from_schmidt(schmidt_vals, alpha):
    """Compute Rényi entropy S_α from Schmidt values."""
    p = schmidt_vals**2
    p = p[p > 1e-30]
    if alpha == 1:
        return float(-np.sum(p * np.log(p)))
    elif alpha == np.inf:
        return float(-np.log(np.max(p)))
    else:
        return float(np.log(np.sum(p**alpha)) / (1.0 - alpha))

# Complex CFT Re(c) for q=7
alpha_cft = np.arccosh(np.sqrt(q) / 2)
Re_c = 1 + 6 * alpha_cft**2 / (np.pi**2 + alpha_cft**2)
results['Re_c'] = float(Re_c)

alphas = [0.5, 1, 2, 3, 5, 10, np.inf]
alpha_labels = ['0.5', '1', '2', '3', '5', '10', 'inf']

print(f"Sprint 086b: DMRG Rényi entropies for S_q Potts q={q}")
print(f"g_c = 1/{q} = {gc:.6f}, Re(c) = {Re_c:.4f}")
print("=" * 70, flush=True)

# q=7: d=7 per site. n=8 chi=56 took 172s (Sprint 081 at q=6 d=6).
# Be conservative: n=6 chi=40, n=8 chi=50, n=10 chi=40, n=12 chi=30
runs = [(6, 40), (8, 50), (10, 40), (12, 30)]

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
    chi_used = max(psi.chi)

    print(f"  DMRG: {dt_dmrg:.1f}s, E0={E0:.8f}, chi_used={chi_used}", flush=True)

    if dt_dmrg > 250:
        print(f"  WARNING: Approaching time limit, will skip further sizes.", flush=True)

    # Get Schmidt values at midchain bond and ALL bonds
    mid_bond = n // 2 - 1
    S_mid = {}
    S_profile = {label: [] for label in alpha_labels}

    for bond in range(n - 1):
        sv = psi.get_SL(bond)
        for a, label in zip(alphas, alpha_labels):
            S_val = renyi_from_schmidt(sv, a)
            S_profile[label].append(S_val)
            if bond == mid_bond:
                S_mid[label] = S_val

    # Also get entanglement spectrum at midchain
    sv_mid = psi.get_SL(mid_bond)
    spec_mid = sv_mid**2  # eigenvalues of rho_A
    xi_mid = -np.log(np.clip(spec_mid, 1e-30, None))

    print(f"\n  Midchain bond={mid_bond}, Schmidt rank={len(sv_mid)}")
    print(f"  Top eigenvalues: {[f'{x:.5f}' for x in spec_mid[:8]]}")
    print(f"\n  {'α':>5} {'S_mid':>10}")
    for label in alpha_labels:
        print(f"  {label:>5} {S_mid[label]:10.6f}")

    entry = {
        'n': n, 'chi_max': chi_max, 'chi_used': chi_used,
        'E0': float(E0), 'time_dmrg': dt_dmrg,
        'S_mid': S_mid,
        'S_profile': S_profile,
        'spec_mid_top': [float(x) for x in spec_mid[:20]],
        'xi_mid_top': [float(x) for x in xi_mid[:20]],
    }
    results['data'].append(entry)
    save()
    print(f"  Saved.", flush=True)

    if dt_dmrg > 250:
        print(f"  Stopping: time limit.", flush=True)
        break

# Size-pair c_α extraction using midchain entropy
# Open BC formula: S_α(L/2) = (c/12)(1+1/α) ln(2L/π) + const
# Size pairs: c_α = 12·ΔS / ((1+1/α)·Δln(2L/π)) = 12·ΔS / ((1+1/α)·ln(L₂/L₁))
print(f"\n{'='*70}")
print("SIZE-PAIR c_α EXTRACTION (open BC)")
print(f"{'='*70}")

sizes_done = [d['n'] for d in results['data']]
pairs = []
for i in range(len(sizes_done)):
    for j in range(i+1, len(sizes_done)):
        pairs.append((i, j))

results['size_pairs'] = []
for i, j in pairs:
    n1, n2 = sizes_done[i], sizes_done[j]
    d1, d2 = results['data'][i], results['data'][j]
    print(f"\n  Pair ({n1}, {n2}):")
    delta_ln = np.log(n2 / n1)

    c_alpha_pair = {}
    for a, label in zip(alphas, alpha_labels):
        a_num = np.inf if a == np.inf else a
        dS = d2['S_mid'][label] - d1['S_mid'][label]
        if a_num == np.inf:
            prefactor = 1.0
        else:
            prefactor = 1.0 + 1.0/a_num
        c_a = 12.0 * dS / (prefactor * delta_ln)
        c_alpha_pair[label] = c_a

    print(f"  {'α':>5} {'c_α':>8} {'c_α/Re(c)':>10}")
    best_label = None
    best_dev = 999
    for label in alpha_labels:
        ratio = c_alpha_pair[label] / Re_c
        dev = abs(ratio - 1.0)
        if dev < best_dev:
            best_dev = dev
            best_label = label
        print(f"  {label:>5} {c_alpha_pair[label]:8.4f} {ratio:10.4f}")

    print(f"  Best α = {best_label} (|c_α/Re(c) - 1| = {best_dev:.4f})")

    results['size_pairs'].append({
        'n1': n1, 'n2': n2,
        'c_alpha': c_alpha_pair,
        'c_alpha_over_Rec': {label: c_alpha_pair[label] / Re_c for label in alpha_labels},
        'best_alpha': best_label,
        'best_dev': best_dev,
    })
    save()

# Summary
print(f"\n{'='*70}")
print("SUMMARY: c_α / Re(c) from size pairs (q=7, open BC DMRG)")
print(f"{'='*70}")
header = f"{'pair':>10}"
for label in alpha_labels:
    header += f" {'α='+label:>8}"
header += "  best"
print(header)
for sp in results['size_pairs']:
    row = f"({sp['n1']:2d},{sp['n2']:2d})  "
    for label in alpha_labels:
        row += f" {sp['c_alpha_over_Rec'][label]:8.4f}"
    row += f"  α={sp['best_alpha']}"
    print(row)

save()
print(f"\nSaved to results/sprint_086b_dmrg_renyi_q7.json")

from db_utils import record
for d in results['data']:
    for label in alpha_labels:
        record(sprint=86, model='sq_potts', q=q, n=d['n'],
               quantity=f'S_mid_alpha_{label}_dmrg', value=d['S_mid'][label],
               method=f'dmrg_chi{d["chi_used"]}_open',
               notes=f'midchain Renyi entropy')
# Record size-pair c values for the widest pair
if results['size_pairs']:
    sp = results['size_pairs'][-1]  # widest pair
    for label in alpha_labels:
        record(sprint=86, model='sq_potts', q=q, n=sp['n2'],
               quantity=f'c_alpha_{label}_pair', value=sp['c_alpha'][label],
               method=f'dmrg_size_pair_{sp["n1"]}_{sp["n2"]}',
               notes=f'c/Re(c)={sp["c_alpha_over_Rec"][label]:.4f}')
print("Recorded to DB.")
