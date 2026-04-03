#!/usr/bin/env python3
"""Sprint 089a: Extend q=7 entanglement spectrum to n=14,16 via DMRG.

Key question: is the b≈3.0 tail exponent from n=6-12 (4 pts) genuine,
or pre-asymptotic (converging to universal b≈2.0)?

Strategy: DMRG at g_c=1/7, open BC, midchain Schmidt spectrum.
n=12 took 293s with chi=40. Try n=14 with chi=30 first (timing test).
"""
import numpy as np
import json, time

q = 7
gc = 1.0 / q

results = {
    'experiment': '089a_entspec_dmrg_q7_ext',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'gc': gc,
    'data': [],
}

def save():
    with open("results/sprint_089a_entspec_dmrg_q7_ext.json", "w") as f:
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

# Complex CFT Re(c)
alpha_cft = np.arccosh(np.sqrt(q) / 2)
Re_c = 1 + 6 * alpha_cft**2 / (np.pi**2 + alpha_cft**2)
results['Re_c'] = float(Re_c)

print(f"Sprint 089a: Extend q={q} DMRG entanglement spectrum to n=14,16")
print(f"g_c = 1/{q} = {gc:.6f}, Re(c) = {Re_c:.4f}")
print("=" * 70, flush=True)

# Conservative approach: start with n=14 chi=30, monitor timing
# n=12 chi=40 took 293s. Reducing chi to 30 and using fewer sweeps.
runs = [(14, 30), (16, 25)]

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
        'mixer': True, 'max_E_err': 1e-9,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-11},
        'max_sweeps': 20,
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

    # Stop if running too long
    if total_time > 400:
        print(f"  Stopping: cumulative time {total_time:.0f}s > 400s.", flush=True)
        break

# Also load previous q=7 data for combined analysis
print(f"\n{'='*70}")
print("Loading previous q=7 data (n=6,8,10,12) from Sprint 087b...")
try:
    with open("results/sprint_087b_entspec_dmrg_q7.json") as f:
        prev = json.load(f)
    prev_data = prev['data']
    all_ns = [d['n'] for d in prev_data] + [d['n'] for d in results['data']]
    all_wtails = [d['w_tail'] for d in prev_data] + [d['w_tail'] for d in results['data']]

    print(f"\nCombined q=7 data (n={min(all_ns)}-{max(all_ns)}):")
    print(f"{'n':>4} {'w_tail':>12} {'ln(n)':>8} {'ln(w_tail)':>12}")
    for n_val, wt in zip(all_ns, all_wtails):
        print(f"{n_val:4d} {wt:12.6f} {np.log(n_val):8.3f} {np.log(max(wt,1e-30)):12.3f}")

    # Power-law fit: log-log
    from scipy.optimize import curve_fit
    def log_fit(x, a, b):
        return a + b * x

    mask = np.array(all_wtails) > 0
    ln_n = np.log(np.array(all_ns)[mask])
    ln_wt = np.log(np.array(all_wtails)[mask])

    popt, pcov = curve_fit(log_fit, ln_n, ln_wt)
    ln_A, b = popt
    A = np.exp(ln_A)
    b_err = np.sqrt(pcov[1,1])

    predicted = log_fit(ln_n, *popt)
    ss_res = np.sum((ln_wt - predicted)**2)
    ss_tot = np.sum((ln_wt - np.mean(ln_wt))**2)
    R2 = 1 - ss_res / ss_tot

    print(f"\nPower-law fit (all {len(ln_n)} points): w_tail = {A:.2e} × n^{b:.3f}")
    print(f"  b = {b:.3f} ± {b_err:.3f}, R² = {R2:.4f}")
    print(f"  Previous (n=6-12, 4pts): b ≈ 2.97")
    print(f"  Change: {b:.3f} vs 2.97 → {'converging toward 2.0' if b < 2.97 else 'stable at ~3'}")

    # Also fit just n=8-max (dropping n=6 which may be in crossover)
    mask2 = np.array(all_ns) >= 8
    if np.sum(mask2) >= 3:
        ln_n2 = np.log(np.array(all_ns)[mask2])
        ln_wt2 = np.log(np.array(all_wtails)[mask2])
        popt2, pcov2 = curve_fit(log_fit, ln_n2, ln_wt2)
        b2 = popt2[1]
        b2_err = np.sqrt(pcov2[1,1])
        print(f"  Fit excluding n=6 (n≥8): b = {b2:.3f} ± {b2_err:.3f}")

    results['combined_fit'] = {
        'A': float(A), 'b': float(b), 'b_err': float(b_err), 'R2': float(R2),
        'n_points': int(len(ln_n)),
        'ns': [int(x) for x in all_ns],
        'w_tails': [float(x) for x in all_wtails],
    }
    save()

except Exception as e:
    print(f"  Could not load previous data: {e}")

# Summary
print(f"\n{'='*70}")
print("SUMMARY 089a")
print(f"{'='*70}")
print(f"{'n':>4} {'chi':>4} {'λ_max':>8} {'w_mult':>8} {'w_tail':>10} {'%S(l0)':>8} {'%S(l1)':>8} {'%S(tail)':>8} {'Δξ':>8}")
for d in results['data']:
    print(f"{d['n']:4d} {d['chi_used']:4d} {d['lam_max']:8.5f} {d['w_mult']:8.5f} {d['w_tail']:10.6f} "
          f"{d['S_lev0_frac']:8.4f} {d['S_lev1_frac']:8.4f} {d['S_tail_frac']:8.4f} {d['ent_gap']:8.4f}")

save()
print(f"\nSaved to results/sprint_089a_entspec_dmrg_q7_ext.json")

from db_utils import record
for d in results['data']:
    record(sprint=89, model='sq_potts', q=q, n=d['n'],
           quantity='w_tail', value=d['w_tail'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='tail weight (levels >= q)')
    record(sprint=89, model='sq_potts', q=q, n=d['n'],
           quantity='S_tail_frac', value=d['S_tail_frac'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='tail entropy fraction')
    record(sprint=89, model='sq_potts', q=q, n=d['n'],
           quantity='lam_max', value=d['lam_max'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='largest Schmidt eigenvalue')
    record(sprint=89, model='sq_potts', q=q, n=d['n'],
           quantity='ent_gap', value=d['ent_gap'],
           method=f'dmrg_chi{d["chi_used"]}_open',
           notes='entanglement gap')
print("Recorded to DB.")
