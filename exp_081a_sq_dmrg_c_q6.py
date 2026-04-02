#!/usr/bin/env python3
"""Sprint 081a: DMRG c_eff(q=6) at n=10,12 g_c=1/6 (S_q Potts).

Redesigned: q=6 DMRG is very slow (d=6 local dim). n=8 chi=40 took 172s.
Strategy: chi=40 for n=10 (~250s), chi=30 for n=12 (~200s estimate).
Exact diag already gives c_eff at n=6,7,8 (Sprint 080). DMRG extends to n=10,12.

Key question: Does c_eff stay at ~1.15 or drop toward Re(c)=1.253?
"""
import numpy as np
import json, time

q = 6
gc = 1.0 / q  # exact from self-duality

results = {
    'experiment': '081a_sq_dmrg_c_q6',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'gc': gc,
    'data': [],
    'exact_diag_reference': {
        'n6': {'c_eff': 1.146, 'source': 'Sprint 080'},
        'n7': {'c_eff': 1.148, 'source': 'Sprint 080'},
        'n8': {'c_eff': 1.147, 'source': 'Sprint 080 exact + DMRG verified'},
    },
}

def save():
    with open("results/sprint_081a_sq_dmrg_c_q6.json", "w") as f:
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

def dmrg_entropy(n, q_val, g, chi_max=40):
    model = SqPottsChain({'L': n, 'q': q_val, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n)
    init = [np.random.randint(q_val) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 30,
    })
    E0, _ = eng.run()
    S_profile = [float(psi.entanglement_entropy()[i]) for i in range(n - 1)]
    chi_used = max(psi.chi)
    return float(E0), S_profile, chi_used

def fit_cc(L, S_profile):
    x_vals = np.arange(1, L)
    chord = np.log((2 * L / np.pi) * np.sin(np.pi * x_vals / L))
    S_arr = np.array(S_profile)
    n_pts = len(x_vals)
    lo = n_pts // 4
    hi = 3 * n_pts // 4
    if hi - lo < 2:
        lo = 0; hi = n_pts
    A = np.vstack([chord[lo:hi], np.ones(hi - lo)]).T
    (slope, s0), _, _, _ = np.linalg.lstsq(A, S_arr[lo:hi], rcond=None)
    c = 6 * slope
    S_pred = slope * chord[lo:hi] + s0
    ss_res = np.sum((S_arr[lo:hi] - S_pred) ** 2)
    ss_tot = np.sum((S_arr[lo:hi] - np.mean(S_arr[lo:hi])) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    return float(c), float(R2), float(s0)

# Complex CFT Re(c) for q=6
alpha = np.arccosh(np.sqrt(q) / 2)
Re_c = 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)
print(f"Sprint 081a: S_q Potts DMRG c_eff(q={q}) at g_c=1/{q}={gc:.6f}")
print(f"Complex CFT: Re(c) = {Re_c:.4f}")
print(f"Exact diag (Sprint 080): c_eff(n=6)=1.146, c_eff(n=7)=1.148, c_eff(n=8)=1.147")
print("=" * 60, flush=True)
results['Re_c'] = float(Re_c)

# Run n=10 first (highest priority), then n=12 if time allows
runs = [(10, 40), (12, 30)]

for n, chi in runs:
    print(f"\n--- n={n}, chi_max={chi} ---", flush=True)
    t0 = time.time()
    E0, S, chi_used = dmrg_entropy(n, q, gc, chi_max=chi)
    dt = time.time() - t0

    c_eff, R2, s0 = fit_cc(n, S)
    S_mid = S[n // 2 - 1]
    ratio = c_eff / Re_c

    print(f"  {dt:.1f}s, chi_used={chi_used}, E0={E0:.8f}", flush=True)
    print(f"  S_mid = {S_mid:.4f}", flush=True)
    print(f"  c_eff = {c_eff:.4f}, R² = {R2:.6f}", flush=True)
    print(f"  c_eff/Re(c) = {ratio:.4f}", flush=True)

    entry = {
        'n': n, 'chi_max': chi, 'chi_used': chi_used,
        'E0': E0, 'S_profile': S, 'S_mid': S_mid,
        'c_eff': c_eff, 'R2': R2, 's0': s0,
        'ratio_c_Re_c': ratio, 'time_s': dt,
    }
    results['data'].append(entry)
    save()
    print(f"  Saved.", flush=True)

    if dt > 200:
        print(f"  Took {dt:.0f}s, skipping further sizes.", flush=True)
        break

# Summary with exact diag data included
print(f"\n{'='*60}", flush=True)
print("SUMMARY: c_eff(n) for q=6 S_q Potts (exact diag + DMRG)", flush=True)
print(f"{'='*60}", flush=True)
print(f"{'n':>4} {'method':>10} {'c_eff':>8} {'c/Re(c)':>8}", flush=True)
for n, c in [(6, 1.146), (7, 1.148), (8, 1.147)]:
    print(f"{n:4d} {'exact_diag':>10} {c:8.4f} {c/Re_c:8.4f}", flush=True)
for d in results['data']:
    print(f"{d['n']:4d} {'dmrg':>10} {d['c_eff']:8.4f} {d['ratio_c_Re_c']:8.4f}", flush=True)

# Drift analysis
all_c = [(6, 1.146), (7, 1.148), (8, 1.147)] + [(d['n'], d['c_eff']) for d in results['data']]
all_c.sort()
if len(all_c) >= 2:
    c_first = all_c[0][1]
    c_last = all_c[-1][1]
    drift_pct = 100 * (c_last - c_first) / c_first
    print(f"\nc_eff drift n={all_c[0][0]}→{all_c[-1][0]}: {c_first:.4f} → {c_last:.4f} ({drift_pct:+.2f}%)", flush=True)
    results['c_eff_drift_pct'] = drift_pct
    results['all_c_vs_n'] = all_c

    if abs(drift_pct) < 5:
        print("  → c_eff STABLE: walking extends beyond accessible sizes", flush=True)
        results['walking_status'] = 'extends'
    elif drift_pct < -5:
        print("  → c_eff DROPPING: walking breaking down", flush=True)
        results['walking_status'] = 'breaking'
    else:
        print("  → c_eff RISING: FSS overshoot", flush=True)
        results['walking_status'] = 'overshoot'

save()
print(f"\nFinal save to results/sprint_081a_sq_dmrg_c_q6.json", flush=True)

from db_utils import record
for d in results['data']:
    record(sprint=81, model='sq_potts', q=q, n=d['n'],
           quantity='c_eff', value=d['c_eff'],
           method=f'dmrg_chi{d["chi_used"]}',
           notes=f'gc=1/{q}, CC fit R2={d["R2"]:.4f}')
    record(sprint=81, model='sq_potts', q=q, n=d['n'],
           quantity='S_mid', value=d['S_mid'],
           method=f'dmrg_chi{d["chi_used"]}',
           notes=f'gc=1/{q}, midchain entropy')
print("Recorded to DB.", flush=True)
