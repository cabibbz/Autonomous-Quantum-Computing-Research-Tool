#!/usr/bin/env python3
"""Sprint 086a: DMRG Rényi entropies at q=5, n=8,10,12 (open BC, S_q Potts).

Extract entanglement spectrum from DMRG MPS at each bond, compute S_α for
multiple α, fit Calabrese-Cardy to get c_α. Compare to exact diag (Sprint 085).

Key question: Does α=1 remain the most accurate Rényi index for walking q=5
at larger system sizes?
"""
import numpy as np
import json, time

q = 5
gc = 1.0 / q

results = {
    'experiment': '086a_dmrg_renyi_q5',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'gc': gc,
    'data': [],
}

def save():
    with open("results/sprint_086a_dmrg_renyi_q5.json", "w") as f:
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
    """Compute Rényi entropy S_α from Schmidt values λ_i.
    S_α = (1/(1-α)) ln(Σ λ_i^{2α}) for α≠1.
    S_1 = -Σ λ_i² ln(λ_i²) (von Neumann).
    """
    p = schmidt_vals**2  # probabilities
    p = p[p > 1e-30]
    if alpha == 1:
        return float(-np.sum(p * np.log(p)))
    elif alpha == np.inf:
        return float(-np.log(np.max(p)))
    else:
        return float(np.log(np.sum(p**alpha)) / (1.0 - alpha))

def fit_cc_renyi(L, S_profile_alpha, alpha):
    """Fit Calabrese-Cardy for open BC: S_α(ℓ) = (c_α/6)(1+1/α) ln[(2L/π)sin(πℓ/L)] + const.
    Returns c_α, R², const."""
    x_vals = np.arange(1, L)  # bond positions ℓ = 1..L-1
    chord = np.log((2 * L / np.pi) * np.sin(np.pi * x_vals / L))
    S_arr = np.array(S_profile_alpha)

    # Use central half for fit (avoid boundary effects)
    n_pts = len(x_vals)
    lo = n_pts // 4
    hi = 3 * n_pts // 4
    if hi - lo < 3:
        lo = max(0, n_pts // 2 - 2)
        hi = min(n_pts, n_pts // 2 + 2)

    A = np.vstack([chord[lo:hi], np.ones(hi - lo)]).T
    (slope, s0), _, _, _ = np.linalg.lstsq(A, S_arr[lo:hi], rcond=None)

    if alpha == np.inf:
        prefactor = 1.0 / 6.0  # (1+1/∞)/6 = 1/6
    else:
        prefactor = (1.0 + 1.0/alpha) / 6.0
    c_alpha = slope / prefactor

    S_pred = slope * chord[lo:hi] + s0
    ss_res = np.sum((S_arr[lo:hi] - S_pred)**2)
    ss_tot = np.sum((S_arr[lo:hi] - np.mean(S_arr[lo:hi]))**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 1e-30 else 0

    return float(c_alpha), float(R2), float(s0)

# Complex CFT Re(c) for q=5
alpha_cft = np.arccosh(np.sqrt(q) / 2)
Re_c = 1 + 6 * alpha_cft**2 / (np.pi**2 + alpha_cft**2)
results['Re_c'] = float(Re_c)

alphas = [0.5, 1, 2, 3, 5, 10, np.inf]
alpha_labels = ['0.5', '1', '2', '3', '5', '10', 'inf']

print(f"Sprint 086a: DMRG Rényi entropies for S_q Potts q={q}")
print(f"g_c = 1/{q} = {gc:.6f}, Re(c) = {Re_c:.4f}")
print("=" * 70, flush=True)

runs = [(8, 60), (10, 80), (12, 60)]

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

    # Get entanglement spectrum at each bond and compute Rényi entropies
    S_profiles = {label: [] for label in alpha_labels}

    for bond in range(n - 1):
        sv = psi.get_SL(bond)  # Schmidt values at bond
        for a, label in zip(alphas, alpha_labels):
            S_profiles[label].append(renyi_from_schmidt(sv, a))

    # Fit CC for each α
    c_alpha = {}
    R2_alpha = {}
    for a, label in zip(alphas, alpha_labels):
        c_a, R2, _ = fit_cc_renyi(n, S_profiles[label], a)
        c_alpha[label] = c_a
        R2_alpha[label] = R2

    # Print results
    print(f"\n  {'α':>5} {'c_α':>8} {'c_α/Re(c)':>10} {'R²':>10} {'S_mid':>8}")
    best_label = None
    best_dev = 999
    for label in alpha_labels:
        ratio = c_alpha[label] / Re_c
        S_mid = S_profiles[label][n // 2 - 1]
        dev = abs(ratio - 1.0)
        if dev < best_dev:
            best_dev = dev
            best_label = label
        marker = " *" if label == best_label else ""
        print(f"  {label:>5} {c_alpha[label]:8.4f} {ratio:10.4f} {R2_alpha[label]:10.6f} {S_mid:8.4f}{marker}")

    print(f"\n  Best α = {best_label} (|c_α/Re(c) - 1| = {best_dev:.4f})")

    entry = {
        'n': n, 'chi_max': chi_max, 'chi_used': chi_used,
        'E0': float(E0), 'time_dmrg': dt_dmrg,
        'c_alpha': c_alpha, 'R2_alpha': R2_alpha,
        'c_alpha_over_Rec': {label: c_alpha[label] / Re_c for label in alpha_labels},
        'S_mid': {label: S_profiles[label][n // 2 - 1] for label in alpha_labels},
        'best_alpha': best_label,
        'best_dev': best_dev,
    }
    results['data'].append(entry)
    save()
    print(f"  Saved.", flush=True)

# Summary
print(f"\n{'='*70}")
print("SUMMARY: c_α / Re(c) for q=5 (DMRG, open BC)")
print(f"{'='*70}")
header = f"{'n':>4} {'chi':>4}"
for label in alpha_labels:
    header += f" {'α='+label:>8}"
header += "  best"
print(header)
for d in results['data']:
    row = f"{d['n']:4d} {d['chi_used']:4d}"
    for label in alpha_labels:
        row += f" {d['c_alpha_over_Rec'][label]:8.4f}"
    row += f"  α={d['best_alpha']}"
    print(row)

# Compare exact diag (periodic, Sprint 085) vs DMRG (open)
print(f"\nExact diag (periodic, Sprint 085) at n=8:")
print(f"  c₁/Re(c)=0.999, c₂/Re(c)=1.126, c₃/Re(c)=1.132, best α=1")

save()
print(f"\nSaved to results/sprint_086a_dmrg_renyi_q5.json")

from db_utils import record
for d in results['data']:
    for label in alpha_labels:
        record(sprint=86, model='sq_potts', q=q, n=d['n'],
               quantity=f'c_alpha_{label}_dmrg', value=d['c_alpha'][label],
               method=f'dmrg_chi{d["chi_used"]}_open',
               notes=f'c/Re(c)={d["c_alpha_over_Rec"][label]:.4f}, R2={d["R2_alpha"][label]:.4f}')
print("Recorded to DB.")
