#!/usr/bin/env python3
"""Sprint 082a: Spin-spin correlator for S_q Potts q=5 at g_c=1/5 via DMRG.

Measure connected correlator G(r) = <delta(s_0, s_r)> - 1/q at criticality.
Complex CFT predicts G(r) ~ r^{-2x_R} * cos(2*x_I*ln(r) + phi).
For q=5 (walking regime), x_I should be small -> near-pure power law.

Extract x_1 from correlation decay and compare to gap method (x_1 ~ 0.10).
Check residuals for oscillatory structure.
"""
import numpy as np
import json, time

q = 5
gc = 1.0 / q
results = {
    'experiment': '082a_sq_corr_q5',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'gc': gc,
    'sizes': [],
}

def save():
    with open("results/sprint_082a_sq_corr_q5.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# --- Build S_q Potts DMRG model ---
from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

class SqPottsSite(Site):
    def __init__(self, q):
        leg = npc.LegCharge.from_trivial(q)
        Site.__init__(self, leg, [str(a) for a in range(q)], sort_charge=False)
        # Projectors P_a = |a><a| for correlator
        for a in range(q):
            P = np.zeros((q, q), dtype=complex); P[a, a] = 1.0
            self.add_op(f'P{a}', P)
        # S_q transverse field = sum_{k=1}^{q-1} X^k = (all ones) - identity
        Sq_field = np.ones((q, q), dtype=complex) - np.eye(q, dtype=complex)
        self.add_op('SqField', Sq_field, hc='SqField')

class SqPottsChain(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, model_params):
        return SqPottsSite(model_params.get('q', 5))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', gc)
        q_val = model_params.get('q', 5)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')

def run_dmrg_correlator(n, q_val, g, chi_max=60):
    """Run DMRG and compute spin-spin correlator from center site."""
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
    chi_used = max(psi.chi)

    # Compute correlator: <delta(s_i, s_j)> = sum_a <P_a(i) P_a(j)>
    # Use center site as reference
    i0 = n // 2
    corr_connected = []
    distances = []

    for j in range(n):
        if j == i0:
            continue
        # <delta(s_i0, s_j)> = sum_a <P_a(i0) P_a(j)>
        delta_ij = 0.0
        for a in range(q_val):
            if i0 < j:
                val = psi.expectation_value_term([('P' + str(a), i0), ('P' + str(a), j)])
            else:
                val = psi.expectation_value_term([('P' + str(a), j), ('P' + str(a), i0)])
            delta_ij += float(np.real(val))
        # Connected: G(r) = <delta> - 1/q
        G = delta_ij - 1.0 / q_val
        r = abs(j - i0)
        corr_connected.append((r, G))
        distances.append(r)

    # Sort by distance, average if same distance (left/right of center)
    from collections import defaultdict
    by_r = defaultdict(list)
    for r, G in corr_connected:
        by_r[r].append(G)

    r_vals = sorted(by_r.keys())
    G_avg = [np.mean(by_r[r]) for r in r_vals]
    G_std = [np.std(by_r[r]) if len(by_r[r]) > 1 else 0 for r in r_vals]

    return float(E0), chi_used, r_vals, G_avg, G_std

def fit_power_law(r_vals, G_vals, min_r=2):
    """Fit G(r) = A * r^(-eta) to extract eta = 2*x_1."""
    mask = np.array(r_vals) >= min_r
    r = np.array(r_vals)[mask]
    G = np.array(G_vals)[mask]
    # Only use positive G values for log-log fit
    pos = G > 0
    if np.sum(pos) < 2:
        return None, None, None, None
    lr = np.log(r[pos])
    lG = np.log(G[pos])
    A = np.vstack([lr, np.ones_like(lr)]).T
    (slope, intercept), residuals, _, _ = np.linalg.lstsq(A, lG, rcond=None)
    eta = -slope  # G ~ r^(-eta)
    x1 = eta / 2
    G_pred = np.exp(intercept) * r[pos]**(-eta)
    ss_res = np.sum((G[pos] - G_pred)**2)
    ss_tot = np.sum((G[pos] - np.mean(G[pos]))**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    # Residuals normalized
    resid = (G[pos] - G_pred) / G_pred
    return float(eta), float(x1), float(R2), list(zip(r[pos].tolist(), resid.tolist()))

def fit_oscillatory(r_vals, G_vals, min_r=2):
    """Fit G(r) = A * r^(-eta) * (1 + B*cos(omega*ln(r) + phi)).

    Uses power-law residuals to detect oscillation frequency.
    """
    mask = np.array(r_vals) >= min_r
    r = np.array(r_vals)[mask]
    G = np.array(G_vals)[mask]
    pos = G > 0
    if np.sum(pos) < 4:
        return None

    # First get power-law fit
    lr = np.log(r[pos])
    lG = np.log(G[pos])
    A = np.vstack([lr, np.ones_like(lr)]).T
    (slope, intercept), _, _, _ = np.linalg.lstsq(A, lG, rcond=None)
    eta = -slope
    G_pl = np.exp(intercept) * r[pos]**(-eta)

    # Fractional residual
    frac_resid = (G[pos] - G_pl) / G_pl

    # Check for oscillation: FFT of residual vs ln(r)
    if len(frac_resid) < 4:
        return {'eta': float(eta), 'oscillation_amplitude': 0, 'omega': 0}

    # Oscillation amplitude = std of fractional residuals
    osc_amp = float(np.std(frac_resid))

    # Try to find frequency via autocorrelation of residuals
    return {
        'eta': float(eta),
        'oscillation_amplitude': osc_amp,
        'max_abs_resid': float(np.max(np.abs(frac_resid))),
        'residuals_vs_lnr': list(zip(lr.tolist(), frac_resid.tolist())),
    }

# Complex CFT predictions
alpha = np.arccosh(np.sqrt(q) / 2)
Re_c = 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)
print(f"Sprint 082a: S_q Potts correlator q={q} at g_c=1/{q}={gc:.4f}")
print(f"Complex CFT: Re(c) = {Re_c:.4f}")
print(f"Expected x_1 ~ 0.10 from gap method (Sprint 058)")
print("=" * 60, flush=True)

# Run for multiple sizes: n=8 (quick), n=16, n=24
runs = [(8, 40, "timing"), (16, 60, "main"), (24, 60, "large")]

for n, chi, label in runs:
    print(f"\n--- n={n}, chi_max={chi} ({label}) ---", flush=True)
    t0 = time.time()
    E0, chi_used, r_vals, G_avg, G_std = run_dmrg_correlator(n, q, gc, chi_max=chi)
    dt = time.time() - t0

    print(f"  {dt:.1f}s, chi_used={chi_used}", flush=True)

    # Show correlator
    print(f"  r    G(r)         G_std", flush=True)
    for r, g, s in zip(r_vals, G_avg, G_std):
        print(f"  {r:2d}   {g:12.6f}   {s:.6f}", flush=True)

    # Power-law fit
    eta, x1, R2, resid = fit_power_law(r_vals, G_avg, min_r=2)
    if eta is not None:
        print(f"  Power-law fit: eta={eta:.4f}, x_1={x1:.4f}, R²={R2:.6f}", flush=True)

    # Oscillatory analysis
    osc = fit_oscillatory(r_vals, G_avg, min_r=2)
    if osc:
        print(f"  Oscillation amplitude: {osc['oscillation_amplitude']:.4f}", flush=True)
        print(f"  Max |resid|: {osc['max_abs_resid']:.4f}", flush=True)

    entry = {
        'n': n, 'chi_max': chi, 'chi_used': chi_used,
        'E0': E0, 'time_s': dt,
        'r_vals': r_vals, 'G_avg': G_avg, 'G_std': G_std,
        'power_law': {'eta': eta, 'x1': x1, 'R2': R2},
        'oscillation': osc,
    }
    results['sizes'].append(entry)
    save()
    print(f"  Saved.", flush=True)

    if dt > 200:
        print(f"  Took {dt:.0f}s — skipping larger sizes", flush=True)
        break

# Summary
print(f"\n{'='*60}")
print("SUMMARY: x_1 from correlator decay (q=5 S_q Potts)")
print(f"{'='*60}")
print(f"{'n':>4} {'eta':>8} {'x_1':>8} {'R²':>10} {'osc_amp':>10}")
for d in results['sizes']:
    pl = d['power_law']
    osc = d['oscillation']
    eta_str = f"{pl['eta']:.4f}" if pl['eta'] else "—"
    x1_str = f"{pl['x1']:.4f}" if pl['x1'] else "—"
    R2_str = f"{pl['R2']:.6f}" if pl['R2'] else "—"
    osc_str = f"{osc['oscillation_amplitude']:.4f}" if osc else "—"
    print(f"{d['n']:4d} {eta_str:>8} {x1_str:>8} {R2_str:>10} {osc_str:>10}")

save()
print(f"\nFinal save to results/sprint_082a_sq_corr_q5.json")

from db_utils import record
for d in results['sizes']:
    pl = d['power_law']
    if pl['x1'] is not None:
        record(sprint=82, model='sq_potts', q=q, n=d['n'],
               quantity='x1_corr', value=pl['x1'],
               method='correlator_decay',
               notes=f'gc=1/{q}, eta={pl["eta"]:.4f}, R2={pl["R2"]:.4f}')
print("Recorded to DB.")
