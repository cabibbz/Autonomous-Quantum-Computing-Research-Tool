#!/usr/bin/env python3
"""Sprint 079b: Central charge c for S_q Potts at q=10 via DMRG entropy profile.

Use Calabrese-Cardy formula at exact g_c = 1/10 for open BC chain.
q=10 has d=10 per site so chi needs to be larger for convergence.
Time each size and stop if approaching 200s.
"""
import numpy as np
import json, time

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

results = {
    'experiment': '079b_sq_potts_c_q10',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': 10, 'gc': 0.1,
}

def save():
    with open("results/sprint_079b_sq_potts_c_q10.json", "w") as f:
        json.dump(results, f, indent=2, default=str)


def complex_cft_c(q):
    """Compute complex CFT central charge for q-state Potts model."""
    if q <= 4:
        g = np.arccos(np.sqrt(q) / 2) / np.pi
        c = 1 - 6 * (1 - g)**2 / g
        return complex(c, 0), complex(g, 0)
    else:
        theta = np.arccosh((q - 2) / 2) / (2 * np.pi)
        g = 0.5 + 1j * theta
        c = 1 - 6 * (1 - g)**2 / g
        return c, g


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
        return SqPottsSite(model_params.get('q', 10))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 0.1)
        q_val = model_params.get('q', 10)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')


def dmrg_entropy_profile(n, q, g, chi_max=80):
    """Get ground state entropy profile via DMRG."""
    model = SqPottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n + q * 100)
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 40,
    })
    E0, _ = eng.run()
    S = [float(psi.entanglement_entropy()[i]) for i in range(n - 1)]
    chi_actual = max(psi.chi)
    return float(E0), S, chi_actual


def fit_cc(L, S_profile):
    """Fit Calabrese-Cardy formula, central half only."""
    x_vals = np.arange(1, L)
    chord = np.log((2 * L / np.pi) * np.sin(np.pi * x_vals / L))
    S_arr = np.array(S_profile)

    n_pts = len(x_vals)
    lo = n_pts // 4
    hi = 3 * n_pts // 4
    chord_c = chord[lo:hi]
    S_c = S_arr[lo:hi]

    A = np.vstack([chord_c, np.ones_like(chord_c)]).T
    (slope, s0), residuals, _, _ = np.linalg.lstsq(A, S_c, rcond=None)
    c = 6 * slope

    S_pred = slope * chord_c + s0
    ss_res = np.sum((S_c - S_pred) ** 2)
    ss_tot = np.sum((S_c - np.mean(S_c)) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    # Full range
    A_full = np.vstack([chord, np.ones_like(chord)]).T
    (slope_f, s0_f), _, _, _ = np.linalg.lstsq(A_full, S_arr, rcond=None)
    c_full = 6 * slope_f
    S_pred_f = slope_f * chord + s0_f
    ss_res_f = np.sum((S_arr - S_pred_f) ** 2)
    ss_tot_f = np.sum((S_arr - np.mean(S_arr)) ** 2)
    R2_full = 1 - ss_res_f / ss_tot_f if ss_tot_f > 0 else 0

    return float(c), float(s0), float(R2), float(c_full), float(R2_full)


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 079b: Central charge c for S_q Potts at q=10", flush=True)
print("=" * 60, flush=True)

# Complex CFT prediction
c_cft, g_cft = complex_cft_c(10)
print(f"\nComplex CFT prediction for q=10:", flush=True)
print(f"  g = {g_cft:.6f}", flush=True)
print(f"  c = {c_cft:.6f}", flush=True)
print(f"  Re(c) = {c_cft.real:.6f}, Im(c) = {c_cft.imag:.6f}", flush=True)
results['complex_cft'] = {'Re_c': c_cft.real, 'Im_c': c_cft.imag,
                           'Re_g': g_cft.real, 'Im_g': g_cft.imag}
save()

# DMRG entropy profile for q=10
q = 10
gc = 0.1
print(f"\n--- q=10, g_c = 1/10 = {gc} ---", flush=True)

# q=10 has d=10 per site — needs higher chi but slower convergence
# Start small to test timing
q10_data = {}
for n in [8, 12, 16, 20]:
    chi = min(40 + n * 3, 120)
    print(f"\n  n={n}, chi_max={chi}...", flush=True)
    t0 = time.time()
    E0, S, chi_actual = dmrg_entropy_profile(n, q, gc, chi_max=chi)
    dt = time.time() - t0
    c_cen, s0, R2, c_full, R2_full = fit_cc(n, S)
    S_mid = S[n // 2 - 1]
    print(f"  {dt:.1f}s, chi_actual={chi_actual}, E0={E0:.8f}", flush=True)
    print(f"  c(central)={c_cen:.4f} R²={R2:.6f}, c(full)={c_full:.4f} R²={R2_full:.6f}", flush=True)
    print(f"  S_mid={S_mid:.4f}", flush=True)

    q10_data[f'n{n}'] = {
        'n': n, 'chi_max': chi, 'chi_actual': chi_actual,
        'c_central': c_cen, 'c_full': c_full,
        'R2_central': R2, 'R2_full': R2_full,
        's0': s0, 'S_mid': S_mid, 'S_profile': S,
        'E0': E0, 'time_s': dt
    }
    results['q10'] = q10_data
    save()

    if dt > 200:
        print(f"  Approaching time limit, stopping.", flush=True)
        break

# Summary
print(f"\n{'='*60}", flush=True)
print("SUMMARY: c(q=10) from DMRG entropy profile", flush=True)
print(f"{'='*60}", flush=True)
print(f"Complex CFT prediction: Re(c) = {c_cft.real:.4f}", flush=True)
print(f"\n{'n':>4} {'chi':>4} {'c(central)':>12} {'R²':>10} {'c(full)':>10} {'time':>8}", flush=True)
for k, v in q10_data.items():
    print(f"{v['n']:4d} {v['chi_actual']:4d} {v['c_central']:12.4f} {v['R2_central']:10.6f} {v['c_full']:10.4f} {v['time_s']:8.1f}s", flush=True)

c_vals = [v['c_central'] for v in q10_data.values()]
if len(c_vals) >= 2:
    print(f"\nc convergence: {' -> '.join(f'{c:.4f}' for c in c_vals)}", flush=True)
    print(f"Last value: {c_vals[-1]:.4f}, CFT: {c_cft.real:.4f}, ratio: {c_vals[-1]/c_cft.real:.4f}", flush=True)

save()
print("\nSaved to results/sprint_079b_sq_potts_c_q10.json", flush=True)

from db_utils import record
for k, v in q10_data.items():
    record(sprint=79, model='sq_potts', q=10, n=v['n'],
           quantity='c_eff', value=v['c_central'],
           method='dmrg_CC_central', notes=f'gc=1/10, chi={v["chi_actual"]}, R2={v["R2_central"]:.4f}')
print("Recorded to DB.", flush=True)
