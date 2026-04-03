#!/usr/bin/env python3
"""Sprint 079a: Central charge c for S_q Potts at q=7 via DMRG entropy profile.

Use Calabrese-Cardy formula S(x) = (c/6)*ln[(2L/pi)*sin(pi*x/L)] + const
at exact g_c = 1/7 for open BC chain.

Also compute complex CFT prediction: Q = 4*cos^2(pi*g), c = 1 - 6*(1-g)^2/g.
For q>4, g is complex: g = 1/2 + i*theta where cosh(2*pi*theta) = (q-2)/2.
"""
import numpy as np
import json, time

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

results = {
    'experiment': '079a_sq_potts_c_q7',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': 7, 'gc': 1.0/7,
}

def save():
    with open("results/sprint_079a_sq_potts_c_q7.json", "w") as f:
        json.dump(results, f, indent=2, default=str)


def complex_cft_c(q):
    """Compute complex CFT central charge for q-state Potts model.

    Coulomb gas: Q = 4*cos^2(pi*g), c = 1 - 6*(1-g)^2/g
    For q>4: g = 1/2 + i*theta where cosh(2*pi*theta) = (q-2)/2
    """
    if q <= 4:
        # Real g
        g = np.arccos(np.sqrt(q) / 2) / np.pi
        c = 1 - 6 * (1 - g)**2 / g
        return complex(c, 0), complex(g, 0)
    else:
        # Complex g: cosh(2*pi*theta) = (q-2)/2
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
        return SqPottsSite(model_params.get('q', 7))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 1.0/7)
        q_val = model_params.get('q', 7)
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
    """Fit Calabrese-Cardy formula to extract c.
    S(x) = (c/6)*ln[(2L/pi)*sin(pi*x/L)] + s0
    Use central half only (avoid boundary effects).
    """
    x_vals = np.arange(1, L)
    chord = np.log((2 * L / np.pi) * np.sin(np.pi * x_vals / L))
    S_arr = np.array(S_profile)

    # Use central half to avoid boundary effects
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

    # Also fit full range for comparison
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
print("Sprint 079a: Central charge c for S_q Potts at q=7", flush=True)
print("=" * 60, flush=True)

# Complex CFT prediction
c_cft, g_cft = complex_cft_c(7)
print(f"\nComplex CFT prediction for q=7:", flush=True)
print(f"  g = {g_cft:.6f}", flush=True)
print(f"  c = {c_cft:.6f}", flush=True)
print(f"  Re(c) = {c_cft.real:.6f}, Im(c) = {c_cft.imag:.6f}", flush=True)
results['complex_cft'] = {'Re_c': c_cft.real, 'Im_c': c_cft.imag,
                           'Re_g': g_cft.real, 'Im_g': g_cft.imag}
save()

# Also compute for q=5 as cross-check
c5, g5 = complex_cft_c(5)
print(f"\nCross-check q=5: Re(c) = {c5.real:.6f} (known: ~1.138)", flush=True)

# DMRG entropy profile for q=7
q = 7
gc = 1.0 / 7
print(f"\n--- q=7, g_c = 1/7 = {gc:.6f} ---", flush=True)

q7_data = {}
for n in [8, 12, 16, 20, 24]:
    chi = min(40 + n * 2, 120)
    print(f"\n  n={n}, chi_max={chi}...", flush=True)
    t0 = time.time()
    E0, S, chi_actual = dmrg_entropy_profile(n, q, gc, chi_max=chi)
    dt = time.time() - t0
    c_cen, s0, R2, c_full, R2_full = fit_cc(n, S)
    S_mid = S[n // 2 - 1]
    print(f"  {dt:.1f}s, chi_actual={chi_actual}, E0={E0:.8f}", flush=True)
    print(f"  c(central)={c_cen:.4f} R²={R2:.6f}, c(full)={c_full:.4f} R²={R2_full:.6f}", flush=True)
    print(f"  S_mid={S_mid:.4f}", flush=True)

    q7_data[f'n{n}'] = {
        'n': n, 'chi_max': chi, 'chi_actual': chi_actual,
        'c_central': c_cen, 'c_full': c_full,
        'R2_central': R2, 'R2_full': R2_full,
        's0': s0, 'S_mid': S_mid, 'S_profile': S,
        'E0': E0, 'time_s': dt
    }
    results['q7'] = q7_data
    save()

    if dt > 200:
        print(f"  Approaching time limit, stopping.", flush=True)
        break

# Summary
print(f"\n{'='*60}", flush=True)
print("SUMMARY: c(q=7) from DMRG entropy profile", flush=True)
print(f"{'='*60}", flush=True)
print(f"Complex CFT prediction: Re(c) = {c_cft.real:.4f}", flush=True)
print(f"\n{'n':>4} {'chi':>4} {'c(central)':>12} {'R²':>10} {'c(full)':>10} {'time':>8}", flush=True)
for k, v in q7_data.items():
    print(f"{v['n']:4d} {v['chi_actual']:4d} {v['c_central']:12.4f} {v['R2_central']:10.6f} {v['c_full']:10.4f} {v['time_s']:8.1f}s", flush=True)

# Convergence check
c_vals = [v['c_central'] for v in q7_data.values()]
if len(c_vals) >= 2:
    print(f"\nc convergence: {' -> '.join(f'{c:.4f}' for c in c_vals)}", flush=True)
    print(f"Last value: {c_vals[-1]:.4f}, CFT: {c_cft.real:.4f}, ratio: {c_vals[-1]/c_cft.real:.4f}", flush=True)

save()
print("\nSaved to results/sprint_079a_sq_potts_c_q7.json", flush=True)

from db_utils import record
for k, v in q7_data.items():
    record(sprint=79, model='sq_potts', q=7, n=v['n'],
           quantity='c_eff', value=v['c_central'],
           method='dmrg_CC_central', notes=f'gc=1/7, chi={v["chi_actual"]}, R2={v["R2_central"]:.4f}')
print("Recorded to DB.", flush=True)
