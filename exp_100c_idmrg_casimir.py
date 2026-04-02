#!/usr/bin/env python3
"""Sprint 100c: iDMRG bulk energy ε_∞ + periodic exact diag Casimir comparison.

Use iDMRG to get precise ε_∞ for q=2,5,7 S_q Potts at g_c=1/q.
Then combine with existing periodic-BC exact diag E₀(N) data (Sprint 099a)
to extract c_Casimir with ε_∞ fixed (eliminates largest fit parameter).

Periodic formula: E₀(N)/N = ε_∞ - πvc/(6N²)
With ε_∞ from iDMRG: slope of (E₀/N - ε_∞) vs 1/N² gives vc directly.
"""
import numpy as np
import json, time
from scipy.optimize import curve_fit

results = {
    'experiment': '100c_idmrg_casimir',
    'sprint': 100,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_100c_idmrg_casimir.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# DMRG setup
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

class SqPottsSite(Site):
    def __init__(self, q_val):
        leg = npc.LegCharge.from_trivial(q_val)
        Site.__init__(self, leg, [str(a) for a in range(q_val)], sort_charge=False)
        for a in range(q_val):
            P = np.zeros((q_val, q_val), dtype=complex); P[a, a] = 1.0
            self.add_op(f'P{a}', P)
        Sq_field = np.ones((q_val, q_val), dtype=complex) - np.eye(q_val, dtype=complex)
        self.add_op('SqField', Sq_field, hc='SqField')

class SqPottsInfinite(CouplingMPOModel):
    def init_sites(self, model_params):
        return SqPottsSite(model_params.get('q', 5))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 0.2)
        q_val = model_params.get('q', 5)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')

def run_idmrg(q_val, g, chi_max=100):
    """Run infinite DMRG to get bulk energy density ε_∞."""
    model = SqPottsInfinite({
        'L': 2, 'q': q_val, 'J': 1.0, 'g': g,
        'bc_MPS': 'infinite', 'bc_x': 'periodic',
    })
    np.random.seed(42 + q_val)
    init = [np.random.randint(q_val) for _ in range(2)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='infinite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-14,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 40,
    })
    E0, _ = eng.run()
    # E0 is total energy per unit cell (2 sites)
    eps_inf = float(E0) / 2.0
    chi_used = max(psi.chi)
    S_inf = float(psi.entanglement_entropy()[0])
    return eps_inf, chi_used, S_inf

# Reference values
ref = {
    2: {'rec': 0.500, 'v': 1.02, 'x_sigma': 0.123},
    5: {'rec': 1.138, 'v': 0.75, 'x_sigma': 0.136},
    7: {'rec': 1.351, 'v': 0.68, 'x_sigma': 0.132},
}

# Periodic exact diag data from Sprint 099a
periodic_data = {
    2: {
        'N': [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
        'source': 'sprint_099a',
    },
    5: {
        'N': [4, 5, 6, 7, 8, 9, 10],
        'source': 'sprint_099a',
    },
    7: {
        'N': [4, 5, 6, 7, 8],
        'source': 'sprint_099a',
    },
}

# Load the periodic E₀ data from Sprint 099a
try:
    with open("results/sprint_099a_casimir_dense.json") as f:
        s099 = json.load(f)
except FileNotFoundError:
    s099 = None

print("Sprint 100c: iDMRG ε_∞ + Periodic Casimir Comparison")
print("=" * 70, flush=True)

# Step 1: iDMRG for each q
for q_val in [2, 5, 7]:
    gc = 1.0 / q_val
    print(f"\n--- iDMRG: q={q_val}, g_c={gc:.6f} ---")

    for chi in [60, 100, 150]:
        print(f"  chi={chi}: ", end='', flush=True)
        t0 = time.time()
        eps_inf, chi_used, S_inf = run_idmrg(q_val, gc, chi_max=chi)
        t_e = time.time() - t0
        print(f"ε_∞ = {eps_inf:.14f}, S = {S_inf:.6f}, chi_used={chi_used}, t={t_e:.1f}s")

        if str(q_val) not in results['data']:
            results['data'][str(q_val)] = {'idmrg': [], 'gc': gc}
        results['data'][str(q_val)]['idmrg'].append({
            'chi_max': chi, 'chi_used': chi_used,
            'eps_inf': eps_inf, 'S_inf': S_inf, 'time_s': t_e,
        })
        save()

        if t_e > 120:
            print("  Time limit, using this chi")
            break

# Step 2: Load periodic E₀ data and combine with iDMRG ε_∞
print(f"\n{'='*70}")
print("PERIODIC BC CASIMIR WITH iDMRG ε_∞")
print(f"{'='*70}")

if s099:
    for q_val in [2, 5, 7]:
        q_str = str(q_val)
        r = ref[q_val]
        d = results['data'][q_str]

        # Best ε_∞ = highest chi
        eps_inf = d['idmrg'][-1]['eps_inf']
        gc = d['gc']

        # Get periodic E₀(N) from Sprint 099a
        q_key = q_str
        if q_key in s099.get('data', {}):
            q_data = s099['data'][q_key]
            N_arr = np.array(q_data['sizes'], dtype=float)
            E0_arr = np.array(q_data['E0_raw'])
        else:
            print(f"\nq={q_val}: No Sprint 099a data found")
            continue

        print(f"\nq={q_val}: ε_∞(iDMRG) = {eps_inf:.12f}")
        print(f"  Periodic E₀ at N = {[int(x) for x in N_arr]}")

        # E₀(N)/N - ε_∞ should scale as -πvc/(6N²)
        y = E0_arr / N_arr - eps_inf
        x = 1.0 / N_arr**2

        # 1-param fit: y = slope * x (intercept should be zero if ε_∞ is exact)
        slope_1p = np.dot(y, x) / np.dot(x, x)
        vc_1p = -slope_1p * 6 / np.pi
        c_1p = vc_1p / r['v']
        y_pred_1p = slope_1p * x
        R2_1p = 1 - np.sum((y - y_pred_1p)**2) / np.sum((y - np.mean(y))**2)

        print(f"  1-param (slope only): vc={vc_1p:.6f}, c={c_1p:.4f}, c/Rec={c_1p/r['rec']:.4f}, R²={R2_1p:.8f}")

        # 2-param fit: y = slope * x + intercept (intercept = finite-size ε_∞ correction)
        A = np.vstack([x, np.ones_like(x)]).T
        (slope_2p, intercept), _, _, _ = np.linalg.lstsq(A, y, rcond=None)
        vc_2p = -slope_2p * 6 / np.pi
        c_2p = vc_2p / r['v']
        y_pred_2p = slope_2p * x + intercept
        R2_2p = 1 - np.sum((y - y_pred_2p)**2) / np.sum((y - np.mean(y))**2)

        print(f"  2-param (+ intercept): vc={vc_2p:.6f}, c={c_2p:.4f}, c/Rec={c_2p/r['rec']:.4f}")
        print(f"    intercept = {intercept:.2e} (ε_∞ correction), R²={R2_2p:.8f}")

        # 2-param with 1/N⁴: y = slope*x + d*x²
        A3 = np.vstack([x, x**2]).T
        (slope_3, slope4), _, _, _ = np.linalg.lstsq(A3, y, rcond=None)
        vc_3 = -slope_3 * 6 / np.pi
        c_3 = vc_3 / r['v']
        y_pred_3 = slope_3 * x + slope4 * x**2
        R2_3 = 1 - np.sum((y - y_pred_3)**2) / np.sum((y - np.mean(y))**2)

        print(f"  2-param (+ 1/N⁴): vc={vc_3:.6f}, c={c_3:.4f}, c/Rec={c_3/r['rec']:.4f}, R²={R2_3:.8f}")

        # Pairwise from consecutive sizes
        print(f"  Pairwise vc (consecutive):")
        for i in range(len(N_arr) - 1):
            N1, N2 = N_arr[i], N_arr[i+1]
            dy = y[i+1] - y[i]
            dx = x[i+1] - x[i]
            vc_pair = -dy/dx * 6 / np.pi
            c_pair = vc_pair / r['v']
            print(f"    ({int(N1)},{int(N2)}): vc={vc_pair:.5f}, c={c_pair:.4f}, c/Rec={c_pair/r['rec']:.4f}")

        # Residuals
        print(f"  Residual (E₀/N - ε_∞) vs -πvc/(6N²):")
        for i in range(len(N_arr)):
            casimir = -np.pi * vc_2p / (6 * N_arr[i]**2)
            resid = y[i] - casimir
            print(f"    N={int(N_arr[i]):3d}: E/N-ε_∞={y[i]:.8e}, Casimir={casimir:.8e}, resid={resid:.2e}")

        d['casimir_1p'] = {'vc': float(vc_1p), 'c': float(c_1p), 'c_rec': float(c_1p/r['rec']), 'R2': float(R2_1p)}
        d['casimir_2p'] = {'vc': float(vc_2p), 'c': float(c_2p), 'c_rec': float(c_2p/r['rec']),
                          'R2': float(R2_2p), 'intercept': float(intercept)}
        d['casimir_2p_N4'] = {'vc': float(vc_3), 'c': float(c_3), 'c_rec': float(c_3/r['rec']), 'R2': float(R2_3)}
        save()

# Grand summary
print(f"\n{'='*70}")
print("SUMMARY: Periodic-BC Casimir with iDMRG ε_∞")
print(f"{'='*70}")
print(f"{'q':>3} {'ε_∞(iDMRG)':>16} {'c(2p+int)':>10} {'c/Rec':>7} {'c(2p+N4)':>10} {'c/Rec':>7}")
for q_val in [2, 5, 7]:
    d = results['data'][str(q_val)]
    eps = d['idmrg'][-1]['eps_inf']
    c2 = d.get('casimir_2p', {})
    c3 = d.get('casimir_2p_N4', {})
    if c2 and c3:
        print(f"{q_val:3d} {eps:16.12f} {c2['c']:10.4f} {c2['c_rec']:7.4f} {c3['c']:10.4f} {c3['c_rec']:7.4f}")

save()
from db_utils import record
for q_val in [2, 5, 7]:
    d = results['data'][str(q_val)]
    eps = d['idmrg'][-1]['eps_inf']
    record(sprint=100, model='sq_potts', q=q_val, n=999,
           quantity='eps_inf_idmrg', value=eps,
           method='idmrg',
           notes=f"chi={d['idmrg'][-1]['chi_used']}")
    if 'casimir_2p_N4' in d:
        c = d['casimir_2p_N4']
        record(sprint=100, model='sq_potts', q=q_val, n=999,
               quantity='c_casimir_periodic_idmrg', value=c['c'],
               method='periodic_casimir_idmrg_eps',
               notes=f"c/Rec={c['c_rec']:.4f}, R2={c['R2']:.8f}")

print("\nRecorded to DB.")
