#!/usr/bin/env python3
"""Sprint 100a (part 2): Extend DMRG Casimir q=2 to n=24,30 with relaxed convergence.

Previous run got n=8-20. The n=24 chi=120 run with max_E_err=1e-12 hung.
Fix: relax to 1e-10, fewer sweeps. Also do analysis on all data.
"""
import numpy as np
import json, time
from scipy.optimize import curve_fit

q = 2
gc = 0.5
v_q2 = 1.02  # from Sprint 082

# Load existing results
with open("results/sprint_100a_dmrg_casimir_q2.json") as f:
    results = json.load(f)

def save():
    with open("results/sprint_100a_dmrg_casimir_q2.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
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

class SqPottsChain(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, model_params):
        return SqPottsSite(model_params.get('q', 2))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', gc)
        q_val = model_params.get('q', 2)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')

def run_dmrg(n, q_val, g, chi_max=80, max_sweeps=20):
    model = SqPottsChain({'L': n, 'q': q_val, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n)
    init = [np.random.randint(q_val) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': max_sweeps,
    })
    E0, _ = eng.run()
    chi_used = max(psi.chi)
    return float(E0), chi_used

print("Sprint 100a2: Extending DMRG q=2 to n=24,30")
print("=" * 70, flush=True)

# Extend with new sizes
new_sizes = [(24, 80, 20), (30, 100, 20)]  # (n, chi_max, max_sweeps)

for n, chi, sweeps in new_sizes:
    # Skip if already have this size
    existing = [d for d in results['data_dmrg'] if d['n'] == n]
    if existing:
        print(f"  n={n} already exists, skipping")
        continue

    print(f"  n={n}, chi_max={chi}, max_sweeps={sweeps} ... ", end='', flush=True)
    t0 = time.time()
    E0_dmrg, chi_used = run_dmrg(n, q, gc, chi_max=chi, max_sweeps=sweeps)
    t_elapsed = time.time() - t0
    print(f"E0={E0_dmrg:.12f}, chi_used={chi_used}, t={t_elapsed:.1f}s")

    results['data_dmrg'].append({
        'n': n, 'E0': E0_dmrg, 'chi_max': chi, 'chi_used': chi_used,
        'time_s': t_elapsed,
    })
    save()

# Sort by n
results['data_dmrg'].sort(key=lambda x: x['n'])

# --- Full analysis on all data ---
print("\n--- Casimir fit analysis (all DMRG data) ---")
N_arr = np.array([d['n'] for d in results['data_dmrg']], dtype=float)
E0_arr = np.array([d['E0'] for d in results['data_dmrg']])

print(f"Sizes: {[int(x) for x in N_arr]}")
print(f"E0 values: {[f'{x:.10f}' for x in E0_arr]}")

# Open BC: E0 = eps_inf * N + e_s + A/N + B/N^2
# A = -pi*v*c/24

# 3-parameter fit: E0 = a*N + b + c_coeff/N
def model_3p(N, a, b, cc):
    return a * N + b + cc / N

popt3, pcov3 = curve_fit(model_3p, N_arr, E0_arr, p0=[-1.27, 0, 0.01])
eps_inf3, e_s3, A3 = popt3
y_pred3 = model_3p(N_arr, *popt3)
ss_tot = np.sum((E0_arr - np.mean(E0_arr))**2)
ss_res3 = np.sum((E0_arr - y_pred3)**2)
R2_3 = 1 - ss_res3 / ss_tot

c_3p = -A3 * 24 / (np.pi * v_q2)
print(f"\n3-param fit: E0 = {eps_inf3:.10f}*N + {e_s3:.8f} + {A3:.8f}/N")
print(f"  R² = {R2_3:.14f}")
print(f"  c = {c_3p:.6f} (known: 0.500, ratio: {c_3p/0.5:.4f})")

# 4-parameter fit: + B/N^2
def model_4p(N, a, b, cc, d):
    return a * N + b + cc / N + d / N**2

popt4, pcov4 = curve_fit(model_4p, N_arr, E0_arr, p0=[*popt3, 0])
eps_inf4, e_s4, A4, B4 = popt4
y_pred4 = model_4p(N_arr, *popt4)
ss_res4 = np.sum((E0_arr - y_pred4)**2)
R2_4 = 1 - ss_res4 / ss_tot

c_4p = -A4 * 24 / (np.pi * v_q2)
print(f"\n4-param fit: E0 = {eps_inf4:.10f}*N + {e_s4:.8f} + {A4:.8f}/N + {B4:.8f}/N²")
print(f"  R² = {R2_4:.14f}")
print(f"  c = {c_4p:.6f} (known: 0.500, ratio: {c_4p/0.5:.4f})")
print(f"  B/A = {B4/A4:.4f} (measures subleading correction strength)")

# Pairwise extraction
print("\n--- Pairwise c from consecutive sizes ---")
# Subtract bulk: E_sub = E0 - eps_inf*N = e_s + A/N + ...
E0_sub = E0_arr - eps_inf4 * N_arr
pairwise = []
for i in range(len(N_arr) - 1):
    N1, N2 = N_arr[i], N_arr[i+1]
    e1, e2 = E0_sub[i], E0_sub[i+1]
    A_pair = (e2 - e1) / (1/N2 - 1/N1)
    c_pair = -A_pair * 24 / (np.pi * v_q2)
    pairwise.append({
        'pair': f'({int(N1)},{int(N2)})',
        'A_pair': float(A_pair),
        'c_pair': float(c_pair),
        'c_over_exact': float(c_pair / 0.5),
    })
    print(f"  ({int(N1):2d},{int(N2):2d}): A={A_pair:.8f}, c={c_pair:.4f}, c/0.5={c_pair/0.5:.4f}")

# Residuals
print("\n--- Residuals (4-param) ---")
for i in range(len(N_arr)):
    resid = E0_arr[i] - y_pred4[i]
    print(f"  n={int(N_arr[i]):3d}: resid={resid:+.2e}")

# Save everything
results['fit_3param'] = {
    'eps_inf': float(eps_inf3), 'e_s': float(e_s3),
    'A': float(A3), 'R2': float(R2_3), 'c_implied': float(c_3p),
}
results['fit_4param'] = {
    'eps_inf': float(eps_inf4), 'e_s': float(e_s4),
    'A': float(A4), 'B': float(B4), 'R2': float(R2_4), 'c_implied': float(c_4p),
}
results['pairwise'] = pairwise
results['v_used'] = v_q2
save()

# Record to DB
from db_utils import record
record(sprint=100, model='sq_potts', q=2, n=int(max(N_arr)),
       quantity='c_casimir_open', value=c_4p,
       method='dmrg_casimir_open_bc',
       notes=f"4p fit, v={v_q2}, R2={R2_4:.12f}, "
             f"N=8-{int(max(N_arr))}, known c=0.5, c/0.5={c_4p/0.5:.4f}")

print(f"\nSaved. Recorded to DB.")
print(f"\n{'='*70}")
print(f"VERDICT: DMRG open-BC Casimir extraction for q=2")
print(f"  c(4-param) = {c_4p:.4f}, deviation from 0.500: {abs(c_4p-0.5)/0.5*100:.2f}%")
print(f"  Pairwise at largest sizes converges to: {pairwise[-1]['c_pair']:.4f}")
print(f"{'='*70}")
