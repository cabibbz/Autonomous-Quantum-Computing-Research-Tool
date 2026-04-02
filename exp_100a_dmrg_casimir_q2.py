#!/usr/bin/env python3
"""Sprint 100a: DMRG Casimir energy precision test for q=2 S_q Potts (= Ising).

Known answer: c = 0.5 exactly.
Open BC Casimir formula: E0 = eps_inf * N + e_s - pi*v*c/(24*N) + O(1/N^2)

Strategy:
  1. Run finite DMRG at g_c = 1/2 for N = 8, 10, 12, 14, 16, 20, 24, 30, 40, 50
  2. Also get exact diag E0 at N=8,10,12,14 (open BC) as cross-check
  3. Extract c from pairwise (E0 - eps_inf*N - e_s) vs 1/N

For open BC, eps_inf and v are the same as periodic. We need e_s (surface energy)
as a fit parameter. Fit: E0/N = eps_inf + e_s/N - pi*v*c/(24*N^2) + d/N^3
Or better: fit E0 = a*N + b + c_coeff/N + d/N^2.
"""
import numpy as np
import json, time
from scipy.optimize import curve_fit

q = 2
gc = 0.5  # g_c = 1/q

results = {
    'experiment': '100a_dmrg_casimir_q2',
    'sprint': 100,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'gc': gc,
    'known_c': 0.5,
    'data_dmrg': [],
    'data_exact': [],
}

def save():
    with open("results/sprint_100a_dmrg_casimir_q2.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# ---- Exact diag for cross-check (open BC) ----
from scipy.sparse import coo_matrix, csr_matrix, kron, eye
from gpu_utils import eigsh

def build_sq_potts_open(n, q_val, g):
    """Build S_q Potts Hamiltonian on OPEN chain."""
    dim = q_val**n
    all_idx = np.arange(dim, dtype=np.int64)
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q_val
        tmp //= q_val
    powers = q_val ** np.arange(n, dtype=np.int64)

    # Diagonal: coupling -J*delta(s_i, s_{i+1}) for i=0..n-2 (OPEN)
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n - 1):  # open BC: no wrap
        diag_vals -= (digits[:, site] == digits[:, site + 1]).astype(np.float64)

    rows_list = [all_idx]
    cols_list = [all_idx]
    vals_list = [diag_vals]

    # Off-diagonal: field
    for site in range(n):
        pw = powers[site]
        old_digit = digits[:, site]
        for k in range(1, q_val):
            new_digit = (old_digit + k) % q_val
            delta = (new_digit.astype(np.int64) - old_digit.astype(np.int64)) * pw
            new_idx = all_idx + delta
            rows_list.append(new_idx)
            cols_list.append(all_idx)
            vals_list.append(np.full(dim, -g, dtype=np.float64))

    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    vals = np.concatenate(vals_list)
    return csr_matrix(coo_matrix((vals, (rows, cols)), shape=(dim, dim)))

# ---- DMRG setup ----
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

def run_dmrg(n, q_val, g, chi_max=100):
    """Run finite DMRG, return E0 and max bond dimension used."""
    model = SqPottsChain({'L': n, 'q': q_val, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n)
    init = [np.random.randint(q_val) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': 40,
    })
    E0, _ = eng.run()
    chi_used = max(psi.chi)
    return float(E0), chi_used

print("Sprint 100a: DMRG Casimir Precision Test — q=2 (known c=0.5)")
print("=" * 70, flush=True)

# Step 1: Exact diag at small sizes for cross-check
print("\n--- Exact diagonalization (open BC) ---")
exact_sizes = [8, 10, 12, 14]
for n in exact_sizes:
    dim = q**n
    print(f"  n={n}, dim={dim:,} ... ", end='', flush=True)
    t0 = time.time()
    H = build_sq_potts_open(n, q, gc)
    evals, _ = eigsh(H, k=4, which='SA')
    t_elapsed = time.time() - t0
    evals = np.sort(np.real(evals))
    E0 = float(evals[0])
    gap = float(evals[1] - evals[0])
    print(f"E0={E0:.12f}, gap={gap:.8f}, t={t_elapsed:.1f}s")
    results['data_exact'].append({'n': n, 'E0': E0, 'gap': gap})
    save()

# Step 2: DMRG at various sizes
print("\n--- DMRG (open BC) ---")
dmrg_sizes = [8, 10, 12, 14, 16, 20, 24, 30, 40, 50]
chi_schedule = {8: 60, 10: 60, 12: 80, 14: 80, 16: 100, 20: 100,
                24: 120, 30: 120, 40: 150, 50: 150}

for n in dmrg_sizes:
    chi = chi_schedule[n]
    print(f"  n={n}, chi_max={chi} ... ", end='', flush=True)
    t0 = time.time()
    E0_dmrg, chi_used = run_dmrg(n, q, gc, chi_max=chi)
    t_elapsed = time.time() - t0

    # Compare with exact if available
    exact_entry = next((e for e in results['data_exact'] if e['n'] == n), None)
    if exact_entry:
        err = abs(E0_dmrg - exact_entry['E0'])
        print(f"E0={E0_dmrg:.12f}, chi_used={chi_used}, t={t_elapsed:.1f}s, "
              f"exact_err={err:.2e}")
    else:
        print(f"E0={E0_dmrg:.12f}, chi_used={chi_used}, t={t_elapsed:.1f}s")

    results['data_dmrg'].append({
        'n': n, 'E0': E0_dmrg, 'chi_max': chi, 'chi_used': chi_used,
        'time_s': t_elapsed,
    })
    save()

# Step 3: Fit open-BC Casimir formula
# E0 = eps_inf * N + e_s + A/N + B/N^2
# where A = -pi*v*c/24
print("\n--- Casimir fit (open BC) ---")
N_arr = np.array([d['n'] for d in results['data_dmrg']], dtype=float)
E0_arr = np.array([d['E0'] for d in results['data_dmrg']])

# 3-parameter fit: E0 = a*N + b + c_coeff/N
def model_3p(N, a, b, cc):
    return a * N + b + cc / N

popt3, pcov3 = curve_fit(model_3p, N_arr, E0_arr, p0=[-1.5, 0, 0.01])
eps_inf, e_s, A_coeff = popt3
y_pred3 = model_3p(N_arr, *popt3)
ss_res3 = np.sum((E0_arr - y_pred3)**2)
ss_tot = np.sum((E0_arr - np.mean(E0_arr))**2)
R2_3 = 1 - ss_res3 / ss_tot

print(f"  3-param: eps_inf={eps_inf:.10f}, e_s={e_s:.8f}, A={A_coeff:.8f}")
print(f"  R² = {R2_3:.12f}")

# 4-parameter fit: E0 = a*N + b + c_coeff/N + d/N^2
def model_4p(N, a, b, cc, d):
    return a * N + b + cc / N + d / N**2

popt4, pcov4 = curve_fit(model_4p, N_arr, E0_arr, p0=[popt3[0], popt3[1], popt3[2], 0])
eps_inf4, e_s4, A_coeff4, B_coeff = popt4
y_pred4 = model_4p(N_arr, *popt4)
ss_res4 = np.sum((E0_arr - y_pred4)**2)
R2_4 = 1 - ss_res4 / ss_tot

print(f"  4-param: eps_inf={eps_inf4:.10f}, e_s={e_s4:.8f}, "
      f"A={A_coeff4:.8f}, B={B_coeff:.8f}")
print(f"  R² = {R2_4:.12f}")

# Extract c from A = -pi*v*c/24
# For q=2 at g_c=0.5: v ~ 1.0 (exact for TFI at g_c: v = J*sin(pi/(2*1)) = pi/2...
# actually for S_q Potts at g_c=1/q, v from gap*N/(2*pi*x_sigma)
# Use known v for q=2 Ising: v = pi*J*sin(pi*g_c) at g_c = 1 (TFI convention)
# For our S_q convention, v from Sprint 082: v(q=2) = 1.02
v_q2 = 1.02  # from gap*N/(2*pi*x_sigma), Sprint 082

c_3p = -A_coeff * 24 / (np.pi * v_q2)
c_4p = -A_coeff4 * 24 / (np.pi * v_q2)

print(f"\n  Using v(q=2) = {v_q2}")
print(f"  c (3-param) = {c_3p:.4f}  (known: 0.500, ratio: {c_3p/0.5:.4f})")
print(f"  c (4-param) = {c_4p:.4f}  (known: 0.500, ratio: {c_4p/0.5:.4f})")

# Pairwise extraction: consecutive sizes
print("\n--- Pairwise c extraction ---")
pairwise = []
for i in range(len(N_arr) - 1):
    N1, N2 = N_arr[i], N_arr[i+1]
    E1, E2 = E0_arr[i], E0_arr[i+1]
    # E = a*N + b + A/N  =>  E2-E1 = a*(N2-N1) + A*(1/N2 - 1/N1)
    # Need 3 consecutive points to eliminate a:
    # Or use the 4-param eps_inf and subtract linear part
    pass

# Better: subtract the bulk part and fit the remainder
# (E0 - eps_inf*N) should be  e_s + A/N + ...
E0_sub = E0_arr - eps_inf4 * N_arr
# Pairwise from E0_sub:
for i in range(len(N_arr) - 1):
    N1, N2 = N_arr[i], N_arr[i+1]
    e1, e2 = E0_sub[i], E0_sub[i+1]
    # e = b + A/N  =>  A_pair = (e2 - e1) / (1/N2 - 1/N1)
    A_pair = (e2 - e1) / (1/N2 - 1/N1)
    c_pair = -A_pair * 24 / (np.pi * v_q2)
    pairwise.append({
        'pair': f'({int(N1)},{int(N2)})',
        'A_pair': float(A_pair),
        'c_pair': float(c_pair),
        'c_over_exact': float(c_pair / 0.5),
    })
    print(f"  ({int(N1)},{int(N2)}): A={A_pair:.6f}, c={c_pair:.4f}, c/0.5={c_pair/0.5:.4f}")

# Residuals
print("\n--- Residuals (4-param fit) ---")
for i, n in enumerate([int(x) for x in N_arr]):
    resid = E0_arr[i] - y_pred4[i]
    print(f"  n={n:3d}: E0={E0_arr[i]:.10f}, resid={resid:+.2e}")

results['fit_3param'] = {
    'eps_inf': float(eps_inf), 'e_s': float(e_s),
    'A': float(A_coeff), 'R2': float(R2_3),
    'c_implied': float(c_3p),
}
results['fit_4param'] = {
    'eps_inf': float(eps_inf4), 'e_s': float(e_s4),
    'A': float(A_coeff4), 'B': float(B_coeff), 'R2': float(R2_4),
    'c_implied': float(c_4p),
}
results['pairwise'] = pairwise
results['v_used'] = v_q2
save()

# Record to DB
from db_utils import record
record(sprint=100, model='sq_potts', q=2, n=50,
       quantity='c_casimir_open', value=c_4p,
       method='dmrg_casimir_open_bc',
       notes=f"4-param fit, v={v_q2}, R2={R2_4:.10f}, "
             f"N=8-50, known c=0.5, ratio={c_4p/0.5:.4f}")

print(f"\nSaved to results/sprint_100a_dmrg_casimir_q2.json")
print("Recorded to DB.")
