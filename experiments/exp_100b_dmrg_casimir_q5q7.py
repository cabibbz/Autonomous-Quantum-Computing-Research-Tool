#!/usr/bin/env python3
"""Sprint 100b: DMRG Casimir energy for q=5 and q=7 S_q Potts at large N (open BC).

q=5: Re(c)=1.138, c_eff=1.152 (only 1.2% diff — hard to distinguish)
q=7: Re(c)=1.351, c_eff=1.059 (22% diff — easy to distinguish)

For q=5, DMRG with d=5 is ~10x slower than q=2. Use chi=40-60 and N=6-20.
For q=7, DMRG with d=7 is ~30x slower than q=2. Use chi=30-40 and N=6-14.

Open BC formula: E0 = eps_inf*N + e_s - pi*v*c/(24*N) + B/N^2
"""
import numpy as np
import json, time
from scipy.optimize import curve_fit

results = {
    'experiment': '100b_dmrg_casimir_q5q7',
    'sprint': 100,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_100b_dmrg_casimir_q5q7.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# Exact diag for cross-check
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh

def build_sq_potts_open(n, q_val, g):
    dim = q_val**n
    all_idx = np.arange(dim, dtype=np.int64)
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q_val
        tmp //= q_val
    powers = q_val ** np.arange(n, dtype=np.int64)
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n - 1):
        diag_vals -= (digits[:, site] == digits[:, site + 1]).astype(np.float64)
    rows_list = [all_idx]; cols_list = [all_idx]; vals_list = [diag_vals]
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
    return csr_matrix(coo_matrix((np.concatenate(vals_list),
                                   (np.concatenate(rows_list), np.concatenate(cols_list))),
                                  shape=(dim, dim)))

# DMRG setup
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
        return SqPottsSite(model_params.get('q', 5))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 0.2)
        q_val = model_params.get('q', 5)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')

def run_dmrg(n, q_val, g, chi_max=40, max_sweeps=20):
    model = SqPottsChain({'L': n, 'q': q_val, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n + q_val * 100)
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

# Reference values
ref = {
    5: {'rec': 1.138, 'c_eff': 1.152, 'v': 0.75, 'x_sigma': 0.136},
    7: {'rec': 1.351, 'c_eff': 1.059, 'v': 0.68, 'x_sigma': 0.132},
}

print("Sprint 100b: DMRG Casimir for q=5 and q=7")
print("=" * 70, flush=True)

# ---- q=5 ----
q = 5
gc = 1.0 / q
print(f"\n{'='*70}")
print(f"q={q}, g_c={gc}, Re(c)={ref[q]['rec']}, c_eff={ref[q]['c_eff']}")
print(f"{'='*70}")

# Exact diag at small sizes (open BC)
print("\n--- Exact diag (open BC) ---")
exact_q5 = []
for n in [6, 7, 8]:
    dim = q**n
    if dim > 500000:
        break
    print(f"  n={n}, dim={dim:,} ... ", end='', flush=True)
    t0 = time.time()
    H = build_sq_potts_open(n, q, gc)
    evals, _ = eigsh(H, k=2, which='SA')
    t_e = time.time() - t0
    evals = np.sort(np.real(evals))
    E0 = float(evals[0])
    print(f"E0={E0:.12f}, t={t_e:.1f}s")
    exact_q5.append({'n': n, 'E0': E0})

# GPU exact diag for n=9,10
for n in [9, 10]:
    dim = q**n
    print(f"  n={n}, dim={dim:,} (GPU) ... ", end='', flush=True)
    t0 = time.time()
    H = build_sq_potts_open(n, q, gc)
    evals, _ = eigsh(H, k=2, which='SA')
    t_e = time.time() - t0
    evals = np.sort(np.real(evals))
    E0 = float(evals[0])
    print(f"E0={E0:.12f}, t={t_e:.1f}s")
    exact_q5.append({'n': n, 'E0': E0})
    save()

results['data']['5'] = {'exact': exact_q5, 'dmrg': []}
save()

# DMRG for q=5
print("\n--- DMRG (open BC) ---")
dmrg_q5_schedule = [
    (6, 30, 15), (7, 30, 15), (8, 40, 15), (9, 40, 15), (10, 40, 15),
    (12, 50, 18), (14, 50, 18), (16, 60, 18), (20, 60, 18),
]

for n, chi, sweeps in dmrg_q5_schedule:
    print(f"  n={n}, chi={chi}, sweeps={sweeps} ... ", end='', flush=True)
    t0 = time.time()
    E0, chi_used = run_dmrg(n, q, gc, chi_max=chi, max_sweeps=sweeps)
    t_e = time.time() - t0

    # Cross-check
    exact = next((e for e in exact_q5 if e['n'] == n), None)
    if exact:
        err = abs(E0 - exact['E0'])
        print(f"E0={E0:.12f}, chi={chi_used}, t={t_e:.1f}s, exact_err={err:.2e}")
    else:
        print(f"E0={E0:.12f}, chi={chi_used}, t={t_e:.1f}s")

    results['data']['5']['dmrg'].append({
        'n': n, 'E0': E0, 'chi_max': chi, 'chi_used': chi_used, 'time_s': t_e,
    })
    save()

    if t_e > 250:
        print("  WARNING: approaching time limit, stopping q=5 DMRG")
        break

# ---- q=7 ----
q = 7
gc = 1.0 / q
print(f"\n{'='*70}")
print(f"q={q}, g_c={gc:.6f}, Re(c)={ref[q]['rec']}, c_eff={ref[q]['c_eff']}")
print(f"{'='*70}")

# Exact diag (open BC) — q=7 sizes limited
print("\n--- Exact diag (open BC) ---")
exact_q7 = []
for n in [5, 6, 7]:
    dim = q**n
    print(f"  n={n}, dim={dim:,} ... ", end='', flush=True)
    t0 = time.time()
    H = build_sq_potts_open(n, q, gc)
    evals, _ = eigsh(H, k=2, which='SA')
    t_e = time.time() - t0
    evals = np.sort(np.real(evals))
    E0 = float(evals[0])
    print(f"E0={E0:.12f}, t={t_e:.1f}s")
    exact_q7.append({'n': n, 'E0': E0})

# n=8 via GPU
n = 8
dim = q**n
print(f"  n={n}, dim={dim:,} (GPU) ... ", end='', flush=True)
t0 = time.time()
H = build_sq_potts_open(n, q, gc)
evals, _ = eigsh(H, k=2, which='SA')
t_e = time.time() - t0
evals = np.sort(np.real(evals))
E0 = float(evals[0])
print(f"E0={E0:.12f}, t={t_e:.1f}s")
exact_q7.append({'n': n, 'E0': E0})

results['data']['7'] = {'exact': exact_q7, 'dmrg': []}
save()

# DMRG for q=7
print("\n--- DMRG (open BC) ---")
dmrg_q7_schedule = [
    (6, 30, 15), (7, 30, 15), (8, 30, 15),
    (10, 40, 15), (12, 40, 15), (14, 50, 15),
]

for n, chi, sweeps in dmrg_q7_schedule:
    print(f"  n={n}, chi={chi}, sweeps={sweeps} ... ", end='', flush=True)
    t0 = time.time()
    E0, chi_used = run_dmrg(n, q, gc, chi_max=chi, max_sweeps=sweeps)
    t_e = time.time() - t0

    exact = next((e for e in exact_q7 if e['n'] == n), None)
    if exact:
        err = abs(E0 - exact['E0'])
        print(f"E0={E0:.12f}, chi={chi_used}, t={t_e:.1f}s, exact_err={err:.2e}")
    else:
        print(f"E0={E0:.12f}, chi={chi_used}, t={t_e:.1f}s")

    results['data']['7']['dmrg'].append({
        'n': n, 'E0': E0, 'chi_max': chi, 'chi_used': chi_used, 'time_s': t_e,
    })
    save()

    if t_e > 250:
        print("  WARNING: approaching time limit, stopping q=7 DMRG")
        break

# --- Analysis for both q ---
print(f"\n{'='*70}")
print("ANALYSIS")
print(f"{'='*70}")

for q_str in ['5', '7']:
    q_val = int(q_str)
    r = ref[q_val]
    d = results['data'][q_str]
    dmrg_data = d['dmrg']
    if len(dmrg_data) < 3:
        print(f"\nq={q_val}: insufficient data ({len(dmrg_data)} points)")
        continue

    N_arr = np.array([x['n'] for x in dmrg_data], dtype=float)
    E0_arr = np.array([x['E0'] for x in dmrg_data])

    print(f"\nq={q_val}: {len(N_arr)} points, N={[int(x) for x in N_arr]}")

    # 3-param fit: E0 = a*N + b + A/N
    def model_3p(N, a, b, cc):
        return a * N + b + cc / N

    popt3, _ = curve_fit(model_3p, N_arr, E0_arr)
    y_pred3 = model_3p(N_arr, *popt3)
    ss_tot = np.sum((E0_arr - np.mean(E0_arr))**2)
    R2_3 = 1 - np.sum((E0_arr - y_pred3)**2) / ss_tot
    c_3 = -popt3[2] * 24 / (np.pi * r['v'])

    print(f"  3-param: eps_inf={popt3[0]:.8f}, e_s={popt3[1]:.6f}, A={popt3[2]:.6f}")
    print(f"  R² = {R2_3:.12f}")
    print(f"  c(3p) = {c_3:.4f}, c/Re(c) = {c_3/r['rec']:.4f}, c/c_eff = {c_3/r['c_eff']:.4f}")

    # 4-param if enough data
    if len(N_arr) >= 4:
        def model_4p(N, a, b, cc, dd):
            return a * N + b + cc / N + dd / N**2
        popt4, pcov4 = curve_fit(model_4p, N_arr, E0_arr, p0=[*popt3, 0])
        y_pred4 = model_4p(N_arr, *popt4)
        R2_4 = 1 - np.sum((E0_arr - y_pred4)**2) / ss_tot
        c_4 = -popt4[2] * 24 / (np.pi * r['v'])
        c_4_err = np.sqrt(pcov4[2, 2]) * 24 / (np.pi * r['v'])
        print(f"  4-param: A={popt4[2]:.6f}, B={popt4[3]:.6f}")
        print(f"  R² = {R2_4:.12f}")
        print(f"  c(4p) = {c_4:.4f} ± {c_4_err:.4f}")
        print(f"  c/Re(c) = {c_4/r['rec']:.4f}, c/c_eff = {c_4/r['c_eff']:.4f}")
        d['fit_4param'] = {
            'c_implied': float(c_4), 'c_error': float(c_4_err),
            'c_over_rec': float(c_4/r['rec']),
            'c_over_ceff': float(c_4/r['c_eff']),
            'R2': float(R2_4),
        }

    # Pairwise (using eps_inf from fit)
    eps_best = popt4[0] if len(N_arr) >= 4 else popt3[0]
    E0_sub = E0_arr - eps_best * N_arr
    print(f"\n  Pairwise c extraction:")
    pw_list = []
    for i in range(len(N_arr) - 1):
        N1, N2 = N_arr[i], N_arr[i+1]
        A_pair = (E0_sub[i+1] - E0_sub[i]) / (1/N2 - 1/N1)
        c_pair = -A_pair * 24 / (np.pi * r['v'])
        pw_list.append({'pair': f'({int(N1)},{int(N2)})', 'c': float(c_pair),
                        'c_rec': float(c_pair/r['rec']), 'c_ceff': float(c_pair/r['c_eff'])})
        print(f"    ({int(N1):2d},{int(N2):2d}): c={c_pair:.4f}, "
              f"c/Re(c)={c_pair/r['rec']:.4f}, c/c_eff={c_pair/r['c_eff']:.4f}")

    d['fit_3param'] = {'c_implied': float(c_3), 'c_over_rec': float(c_3/r['rec']), 'R2': float(R2_3)}
    d['pairwise'] = pw_list

save()

# Record to DB
from db_utils import record
for q_str in ['5', '7']:
    q_val = int(q_str)
    d = results['data'][q_str]
    if 'fit_4param' in d:
        f = d['fit_4param']
        n_max = max(x['n'] for x in d['dmrg'])
        record(sprint=100, model='sq_potts', q=q_val, n=n_max,
               quantity='c_casimir_open', value=f['c_implied'],
               error=f['c_error'], method='dmrg_casimir_open_bc',
               notes=f"c/Rec={f['c_over_rec']:.4f}, c/c_eff={f['c_over_ceff']:.4f}")

print(f"\nSaved to results/sprint_100b_dmrg_casimir_q5q7.json")
print("Recorded to DB.")
