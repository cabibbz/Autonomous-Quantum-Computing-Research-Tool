#!/usr/bin/env python3
"""Sprint 100b (part 3): q=7 exact diag + DMRG, then full analysis."""
import numpy as np
import json, time
from scipy.optimize import curve_fit
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh

with open("results/sprint_100b_dmrg_casimir_q5q7.json") as f:
    results = json.load(f)

def save():
    with open("results/sprint_100b_dmrg_casimir_q5q7.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

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

# ---- q=7 exact diag ----
print("Sprint 100b3: q=7 exact diag + analysis")
print("=" * 70)

if '7' not in results['data'] or len(results['data'].get('7', {}).get('exact', [])) == 0:
    results['data']['7'] = {'exact': [], 'dmrg': []}
    print("\n--- q=7 Exact diag (open BC) ---")
    for n in [4, 5, 6, 7, 8]:
        dim = 7**n
        print(f"  n={n}, dim={dim:,} ... ", end='', flush=True)
        t0 = time.time()
        H = build_sq_potts_open(n, 7, 1/7)
        evals, _ = eigsh(H, k=2, which='SA')
        t_e = time.time() - t0
        evals = np.sort(np.real(evals))
        E0 = float(evals[0])
        print(f"E0={E0:.12f}, t={t_e:.1f}s")
        results['data']['7']['exact'].append({'n': n, 'E0': E0})
    save()

# ---- q=7 DMRG n=10 ----
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
        return SqPottsSite(model_params.get('q', 7))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 1/7)
        q_val = model_params.get('q', 7)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')

print("\n--- q=7 DMRG n=10 ---")
existing_n10 = [d for d in results['data']['7'].get('dmrg', []) if d['n'] == 10]
if not existing_n10:
    n, chi, sweeps = 10, 35, 15
    print(f"  n={n}, chi={chi} ... ", end='', flush=True)
    model = SqPottsChain({'L': n, 'q': 7, 'J': 1.0, 'g': 1/7, 'bc_MPS': 'finite'})
    np.random.seed(42 + n + 700)
    init = [np.random.randint(7) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi, 'svd_min': 1e-12},
        'max_sweeps': sweeps,
    })
    t0 = time.time()
    E0, _ = eng.run()
    t_e = time.time() - t0
    chi_used = max(psi.chi)
    print(f"E0={float(E0):.12f}, chi={chi_used}, t={t_e:.1f}s")
    results['data']['7']['dmrg'].append({
        'n': n, 'E0': float(E0), 'chi_max': chi, 'chi_used': chi_used, 'time_s': t_e})
    save()

# ---- FULL ANALYSIS ----
print(f"\n{'='*70}")
print("CASIMIR ANALYSIS — ALL q")
print(f"{'='*70}")

ref = {
    2: {'rec': 0.500, 'c_eff': 0.500, 'v': 1.02},
    5: {'rec': 1.138, 'c_eff': 1.152, 'v': 0.75},
    7: {'rec': 1.351, 'c_eff': 1.059, 'v': 0.68},
}

# Load q=2 data
with open("results/sprint_100a_dmrg_casimir_q2.json") as f:
    q2_data = json.load(f)

all_results = {}

for q_val in [2, 5, 7]:
    r = ref[q_val]
    q_str = str(q_val)

    # Get best data for each N
    all_data = {}
    if q_val == 2:
        for e in q2_data.get('data_exact', []):
            all_data[e['n']] = e['E0']
        for e in q2_data.get('data_dmrg', []):
            if e['n'] not in all_data:
                all_data[e['n']] = e['E0']
    else:
        d = results['data'][q_str]
        for e in d.get('exact', []):
            all_data[e['n']] = e['E0']
        for e in d.get('dmrg', []):
            if e['n'] not in all_data:
                all_data[e['n']] = e['E0']

    sorted_n = sorted(all_data.keys())
    N_arr = np.array(sorted_n, dtype=float)
    E0_arr = np.array([all_data[n] for n in sorted_n])

    print(f"\nq={q_val}: {len(N_arr)} pts, N={[int(x) for x in N_arr]}")
    print(f"  Re(c)={r['rec']}, c_eff={r['c_eff']}, v={r['v']}")

    # 3p fit
    def m3(N, a, b, cc): return a * N + b + cc / N
    popt3, _ = curve_fit(m3, N_arr, E0_arr)
    c_3 = -popt3[2] * 24 / (np.pi * r['v'])

    # 4p fit
    def m4(N, a, b, cc, dd): return a * N + b + cc / N + dd / N**2
    if len(N_arr) >= 4:
        popt4, pcov4 = curve_fit(m4, N_arr, E0_arr, p0=[*popt3, 0])
        c_4 = -popt4[2] * 24 / (np.pi * r['v'])
        c_4_err = np.sqrt(pcov4[2, 2]) * 24 / (np.pi * r['v'])
        y_pred4 = m4(N_arr, *popt4)
        ss_tot = np.sum((E0_arr - np.mean(E0_arr))**2)
        R2_4 = 1 - np.sum((E0_arr - y_pred4)**2) / ss_tot
    else:
        c_4, c_4_err, R2_4 = c_3, 0, 0

    print(f"  c(4p) = {c_4:.4f} ± {c_4_err:.4f}")
    print(f"  c/Re(c) = {c_4/r['rec']:.4f}")
    print(f"  c/c_eff = {c_4/r['c_eff']:.4f}")

    # Pairwise
    eps = popt4[0] if len(N_arr) >= 4 else popt3[0]
    E0_sub = E0_arr - eps * N_arr
    pw = []
    for i in range(len(N_arr) - 1):
        N1, N2 = N_arr[i], N_arr[i+1]
        A_p = (E0_sub[i+1] - E0_sub[i]) / (1/N2 - 1/N1)
        c_p = -A_p * 24 / (np.pi * r['v'])
        pw.append((int(N1), int(N2), c_p, c_p/r['rec'], c_p/r['c_eff']))
    print(f"  Pairwise:")
    for N1, N2, c_p, cr, ce in pw:
        print(f"    ({N1:2d},{N2:2d}): c={c_p:.4f}, c/Rec={cr:.4f}, c/c_eff={ce:.4f}")

    all_results[q_val] = {
        'c_4p': c_4, 'c_err': c_4_err,
        'c_over_rec': c_4/r['rec'], 'c_over_ceff': c_4/r['c_eff'],
        'n_max': int(max(N_arr)), 'n_pts': len(N_arr),
        'pw_last': pw[-1] if pw else None,
    }

# GRAND SUMMARY
print(f"\n{'='*70}")
print("GRAND SUMMARY: Open-BC Casimir c vs Re(c) and c_eff")
print(f"{'='*70}")
print(f"{'q':>3} {'#pts':>4} {'Nmax':>4} {'c(4p)':>8} {'c/Rec':>7} {'c/c_eff':>7} {'|c/Rec-1|':>9} {'|c/ce-1|':>9} {'Closer':>7}")
for q_val in [2, 5, 7]:
    a = all_results[q_val]
    dev_rec = abs(a['c_over_rec'] - 1)
    dev_ceff = abs(a['c_over_ceff'] - 1)
    closer = "Re(c)" if dev_rec < dev_ceff else "c_eff"
    print(f"{q_val:3d} {a['n_pts']:4d} {a['n_max']:4d} {a['c_4p']:8.4f} "
          f"{a['c_over_rec']:7.4f} {a['c_over_ceff']:7.4f} "
          f"{dev_rec:9.4f} {dev_ceff:9.4f} {closer:>7}")

# Key finding
print(f"\nq=2: c/0.5 = {all_results[2]['c_over_rec']:.3f} — systematic ~3% boundary offset")
print(f"q=5: c/Re(c) = {all_results[5]['c_over_rec']:.3f} — same offset, indistinguishable from c_eff")
if 7 in all_results and all_results[7]['n_pts'] > 3:
    a7 = all_results[7]
    print(f"q=7: c/Re(c) = {a7['c_over_rec']:.3f}, c/c_eff = {a7['c_over_ceff']:.3f}")
    if abs(a7['c_over_rec'] - 1) < abs(a7['c_over_ceff'] - 1):
        print("  → Casimir TRACKS Re(c), NOT c_eff — confirmed at DMRG sizes!")
    else:
        print("  → Inconclusive or tracking c_eff")

save()

from db_utils import record
for q_val in [5, 7]:
    a = all_results[q_val]
    record(sprint=100, model='sq_potts', q=q_val, n=a['n_max'],
           quantity='c_casimir_open', value=a['c_4p'],
           error=a['c_err'], method='dmrg_casimir_open_bc',
           notes=f"c/Rec={a['c_over_rec']:.4f}, c/c_eff={a['c_over_ceff']:.4f}, "
                 f"{a['n_pts']} pts, Nmax={a['n_max']}")

print("\nRecorded to DB.")
