"""Sprint 124c: Extend finite DMRG chi_F for S_q q=4 to n=14-20.

Sprint 113b reached n=12 (chi_max=60, open BC, alpha_open=1.51).
This extends to larger sizes with higher bond dimension to check
whether open-BC alpha drifts toward periodic-BC value (1.77).

Also includes q=2 at same sizes as cross-check.
"""
import numpy as np
import json, time, os
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '124c_finite_dmrg_q4_ext',
    'sprint': 124,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q4_dmrg': {},
    'q2_dmrg': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), '..', 'results',
                           'sprint_124c_finite_dmrg_q4_ext.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

# === TeNPy DMRG ===
from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

import warnings
warnings.filterwarnings('ignore')

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
        return SqPottsSite(model_params.get('q', 4))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 0.25)
        q_val = model_params.get('q', 4)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')

def dmrg_ground_state(n, q_val, g, chi_max=80, max_sweeps=30, seed=42):
    model = SqPottsChain({'L': n, 'q': q_val, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(seed + n + int(g * 10000))
    init = [np.random.randint(q_val) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-14},
        'max_sweeps': max_sweeps,
    })
    E0, _ = eng.run()
    return float(E0), psi

def dmrg_chif(n, q_val, g_c, dg=1e-3, chi_max=80):
    """Compute chi_F via DMRG ground state overlap."""
    t0 = time.time()
    E0, psi0 = dmrg_ground_state(n, q_val, g_c, chi_max=chi_max)
    E1, psi1 = dmrg_ground_state(n, q_val, g_c + dg, chi_max=chi_max, seed=137)
    dt = time.time() - t0
    overlap_val = abs(psi0.overlap(psi1))**2
    chi_F_total = (1.0 - overlap_val) / dg**2
    chi_used = max(psi0.chi)
    S_mid = float(psi0.entanglement_entropy()[n//2 - 1])
    return {
        'chi_F': float(chi_F_total / n),
        'overlap': float(overlap_val),
        'E0': float(E0),
        'chi_max': chi_max,
        'chi_used': int(chi_used),
        'S_mid': float(S_mid),
        'time_s': round(dt, 1),
    }

# ============================================================
print("=" * 70)
print("Sprint 124c: Finite DMRG chi_F extension for q=4 and q=2")
print("=" * 70)

dg = 1e-3
total_time = 0

# Prior data from Sprint 113b (open BC):
# q=4: n=6→4.17, n=8→6.43, n=10→9.00, n=12→11.85
# q=2: n=8→0.244, n=10→0.254, ..., n=20→0.283
prior_q4 = {6: 4.170, 8: 6.433, 10: 9.000, 12: 11.851}

# === q=2 first (faster, validates method) ===
q = 2; g_c = 0.5
print(f"\n{'='*60}")
print(f"q={q}, g_c={g_c} — cross-check extension")
print(f"{'='*60}")

q2_schedule = [(14, 60), (16, 60), (20, 80), (24, 80)]
for n, chi_max in q2_schedule:
    if total_time > 250:
        print(f"  Time budget exceeded ({total_time:.0f}s), stopping q=2")
        break
    r = dmrg_chif(n, q, g_c, dg=dg, chi_max=chi_max)
    total_time += r['time_s']
    print(f"  n={n:2d} chi={chi_max:3d}: chi_F={r['chi_F']:.6f}, "
          f"chi_used={r['chi_used']}, S={r['S_mid']:.4f}, t={r['time_s']:.0f}s")
    results['q2_dmrg'][str(n)] = r
    save()
    record(sprint=124, model='sq_potts', q=q, n=n,
           quantity='chi_F_open', value=r['chi_F'], method='dmrg_overlap')

# === q=4 main target ===
q = 4; g_c = 0.25
print(f"\n{'='*60}")
print(f"q={q}, g_c={g_c} — extending to n=14-20")
print(f"Prior (open BC): n=6->4.17, n=8->6.43, n=10->9.00, n=12->11.85")
print(f"Open-BC alpha=1.51 (Sprint 113b)")
print(f"{'='*60}")

# Progressive: start small to estimate timing
q4_schedule = [(14, 80), (16, 100), (18, 120), (20, 150)]
for n, chi_max in q4_schedule:
    if total_time > 400:
        print(f"  Time budget exceeded ({total_time:.0f}s), stopping q=4")
        break
    print(f"  Starting n={n}, chi_max={chi_max}...")
    r = dmrg_chif(n, q, g_c, dg=dg, chi_max=chi_max)
    total_time += r['time_s']
    print(f"  n={n:2d} chi={chi_max:3d}: chi_F={r['chi_F']:.6f}, "
          f"chi_used={r['chi_used']}, S={r['S_mid']:.4f}, t={r['time_s']:.0f}s")
    results['q4_dmrg'][str(n)] = r
    save()
    record(sprint=124, model='sq_potts', q=q, n=n,
           quantity='chi_F_open', value=r['chi_F'], method='dmrg_overlap')

# === Combined analysis with Sprint 113b data ===
print(f"\n{'='*70}")
print("COMBINED ANALYSIS (Sprint 113b + 124c)")
print("=" * 70)

# q=4 combined
all_q4 = dict(prior_q4)
for k, v in results['q4_dmrg'].items():
    all_q4[int(k)] = v['chi_F']

Ns = sorted(all_q4.keys())
chi_Fs = np.array([all_q4[n] for n in Ns])
log_N = np.log(np.array(Ns, dtype=float))
log_chi = np.log(chi_Fs)

print(f"\n  q=4 (open BC): {len(Ns)} sizes, n={Ns}")
for i in range(len(Ns) - 1):
    alpha = (log_chi[i+1] - log_chi[i]) / (log_N[i+1] - log_N[i])
    print(f"    ({Ns[i]:2d},{Ns[i+1]:2d}): alpha = {alpha:.4f}")

if len(Ns) >= 3:
    p = np.polyfit(log_N, log_chi, 1)
    print(f"    Global alpha: {p[0]:.4f}")

    # Fit: chi_F = A * N^alpha * (1 + B/N^2)
    from scipy.optimize import curve_fit
    def power_corr(N, A, alpha, B):
        return A * N**alpha * (1 + B / N**2)
    try:
        popt, pcov = curve_fit(power_corr, np.array(Ns, dtype=float), chi_Fs,
                               p0=[0.5, 1.5, -1])
        perr = np.sqrt(np.diag(pcov))
        print(f"    Power+1/N^2: alpha={popt[1]:.4f}±{perr[1]:.4f}, B={popt[2]:.2f}")
        pred = power_corr(np.array(Ns, dtype=float), *popt)
        ss_res = np.sum((chi_Fs - pred)**2)
        ss_tot = np.sum((chi_Fs - np.mean(chi_Fs))**2)
        print(f"    R^2 = {1 - ss_res/ss_tot:.8f}")
    except Exception as e:
        print(f"    Curve fit failed: {e}")

    # Log-corrected: chi_F = A * N^2 * (ln N)^{-p}
    try:
        def log_corr(N, A, p_lc):
            return A * N**2 * np.log(N)**(-p_lc)
        popt_lc, pcov_lc = curve_fit(log_corr, np.array(Ns, dtype=float), chi_Fs,
                                      p0=[0.1, 1.5])
        pred_lc = log_corr(np.array(Ns, dtype=float), *popt_lc)
        ss_res_lc = np.sum((chi_Fs - pred_lc)**2)
        r2_lc = 1 - ss_res_lc / ss_tot
        print(f"    Log-corrected (alpha=2): p={popt_lc[1]:.4f}, R^2={r2_lc:.8f}")
    except Exception as e:
        print(f"    Log-corrected fit failed: {e}")

results['q4_combined'] = {
    'sizes': Ns,
    'chi_Fs': [float(c) for c in chi_Fs],
}

save()
print(f"\nTotal wall time: {total_time:.0f}s")
print("All results saved.")
