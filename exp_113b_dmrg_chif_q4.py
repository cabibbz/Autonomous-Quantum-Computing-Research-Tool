"""Sprint 113b: DMRG chi_F for q=4 BKT at n=6-20.

Key question: q=4 exact diag (periodic BC) gives alpha=1.77 (n=4-11).
BKT prediction: nu=2/3 -> alpha=1.5 with log corrections.
Does alpha converge downward toward 1.5 at larger sizes?

Also measure q=5 for comparison (walking: alpha=2.09, zero corrections).

Method: DMRG overlap chi_F (open BC), validated in 113a.
"""
import numpy as np
import json, time, os
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '113b_dmrg_chif_q4',
    'sprint': 113,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q4_exact': {},
    'q4_dmrg': {},
    'q5_exact': {},
    'q5_dmrg': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_113b_dmrg_chif_q4.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)

# --- Exact diag (open BC) ---
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

def exact_chif_overlap(n, q, g_c, dg=1e-3):
    H0 = build_sq_potts_open(n, q, g_c)
    H1 = build_sq_potts_open(n, q, g_c + dg)
    _, evecs0 = eigsh(H0, k=1, which='SA')
    _, evecs1 = eigsh(H1, k=1, which='SA')
    psi0 = evecs0[:, 0]
    psi1 = evecs1[:, 0]
    overlap = abs(np.dot(psi0, psi1))**2
    chi_F_total = (1.0 - overlap) / dg**2
    return chi_F_total / n, overlap

# --- DMRG setup ---
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
        return SqPottsSite(model_params.get('q', 4))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 0.25)
        q_val = model_params.get('q', 4)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')

def dmrg_ground_state(n, q_val, g, chi_max=60, max_sweeps=20, seed=42):
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

def dmrg_chif(n, q_val, g_c, dg=1e-3, chi_max=60):
    E0, psi0 = dmrg_ground_state(n, q_val, g_c, chi_max=chi_max)
    E1, psi1 = dmrg_ground_state(n, q_val, g_c + dg, chi_max=chi_max)
    overlap_val = abs(psi0.overlap(psi1))**2
    chi_F_total = (1.0 - overlap_val) / dg**2
    return chi_F_total / n, overlap_val, E0, max(psi0.chi)

# ============================================================
print("=" * 70)
print("Sprint 113b: DMRG chi_F for q=4 (BKT) and q=5 (walking)")
print("=" * 70)

dg = 1e-3

# === q=4 ===
q = 4
g_c = 0.25
print(f"\n{'='*60}")
print(f"q={q}, g_c={g_c}")
print(f"Prior (periodic): alpha=1.77 (n=4-11). BKT prediction: 1.5")
print(f"{'='*60}")

# Exact diag (open BC) at small sizes
print("\n  Exact diag (open BC):")
exact_sizes_q4 = [6, 7, 8, 9]
for n in exact_sizes_q4:
    dim = q**n
    if dim > 500000:
        print(f"  n={n}: dim={dim:,} too large for exact diag, skipping")
        break
    t0 = time.time()
    chi_F, overlap = exact_chif_overlap(n, q, g_c, dg=dg)
    dt = time.time() - t0
    print(f"    n={n:2d} (dim={dim:6d}): chi_F={chi_F:.6f}, F={overlap:.12f}, t={dt:.1f}s")
    results['q4_exact'][str(n)] = {
        'n': n, 'chi_F': float(chi_F), 'overlap': float(overlap), 'time_s': round(dt, 1),
    }
    save()

# DMRG
print("\n  DMRG (open BC):")
dmrg_schedule_q4 = [
    (6, 40), (8, 50), (10, 60), (12, 60),
    (14, 80), (16, 80), (18, 100), (20, 100),
]

for n, chi_max in dmrg_schedule_q4:
    t0 = time.time()
    chi_F, overlap, E0, chi_used = dmrg_chif(n, q, g_c, dg=dg, chi_max=chi_max)
    dt = time.time() - t0

    exact_ref = results['q4_exact'].get(str(n))
    if exact_ref:
        ratio = chi_F / exact_ref['chi_F']
        print(f"    n={n:2d} chi={chi_max:3d}: chi_F={chi_F:.6f}, "
              f"exact_ratio={ratio:.6f}, t={dt:.1f}s")
    else:
        print(f"    n={n:2d} chi={chi_max:3d}: chi_F={chi_F:.6f}, F={overlap:.12f}, "
              f"chi_used={chi_used}, t={dt:.1f}s")

    results['q4_dmrg'][str(n)] = {
        'n': n, 'chi_F': float(chi_F), 'overlap': float(overlap),
        'E0': float(E0), 'chi_max': chi_max, 'chi_used': chi_used,
        'time_s': round(dt, 1),
    }
    save()

    record(sprint=113, model='sq_potts', q=q, n=n,
           quantity='chi_F_open', value=float(chi_F), method='dmrg_overlap')

    if dt > 120:
        print(f"    WARNING: {dt:.0f}s, stopping q=4 DMRG to save time for q=5")
        break

# === q=5 (comparison) ===
q = 5
g_c = 0.2
print(f"\n{'='*60}")
print(f"q={q}, g_c={g_c}")
print(f"Prior (periodic): alpha=2.09. Walking: zero corrections expected.")
print(f"{'='*60}")

# Exact diag (open BC) — small sizes only
print("\n  Exact diag (open BC):")
for n in [6, 7, 8]:
    dim = q**n
    if dim > 500000:
        print(f"  n={n}: dim={dim:,} too large, skipping")
        break
    t0 = time.time()
    chi_F, overlap = exact_chif_overlap(n, q, g_c, dg=dg)
    dt = time.time() - t0
    print(f"    n={n:2d} (dim={dim:6d}): chi_F={chi_F:.6f}, F={overlap:.12f}, t={dt:.1f}s")
    results['q5_exact'][str(n)] = {
        'n': n, 'chi_F': float(chi_F), 'overlap': float(overlap), 'time_s': round(dt, 1),
    }
    save()

# DMRG q=5
print("\n  DMRG (open BC):")
dmrg_schedule_q5 = [
    (6, 40), (8, 60), (10, 80), (12, 80),
    (14, 100), (16, 100),
]

for n, chi_max in dmrg_schedule_q5:
    t0 = time.time()
    chi_F, overlap, E0, chi_used = dmrg_chif(n, q, g_c, dg=dg, chi_max=chi_max)
    dt = time.time() - t0

    exact_ref = results['q5_exact'].get(str(n))
    if exact_ref:
        ratio = chi_F / exact_ref['chi_F']
        print(f"    n={n:2d} chi={chi_max:3d}: chi_F={chi_F:.6f}, "
              f"exact_ratio={ratio:.6f}, t={dt:.1f}s")
    else:
        print(f"    n={n:2d} chi={chi_max:3d}: chi_F={chi_F:.6f}, F={overlap:.12f}, "
              f"chi_used={chi_used}, t={dt:.1f}s")

    results['q5_dmrg'][str(n)] = {
        'n': n, 'chi_F': float(chi_F), 'overlap': float(overlap),
        'E0': float(E0), 'chi_max': chi_max, 'chi_used': chi_used,
        'time_s': round(dt, 1),
    }
    save()

    record(sprint=113, model='sq_potts', q=q, n=n,
           quantity='chi_F_open', value=float(chi_F), method='dmrg_overlap')

    if dt > 120:
        print(f"    WARNING: {dt:.0f}s, stopping q=5 DMRG")
        break

# === Pairwise analysis ===
print(f"\n{'='*70}")
print("PAIRWISE ALPHA ANALYSIS")
print("=" * 70)

for qlabel, data_dict in [("q=4 exact", results['q4_exact']),
                            ("q=4 DMRG", results['q4_dmrg']),
                            ("q=5 exact", results['q5_exact']),
                            ("q=5 DMRG", results['q5_dmrg'])]:
    Ns = sorted([int(k) for k in data_dict.keys()])
    if len(Ns) < 2:
        print(f"\n  {qlabel}: insufficient data ({len(Ns)} points)")
        continue
    log_N = np.log(np.array(Ns, dtype=float))
    log_chi = np.log(np.array([data_dict[str(n)]['chi_F'] for n in Ns]))

    print(f"\n  {qlabel}: sizes = {Ns}")
    pw_alpha = []
    for i in range(len(Ns) - 1):
        alpha = (log_chi[i+1] - log_chi[i]) / (log_N[i+1] - log_N[i])
        pw_alpha.append(alpha)
        print(f"    ({Ns[i]:2d},{Ns[i+1]:2d}): alpha = {alpha:.4f}")

    p = np.polyfit(log_N, log_chi, 1)
    alpha_global = p[0]
    print(f"    Global alpha ({len(Ns)} pts): {alpha_global:.4f}")

    prefix = qlabel.replace(' ', '_').replace('=', '')
    results[f'{prefix}_pairwise'] = {
        'sizes': Ns, 'pw_alpha': [float(a) for a in pw_alpha],
        'alpha_global': float(alpha_global),
    }

save()
print(f"\nResults saved.")
