"""Sprint 113a: Validate DMRG overlap chi_F method on q=2 (exact nu=1, alpha=1.0).

Method: Ground state overlaps at g_c and g_c+dg via DMRG (open BC).
chi_F_total = (1 - |<psi(g)|psi(g+dg)>|^2) / dg^2
chi_F = chi_F_total / N

Cross-validate: exact diag (open BC) at n=8-14, then DMRG extends to n=16-24.
Prior result (periodic BC, spectral): alpha = 1.012 at (n=16,18) pair.
"""
import numpy as np
import json, time, os
from scipy.sparse import coo_matrix, csr_matrix
from gpu_utils import eigsh
from db_utils import record

results = {
    'experiment': '113a_dmrg_chif_q2',
    'sprint': 113,
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'exact': {},
    'dmrg': {},
}

def save():
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'sprint_113a_dmrg_chif_q2.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2, default=str)


# --- Exact diag (open BC) ---
def build_sq_potts_open(n, q_val, g):
    """Build open-BC S_q Potts Hamiltonian."""
    dim = q_val**n
    all_idx = np.arange(dim, dtype=np.int64)
    digits = np.zeros((dim, n), dtype=np.int64)
    tmp = all_idx.copy()
    for site in range(n):
        digits[:, site] = tmp % q_val
        tmp //= q_val
    powers = q_val ** np.arange(n, dtype=np.int64)
    # Coupling: -sum delta(s_i, s_{i+1}) for open BC (n-1 bonds)
    diag_vals = np.zeros(dim, dtype=np.float64)
    for site in range(n - 1):
        diag_vals -= (digits[:, site] == digits[:, site + 1]).astype(np.float64)
    rows_list = [all_idx]; cols_list = [all_idx]; vals_list = [diag_vals]
    # Field: -g sum_i sum_{k=1}^{q-1} |shifted><original|
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
    """Compute chi_F via ground state overlap (exact diag, open BC)."""
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
        return SqPottsSite(model_params.get('q', 2))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 0.5)
        q_val = model_params.get('q', 2)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')


def dmrg_ground_state(n, q_val, g, chi_max=60, max_sweeps=20, seed=42):
    """Run DMRG and return (E0, psi)."""
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


def dmrg_chif_overlap(n, q_val, g_c, dg=1e-3, chi_max=60):
    """Compute chi_F via DMRG ground state overlap."""
    E0, psi0 = dmrg_ground_state(n, q_val, g_c, chi_max=chi_max)
    E1, psi1 = dmrg_ground_state(n, q_val, g_c + dg, chi_max=chi_max)
    overlap_val = abs(psi0.overlap(psi1))**2
    chi_F_total = (1.0 - overlap_val) / dg**2
    return chi_F_total / n, overlap_val, E0, max(psi0.chi)


# ============================================================
q = 2
g_c = 0.5
dg = 1e-3

print("=" * 70)
print(f"Sprint 113a: DMRG chi_F validation for q={q}")
print(f"  g_c = {g_c}, dg = {dg}")
print(f"  Known: alpha = 1.0 (exact nu=1)")
print("=" * 70)

# Phase 1: Exact diag (open BC) at n=8-14
print("\n--- Phase 1: Exact diag (open BC) ---")
exact_sizes = [6, 8, 10, 12, 14]

for n in exact_sizes:
    dim = q**n
    t0 = time.time()
    chi_F, overlap = exact_chif_overlap(n, q, g_c, dg=dg)
    dt = time.time() - t0
    print(f"  n={n:2d} (dim={dim:6d}): chi_F={chi_F:.6f}, F={overlap:.12f}, t={dt:.1f}s")
    results['exact'][str(n)] = {
        'n': n, 'chi_F': float(chi_F), 'overlap': float(overlap), 'time_s': round(dt, 1),
    }
    save()

# Phase 2: DMRG (open BC) at n=8-24
print("\n--- Phase 2: DMRG (open BC) ---")
dmrg_schedule = [
    (8, 40), (10, 40), (12, 50), (14, 50),
    (16, 60), (18, 60), (20, 80), (24, 80),
]

for n, chi_max in dmrg_schedule:
    t0 = time.time()
    chi_F, overlap, E0, chi_used = dmrg_chif_overlap(n, q, g_c, dg=dg, chi_max=chi_max)
    dt = time.time() - t0

    # Cross-check with exact diag where available
    exact_ref = results['exact'].get(str(n))
    if exact_ref:
        ratio = chi_F / exact_ref['chi_F']
        print(f"  n={n:2d} chi={chi_max:3d}: chi_F={chi_F:.6f}, F={overlap:.12f}, "
              f"exact_ratio={ratio:.6f}, t={dt:.1f}s")
    else:
        print(f"  n={n:2d} chi={chi_max:3d}: chi_F={chi_F:.6f}, F={overlap:.12f}, "
              f"chi_used={chi_used}, t={dt:.1f}s")

    results['dmrg'][str(n)] = {
        'n': n, 'chi_F': float(chi_F), 'overlap': float(overlap),
        'E0': float(E0), 'chi_max': chi_max, 'chi_used': chi_used,
        'time_s': round(dt, 1),
    }
    save()

    record(sprint=113, model='sq_potts', q=q, n=n,
           quantity='chi_F_open', value=float(chi_F), method='dmrg_overlap')

    if dt > 120:
        print(f"  WARNING: {dt:.0f}s per point, may hit time limit. Adjusting schedule.")
        # Skip larger sizes if too slow
        break

# Phase 3: Pairwise analysis
print(f"\n{'=' * 70}")
print("PAIRWISE ALPHA ANALYSIS")
print("=" * 70)

for label, data_dict in [("exact", results['exact']), ("dmrg", results['dmrg'])]:
    Ns = sorted([int(k) for k in data_dict.keys()])
    if len(Ns) < 2:
        continue
    log_N = np.log(np.array(Ns, dtype=float))
    log_chi = np.log(np.array([data_dict[str(n)]['chi_F'] for n in Ns]))

    print(f"\n  {label}: sizes = {Ns}")
    pw_alpha = []
    for i in range(len(Ns) - 1):
        alpha = (log_chi[i+1] - log_chi[i]) / (log_N[i+1] - log_N[i])
        pw_alpha.append(alpha)
        print(f"    ({Ns[i]:2d},{Ns[i+1]:2d}): alpha = {alpha:.4f}")

    # Global fit
    p = np.polyfit(log_N, log_chi, 1)
    alpha_global = p[0]
    print(f"    Global alpha ({len(Ns)} pts): {alpha_global:.4f}")
    print(f"    Expected: 1.000 (exact nu=1)")

    results[f'{label}_pairwise'] = {
        'sizes': Ns,
        'pw_alpha': [float(a) for a in pw_alpha],
        'alpha_global': float(alpha_global),
    }

save()
print(f"\nResults saved to results/sprint_113a_dmrg_chif_q2.json")
