#!/usr/bin/env python3
"""Sprint 056b: Measure central charge c(q=7) via exact diag and DMRG.
Exact diag: n=4 (7^4=2401), n=6 (7^6=117649).
DMRG: n=8, n=10 if timing allows (d=7, chi=20-40).
g_c(q=7) = 0.535 from Sprint 052.

Predictions from 056a:
  log: 1.35, power: 1.41, quadratic: 1.00, analytic cont: 0.50-0.70
"""
import numpy as np, json, time
from scipy.sparse import kron, eye, csr_matrix
from scipy.sparse.linalg import eigsh

def potts_hamiltonian(q, n, g, J=1.0):
    """Build Potts H = -J*sum delta(s_i,s_j) - g*sum (X + X†)."""
    d = q
    dim = d**n
    # Operators
    Id = eye(d, format='csr')
    # Projectors |a><a|
    projs = []
    for a in range(q):
        P = csr_matrix((d, d))
        P[a, a] = 1.0
        projs.append(P)
    # X: |s> -> |s+1 mod q>
    X = csr_matrix((d, d))
    for s in range(q):
        X[(s + 1) % q, s] = 1.0
    Xphc = X + X.T  # X + X†

    H = csr_matrix((dim, dim))
    for i in range(n - 1):  # nearest-neighbor
        for a in range(q):
            # -J * P_a(i) * P_a(i+1)
            op = csr_matrix(eye(1))
            for j in range(n):
                if j == i or j == i + 1:
                    op = kron(op, projs[a], format='csr')
                else:
                    op = kron(op, Id, format='csr')
            H -= J * op
    for i in range(n):  # transverse field
        op = csr_matrix(eye(1))
        for j in range(n):
            if j == i:
                op = kron(op, Xphc, format='csr')
            else:
                op = kron(op, Id, format='csr')
        H -= g * op
    return H

def entanglement_entropy(psi, n, d, l):
    """Entropy of first l sites."""
    dim_A = d**l
    dim_B = d**(n - l)
    rho_A = psi.reshape(dim_A, dim_B)
    s = np.linalg.svd(rho_A, compute_uv=False)
    s2 = s[s > 1e-15]**2
    return float(-np.sum(s2 * np.log(s2)))

def extract_c_profile(S_list, n):
    """Central charge from entropy profile S(l) vs chord distance."""
    ls = np.arange(1, n)
    l_chord = (2 * n / np.pi) * np.sin(np.pi * ls / n)
    ln_chord = np.log(l_chord)
    S = np.array(S_list)
    quarter = max(1, n // 4)
    three_quarter = min(n - 1, 3 * n // 4)
    mask = (ls >= quarter) & (ls <= three_quarter)
    if mask.sum() < 2:
        mask = np.ones(len(ls), dtype=bool)
    A = np.vstack([ln_chord[mask], np.ones(mask.sum())]).T
    slope, _ = np.linalg.lstsq(A, S[mask], rcond=None)[0]
    return float(6 * slope)

def extract_c_fss(S_half_list, n_list):
    """Pairwise central charge from S(n/2) at consecutive sizes."""
    c_pairs = []
    for i in range(len(n_list) - 1):
        n1, n2 = n_list[i], n_list[i + 1]
        S1, S2 = S_half_list[i], S_half_list[i + 1]
        c = 6 * (S2 - S1) / (np.log(n2) - np.log(n1))
        c_pairs.append(float(c))
    return c_pairs

# === Exact diag experiments ===
q = 7
g_c = 0.535
data = []

for n in [4, 6]:
    t0 = time.time()
    dim = q**n
    print(f"q={q}, n={n}: dim={dim}, building H...", flush=True)
    H = potts_hamiltonian(q, n, g_c)
    print(f"  H built ({time.time()-t0:.1f}s), finding ground state...", flush=True)
    E0, psi0 = eigsh(H, k=1, which='SA')
    E0 = float(E0[0])
    psi = psi0[:, 0]
    print(f"  E0={E0:.6f} ({time.time()-t0:.1f}s), computing entropies...", flush=True)

    S_profile = []
    for l in range(1, n):
        S_profile.append(entanglement_entropy(psi, n, q, l))

    S_half = S_profile[n // 2 - 1]
    c_prof = extract_c_profile(S_profile, n)
    dt = time.time() - t0
    print(f"  n={n}: S_half={S_half:.6f}, c_profile={c_prof:.4f}, time={dt:.1f}s", flush=True)
    data.append({
        'method': 'exact_diag', 'n': n, 'dim': dim,
        'E0': E0, 'S_half': S_half, 'S_profile': S_profile,
        'c_profile': c_prof, 'time': dt
    })
    # Save incrementally
    results = {'experiment': '056b', 'model': 'Potts', 'q': q, 'g_c': g_c,
               'predictions': {'log': 1.35, 'power': 1.41, 'quadratic': 1.00},
               'data': data}
    with open('results/exp_056b.json', 'w') as f:
        json.dump(results, f, indent=2)

# FSS from exact diag
n_list = [d['n'] for d in data]
S_list = [d['S_half'] for d in data]
c_fss = extract_c_fss(S_list, n_list)
print(f"\nFSS pairwise c: {c_fss}")

# === DMRG at n=8 ===
print(f"\n--- DMRG n=8 ---", flush=True)
try:
    import warnings; warnings.filterwarnings('ignore')
    from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
    from tenpy.networks.site import Site
    from tenpy.linalg import np_conserved as npc
    from tenpy.networks.mps import MPS
    from tenpy.algorithms import dmrg as dmrg_alg

    class PottsSite(Site):
        def __init__(self, q):
            leg = npc.LegCharge.from_trivial(q)
            Site.__init__(self, leg, [str(a) for a in range(q)], sort_charge=False)
            for a in range(q):
                P = np.zeros((q, q), dtype=complex); P[a, a] = 1.0
                self.add_op(f'P{a}', P)
            X = np.zeros((q, q), dtype=complex)
            for a in range(q):
                X[(a + 1) % q, a] = 1.0
            self.add_op('X', X, hc='Xhc')
            self.add_op('Xhc', X.conj().T, hc='X')
            self.add_op('Xphc', X + X.conj().T, hc='Xphc')

    class PottsChain(CouplingMPOModel, NearestNeighborModel):
        def init_sites(self, model_params):
            return PottsSite(model_params.get('q', 3))
        def init_terms(self, model_params):
            J = model_params.get('J', 1.0); g = model_params.get('g', 1.0); q_val = model_params.get('q', 3)
            for a in range(q_val):
                self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
            self.add_onsite(-g, 0, 'Xphc')

    for n_dmrg, chi_max in [(8, 30), (10, 30)]:
        t0 = time.time()
        model = PottsChain({'L': n_dmrg, 'q': q, 'J': 1.0, 'g': g_c, 'bc_MPS': 'finite'})
        np.random.seed(42 + n_dmrg)
        init = [np.random.randint(q) for _ in range(n_dmrg)]
        psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
        eng = dmrg_alg.TwoSiteDMRGEngine(psi, model, {
            'mixer': True, 'max_E_err': 1e-10,
            'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
            'max_sweeps': 30,
        })
        E0_dmrg, _ = eng.run()
        S_prof_dmrg = [float(s) for s in psi.entanglement_entropy()]
        S_half_dmrg = S_prof_dmrg[n_dmrg // 2 - 1]
        c_prof_dmrg = extract_c_profile(S_prof_dmrg, n_dmrg)
        chi_actual = max(psi.chi)
        dt = time.time() - t0
        print(f"  DMRG n={n_dmrg}: S_half={S_half_dmrg:.6f}, c_profile={c_prof_dmrg:.4f}, chi={chi_actual}, time={dt:.1f}s", flush=True)
        data.append({
            'method': 'dmrg', 'n': n_dmrg, 'chi_max': chi_max, 'chi_actual': chi_actual,
            'E0': float(E0_dmrg), 'S_half': S_half_dmrg, 'S_profile': S_prof_dmrg,
            'c_profile': c_prof_dmrg, 'time': dt
        })
        results['data'] = data
        with open('results/exp_056b.json', 'w') as f:
            json.dump(results, f, indent=2)
        if dt > 120:
            print(f"  Skipping larger sizes (too slow)")
            break

except Exception as e:
    print(f"  DMRG failed: {e}")

# Final FSS
n_list_all = [d['n'] for d in data]
S_list_all = [d['S_half'] for d in data]
c_fss_all = extract_c_fss(S_list_all, n_list_all)
results['c_fss_pairs'] = c_fss_all
with open('results/exp_056b.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n=== SUMMARY q={q} ===")
for d in data:
    print(f"  n={d['n']:3d} ({d['method']:10s}): S_half={d['S_half']:.6f}, c_profile={d.get('c_profile', 'N/A')}")
print(f"  FSS pairwise c: {c_fss_all}")
print(f"\nPredictions: log=1.35, power=1.41, quad=1.00, analytic_cont=0.50")
print("Saved results/exp_056b.json")
