#!/usr/bin/env python3
"""Sprint 095b2: BW R² via DMRG for q=5 at large n, open BC.

CRITICAL: chi_max must be >= q^nA at the subsystem cut, otherwise rho_A
is truncated and BW comparison is meaningless. For nA=3 q=5: chi>=125.
For nA=4 q=5: chi>=625 (expensive).

Also fixes basis ordering: MPS contraction gives big-endian index,
but our Hamiltonian operators use little-endian encoding.

Uses open-BC BW envelope via doubling trick (N_eff = 2n).
"""
import numpy as np
import json, time
from db_utils import record

t0 = time.time()
q = 5
gc = 1.0 / q

results = {
    'experiment': '095b2_bw_dmrg_q5',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q, 'gc': gc,
    'data': {},
}

def save():
    with open("results/sprint_095b2_bw_dmrg_q5.json", "w") as f:
        json.dump(results, f, indent=2, default=str)


# --- TeNPy DMRG setup ---
from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

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
        return SqPottsSite(model_params.get('q', q))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', gc)
        q_val = model_params.get('q', q)
        for a in range(q_val):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'SqField')


def build_basis_permutation(nA, q_val):
    """Build permutation: perm[little_endian_idx] = big_endian_idx.

    MPS contraction gives big-endian ordering (site 0 = most significant digit).
    Our operators use little-endian: idx = c[0] + c[1]*q + c[2]*q^2 + ...
    """
    dim = q_val**nA
    perm = np.zeros(dim, dtype=int)
    for idx in range(dim):
        # Decode little-endian
        c = []
        tmp = idx
        for _ in range(nA):
            c.append(tmp % q_val)
            tmp //= q_val
        # Encode big-endian: c[0]*q^{nA-1} + c[1]*q^{nA-2} + ... + c[nA-1]
        be_idx = 0
        for i in range(nA):
            be_idx = be_idx * q_val + c[i]
        perm[idx] = be_idx
    return perm


def extract_rho_A_from_mps(psi, nA, q_val):
    """Extract full rho_A for first nA sites, in little-endian basis."""
    # Contract MPS tensors for sites 0..nA-1
    B0 = psi.get_B(0, 'A')
    B0_np = B0.to_ndarray()  # (1, q, chi_1)
    tensor = B0_np.reshape(q_val, -1)  # (q, chi_1)

    for i in range(1, nA):
        Bi = psi.get_B(i, 'A')
        Bi_np = Bi.to_ndarray()  # (chi_i, q, chi_{i+1})
        tensor = np.tensordot(tensor, Bi_np, axes=([-1], [0]))
        d_phys = tensor.shape[0] * tensor.shape[1]
        tensor = tensor.reshape(d_phys, -1)

    # tensor: (q^nA, chi_nA) in big-endian ordering
    # Permute rows to little-endian ordering (matching our operator basis)
    perm = build_basis_permutation(nA, q_val)
    tensor_le = tensor[perm, :]

    # rho_A = tensor @ tensor^dagger
    rho_A = tensor_le @ tensor_le.conj().T
    tr = np.trace(rho_A).real
    if abs(tr) > 1e-10:
        rho_A = rho_A / tr

    return np.real(rho_A)


def entanglement_hamiltonian(rho_A):
    evals, evecs = np.linalg.eigh(rho_A)
    evals = np.clip(evals, 1e-30, None)
    H_E = -evecs @ np.diag(np.log(evals)) @ evecs.conj().T
    return np.real(H_E)


def potts_delta_bond(q, nA, si, sj):
    dimA = q**nA
    op = np.zeros((dimA, dimA))
    for idx in range(dimA):
        c = []
        tmp = idx
        for _ in range(nA):
            c.append(tmp % q)
            tmp //= q
        if c[si] == c[sj]:
            op[idx, idx] = 1.0
    return op


def clock_field_site(q, nA, site):
    dimA = q**nA
    op = np.zeros((dimA, dimA))
    for idx in range(dimA):
        c = []
        tmp = idx
        for _ in range(nA):
            c.append(tmp % q)
            tmp //= q
        for k in range(1, q):
            c2 = c.copy()
            c2[site] = (c[site] + k) % q
            idx2 = 0
            for i in range(nA - 1, -1, -1):
                idx2 = idx2 * q + c2[i]
            op[idx, idx2] += 1.0
    return op


def bw_envelope_bond_open(n, nA, bond_x):
    """BW envelope for open BC, subsystem at left edge. Doubling trick: N_eff=2n."""
    N_eff = 2 * n
    return (N_eff / np.pi) * np.sin(np.pi * (bond_x + 1) / N_eff) * \
           np.sin(np.pi * (nA - bond_x - 1) / N_eff) / np.sin(np.pi * nA / N_eff)


def bw_envelope_site_open(n, nA, site_x):
    N_eff = 2 * n
    return (N_eff / np.pi) * np.sin(np.pi * (site_x + 0.5) / N_eff) * \
           np.sin(np.pi * (nA - site_x - 0.5) / N_eff) / np.sin(np.pi * nA / N_eff)


def build_H_BW_open(q, n, nA, g):
    dimA = q**nA
    H_BW = np.zeros((dimA, dimA))
    for x in range(nA - 1):
        beta = bw_envelope_bond_open(n, nA, x)
        delta = potts_delta_bond(q, nA, x, x + 1)
        H_BW -= beta * delta
    for x in range(nA):
        beta = bw_envelope_site_open(n, nA, x)
        field = clock_field_site(q, nA, x)
        H_BW -= g * beta * field
    return H_BW


def frobenius_R2(H_target, H_model):
    dim = H_target.shape[0]
    t_mean = np.trace(H_target) / dim
    m_mean = np.trace(H_model) / dim
    H_t = H_target - t_mean * np.eye(dim)
    H_m = H_model - m_mean * np.eye(dim)
    alpha = np.sum(H_t * H_m) / np.sum(H_m * H_m)
    residual = H_t - alpha * H_m
    ss_res = np.sum(residual**2)
    ss_tot = np.sum(H_t**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return R2, alpha


# DMRG with chi >= q^nA at the subsystem cut.
# nA=3: need chi>=125, use 150. nA=4: need chi>=625, try for small n only.
runs_nA3 = [
    (10, 150),   # n=10, quick
    (12, 150),
    (16, 150),
    (20, 180),
    (24, 200),
]

print(f"Sprint 095b2: DMRG BW R² for q={q}, open BC")
print(f"g_c = 1/{q} = {gc:.4f}")
print(f"chi_max >= q^nA = {q**3} for nA=3")
print("=" * 70, flush=True)

# --- nA=3 runs ---
nA = 3
results['data'][f'nA{nA}'] = {}

for n, chi_max in runs_nA3:
    ratio = nA / n
    key = f"n{n}"
    print(f"\n{'='*60}")
    print(f"n={n}, chi_max={chi_max}, nA={nA}, nA/n={ratio:.3f}", flush=True)
    t1 = time.time()

    model = SqPottsChain({'L': n, 'q': q, 'J': 1.0, 'g': gc, 'bc_MPS': 'finite'})
    np.random.seed(42 + n)
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 40,
    })
    E0, _ = eng.run()
    chi_used = max(psi.chi)
    chi_at_cut = psi.chi[nA - 1] if nA - 1 < len(psi.chi) else chi_used
    dt_dmrg = time.time() - t1
    print(f"  DMRG: {dt_dmrg:.1f}s, E0={E0:.8f}, chi_max_used={chi_used}, chi_at_cut={chi_at_cut}", flush=True)

    if dt_dmrg > 250:
        print(f"  SKIP remaining: DMRG too slow")
        break

    # Extract rho_A
    rho_A = extract_rho_A_from_mps(psi, nA, q)
    dimA = q**nA

    # Check entropy
    evals_rho = np.linalg.eigvalsh(rho_A)
    evals_rho = evals_rho[evals_rho > 1e-30]
    S_vN = -np.sum(evals_rho * np.log(evals_rho))
    rank = len(evals_rho)

    # Also get S from Schmidt spectrum for comparison
    sv = psi.get_SL(nA - 1)
    S_schmidt = -np.sum(sv**2 * np.log(np.clip(sv**2, 1e-30, None)))

    print(f"  rho_A: {dimA}x{dimA}, rank={rank}, S_vN(rho)={S_vN:.4f}, S(Schmidt)={S_schmidt:.4f}")

    # Build H_E and H_BW
    H_E = entanglement_hamiltonian(rho_A)
    H_BW = build_H_BW_open(q, n, nA, gc)
    R2, alpha = frobenius_R2(H_E, H_BW)
    one_minus_R2 = 1.0 - R2

    print(f"  BW R²={R2:.8f}, 1-R²={one_minus_R2:.2e}, alpha={alpha:.4f}")

    entry = {
        'nA': nA, 'n': n, 'dimA': dimA, 'bc': 'open',
        'ratio_nA_n': float(ratio), 'chi_max': chi_max,
        'chi_at_cut': int(chi_at_cut),
        'E0': float(E0), 'S_vN': float(S_vN), 'S_schmidt': float(S_schmidt),
        'rank': rank, 'R2': float(R2), 'one_minus_R2': float(one_minus_R2),
        'alpha': float(alpha), 'time': dt_dmrg,
    }
    results['data'][f'nA{nA}'][key] = entry
    record(sprint=95, model='sq_potts', q=q, n=n,
           quantity='bw_R2_dmrg_open', value=R2,
           method=f'095b2_nA{nA}_n{n}')
    save()

# --- nA=4 runs (chi>=625 needed, only try n=10,12) ---
nA = 4
results['data'][f'nA{nA}'] = {}
dimA4 = q**4  # 625

runs_nA4 = [
    (10, 650),
    (12, 650),
]

print(f"\n\n{'='*70}")
print(f"nA=4: need chi>={dimA4}, expensive")

for n, chi_max in runs_nA4:
    remaining = 300 - (time.time() - t0)
    if remaining < 60:
        print(f"  SKIP n={n}: only {remaining:.0f}s remaining")
        continue

    ratio = nA / n
    key = f"n{n}"
    print(f"\nn={n}, chi_max={chi_max}, nA={nA}, nA/n={ratio:.3f}", flush=True)
    t1 = time.time()

    model = SqPottsChain({'L': n, 'q': q, 'J': 1.0, 'g': gc, 'bc_MPS': 'finite'})
    np.random.seed(42 + n)
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 30,
    })
    E0, _ = eng.run()
    chi_used = max(psi.chi)
    chi_at_cut = psi.chi[nA - 1] if nA - 1 < len(psi.chi) else chi_used
    dt_dmrg = time.time() - t1
    print(f"  DMRG: {dt_dmrg:.1f}s, E0={E0:.8f}, chi_used={chi_used}, chi_at_cut={chi_at_cut}")

    rho_A = extract_rho_A_from_mps(psi, nA, q)
    evals_rho = np.linalg.eigvalsh(rho_A)
    evals_rho = evals_rho[evals_rho > 1e-30]
    S_vN = -np.sum(evals_rho * np.log(evals_rho))
    rank = len(evals_rho)

    sv = psi.get_SL(nA - 1)
    S_schmidt = -np.sum(sv**2 * np.log(np.clip(sv**2, 1e-30, None)))

    print(f"  rho_A: {dimA4}x{dimA4}, rank={rank}, S_vN(rho)={S_vN:.4f}, S(Schmidt)={S_schmidt:.4f}")

    H_E = entanglement_hamiltonian(rho_A)
    H_BW = build_H_BW_open(q, n, nA, gc)
    R2, alpha = frobenius_R2(H_E, H_BW)
    one_minus_R2 = 1.0 - R2

    print(f"  BW R²={R2:.8f}, 1-R²={one_minus_R2:.2e}, alpha={alpha:.4f}")

    entry = {
        'nA': nA, 'n': n, 'dimA': dimA4, 'bc': 'open',
        'ratio_nA_n': float(ratio), 'chi_max': chi_max,
        'chi_at_cut': int(chi_at_cut),
        'E0': float(E0), 'S_vN': float(S_vN), 'S_schmidt': float(S_schmidt),
        'rank': rank, 'R2': float(R2), 'one_minus_R2': float(one_minus_R2),
        'alpha': float(alpha), 'time': dt_dmrg,
    }
    results['data'][f'nA{nA}'][key] = entry
    record(sprint=95, model='sq_potts', q=q, n=n,
           quantity='bw_R2_dmrg_open', value=R2,
           method=f'095b2_nA{nA}_n{n}')
    save()

# Summary
print(f"\n\n{'='*70}")
print("SUMMARY: DMRG open-BC BW R² for q=5")
for nA_key in ['nA3', 'nA4']:
    entries = results['data'].get(nA_key, {})
    if not entries:
        continue
    print(f"\n{nA_key}:")
    for key in sorted(entries.keys()):
        e = entries[key]
        print(f"  n={e['n']:>2}, nA/n={e['ratio_nA_n']:.3f}, chi_cut={e['chi_at_cut']}: "
              f"1-R²={e['one_minus_R2']:.2e}, S={e['S_vN']:.3f}, rank={e['rank']}")

save()
print(f"\nTotal time: {time.time()-t0:.1f}s")
print("Done!")
