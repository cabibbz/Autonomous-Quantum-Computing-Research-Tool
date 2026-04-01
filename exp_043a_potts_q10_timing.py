#!/usr/bin/env python3
"""Sprint 043a: q=10 Potts timing test at n=8, single g-point.
Uses DIRECT MPS contraction for rho_ij instead of Gell-Mann correlation reconstruction.
For d=10, Gell-Mann needs 99^2=9801 correlation calls — too slow.
Direct contraction is O(n * chi^2 * d) per pair.
"""
import numpy as np
import numpy.linalg as la
import json, time

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

t_start = time.time()

# ── PottsChain model (Kronecker-delta coupling, any q) ──
class PottsSite(Site):
    def __init__(self, q):
        leg = npc.LegCharge.from_trivial(q)
        Site.__init__(self, leg, [str(a) for a in range(q)], sort_charge=False)
        for a in range(q):
            P = np.zeros((q, q), dtype=complex)
            P[a, a] = 1.0
            self.add_op(f'P{a}', P)
        X = np.zeros((q, q), dtype=complex)
        for a in range(q):
            X[(a + 1) % q, a] = 1.0
        self.add_op('X', X, hc='Xhc')
        self.add_op('Xhc', X.conj().T, hc='X')
        self.add_op('Xphc', X + X.conj().T, hc='Xphc')

class PottsModel(CouplingMPOModel):
    def init_sites(self, model_params):
        return PottsSite(model_params.get('q', 10))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 1.0)
        q = model_params.get('q', 10)
        for a in range(q):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'Xphc')

class PottsChain(PottsModel, NearestNeighborModel):
    pass

# ── Direct MPS contraction for two-site rho ──
def get_mps_tensors(psi, n):
    """Extract MPS tensors as numpy arrays. Shape: (chi_L, d, chi_R)."""
    tensors = []
    for k in range(n):
        B = psi.get_B(k, 'B')
        tensors.append(B.to_ndarray())
    return tensors

def compute_environments(tensors, n):
    """Precompute left and right environments."""
    # Left environments: L[k] has shape (chi_k, chi_k)
    left_envs = [None] * (n + 1)
    chi_0 = tensors[0].shape[0]
    left_envs[0] = np.eye(chi_0, dtype=complex)
    for k in range(n):
        A = tensors[k]
        L = left_envs[k]
        left_envs[k+1] = np.einsum('ac,axb,cyd->bd', L, A, A.conj())

    # Right environments: R[k] has shape (chi_k, chi_k)
    right_envs = [None] * (n + 1)
    chi_n = tensors[n-1].shape[2]
    right_envs[n] = np.eye(chi_n, dtype=complex)
    for k in range(n-1, -1, -1):
        A = tensors[k]
        R = right_envs[k+1]
        right_envs[k] = np.einsum('axb,bd,cyd->ac', A, R, A.conj())

    return left_envs, right_envs

def compute_rho_ij(tensors, left_envs, right_envs, i, j):
    """Compute two-site reduced density matrix rho_ij via MPS contraction.
    Returns rho as (d^2, d^2) matrix.
    """
    d = tensors[i].shape[1]
    L = left_envs[i]
    R = right_envs[j+1]
    A_i = tensors[i]
    A_j = tensors[j]

    # After opening site i: env[si, si', beta_ket, beta_bra]
    env = np.einsum('ac,axb,cyd->xybd', L, A_i, A_i.conj())

    # Transfer through intermediate sites (trace over physical)
    for k in range(i+1, j):
        A_k = tensors[k]
        # Two-step contraction to avoid huge intermediate
        tmp = np.einsum('xybd,bsc->xydsc', env, A_k)
        env = np.einsum('xydsc,dse->xyce', tmp, A_k.conj())

    # Open site j and contract with right environment
    # tmp[x=si, y=si', d=bra_bond_Lj, s=sj, c=ket_bond_Rj]
    tmp = np.einsum('xybd,bsc->xydsc', env, A_j)
    # tmp2[x=si, y=si', s=sj, c=ket_Rj, t=sj', f=bra_Rj]
    tmp2 = np.einsum('xydsc,dtf->xysctf', tmp, A_j.conj())
    # Contract with R[c=ket_Rj, f=bra_Rj]
    rho = np.einsum('xysctf,cf->xyst', tmp2, R)
    # rho[si, si', sj, sj'] — reshape to (d^2, d^2)
    rho = rho.reshape(d*d, d*d)
    rho = (rho + rho.conj().T) / 2
    # Fix negative eigenvalues
    ev = la.eigvalsh(rho)
    if np.min(ev) < -1e-10:
        ev_pos = np.maximum(ev, 0)
        evec = la.eigh(rho)[1]
        rho = evec @ np.diag(ev_pos) @ evec.conj().T
    rho /= np.trace(rho)
    return rho

def entropy(rho):
    ev = la.eigvalsh(rho)
    ev = ev[ev > 1e-15]
    return float(-np.sum(ev * np.log2(ev)))

def mi_cv_potts_direct(n, q, g_J, chi_max=20):
    """Compute MI-CV using direct MPS contraction (fast for large d)."""
    t0 = time.time()
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g_J,
                         'bc_MPS': 'finite'})
    psi = MPS.from_product_state(model.lat.mps_sites(), [0]*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-6,
        'trunc_params': {'chi_max': chi_max}, 'max_sweeps': 6,
    })
    E0, _ = eng.run()
    t_dmrg = time.time() - t0
    print(f"  DMRG: E0={E0:.6f}, time={t_dmrg:.1f}s")

    t1 = time.time()
    tensors = get_mps_tensors(psi, n)
    left_envs, right_envs = compute_environments(tensors, n)

    # All-pairs MI
    mi_vals = []
    d = q
    for i in range(n):
        for j in range(i+1, n):
            rho_ij = compute_rho_ij(tensors, left_envs, right_envs, i, j)
            rho_d = rho_ij.reshape(d, d, d, d)
            rho_i = np.trace(rho_d, axis1=1, axis2=3)
            rho_j = np.trace(rho_d, axis1=0, axis2=2)
            S_i = entropy(rho_i)
            S_j = entropy(rho_j)
            S_ij = entropy(rho_ij)
            mi = S_i + S_j - S_ij
            mi_vals.append(max(float(np.real(mi)), 0))

    t_mi = time.time() - t1
    print(f"  MI computation: {len(mi_vals)} pairs, time={t_mi:.1f}s")

    mi_pos = [m for m in mi_vals if m > 1e-10]
    if len(mi_pos) < 2:
        return 0.0, E0, mi_vals
    cv = float(np.std(mi_pos) / np.mean(mi_pos))
    return cv, E0, mi_vals

# ── Run timing test ──
q = 10
n = 8
g = 0.3  # Ordered phase
print(f"=== q={q} Potts timing test: n={n}, g={g} ===")
print(f"  Direct MPS contraction (no Gell-Mann basis needed)")

t_model = time.time()
# Quick model build test
_m = PottsChain({'L': 2, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
print(f"  Model build (L=2): {time.time()-t_model:.1f}s")

cv, E0, mi_vals = mi_cv_potts_direct(n, q, g, chi_max=10)
t_total = time.time() - t_start

print(f"\nResult: CV={cv:.4f}, E0={E0:.6f}")
print(f"MI: mean={np.mean(mi_vals):.6f}, std={np.std(mi_vals):.6f}, max={np.max(mi_vals):.6f}")
print(f"Total time: {t_total:.1f}s")

results = {
    'sprint': '043a', 'experiment': 'q10_potts_timing',
    'q': q, 'n': n, 'g': g, 'chi_max': 20,
    'cv': cv, 'E0': E0,
    'mi_mean': float(np.mean(mi_vals)),
    'mi_std': float(np.std(mi_vals)),
    'mi_max': float(np.max(mi_vals)),
    'total_time_s': t_total,
    'method': 'direct_mps_contraction',
}
with open('results/sprint_043a_timing.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"Saved to results/sprint_043a_timing.json")
