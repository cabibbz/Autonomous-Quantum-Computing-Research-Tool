#!/usr/bin/env python3
"""Sprint 046d: q=4 Potts MI-CV at low g for n=8,12 to find crossing.

n=16 data shows CV=0.163 at g=0.40, 0.209 at g=0.50.
Need n=8 and n=12 at same points to check if curves cross.
"""
import numpy as np
import numpy.linalg as la
import json, time

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg

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

class PottsModel(CouplingMPOModel):
    def init_sites(self, model_params):
        return PottsSite(model_params.get('q', 4))
    def init_terms(self, model_params):
        J, g, q = model_params.get('J', 1.0), model_params.get('g', 1.0), model_params.get('q', 4)
        for a in range(q):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'Xphc')

class PottsChain(PottsModel, NearestNeighborModel):
    pass

def get_mps_tensors(psi, n):
    return [psi.get_B(k, 'B').to_ndarray() for k in range(n)]

def compute_environments(tensors, n):
    left_envs = [None] * (n + 1)
    left_envs[0] = np.eye(tensors[0].shape[0], dtype=complex)
    for k in range(n):
        A = tensors[k]; L = left_envs[k]
        left_envs[k+1] = np.einsum('ac,axb,cyd->bd', L, A, A.conj())
    right_envs = [None] * (n + 1)
    right_envs[n] = np.eye(tensors[n-1].shape[2], dtype=complex)
    for k in range(n-1, -1, -1):
        A = tensors[k]; R = right_envs[k+1]
        right_envs[k] = np.einsum('axb,bd,cyd->ac', A, R, A.conj())
    return left_envs, right_envs

def compute_rho_ij(tensors, left_envs, right_envs, i, j):
    d = tensors[i].shape[1]
    L, R = left_envs[i], right_envs[j+1]
    env = np.einsum('ac,axb,cyd->xybd', L, tensors[i], tensors[i].conj())
    for k in range(i+1, j):
        A_k = tensors[k]
        tmp = np.einsum('xybd,bsc->xydsc', env, A_k)
        env = np.einsum('xydsc,dse->xyce', tmp, A_k.conj())
    tmp = np.einsum('xybd,bsc->xydsc', env, tensors[j])
    tmp2 = np.einsum('xydsc,dtf->xysctf', tmp, tensors[j].conj())
    rho = np.einsum('xysctf,cf->xyst', tmp2, R).reshape(d*d, d*d)
    rho = (rho + rho.conj().T) / 2
    ev = la.eigvalsh(rho)
    if np.min(ev) < -1e-10:
        ev_pos = np.maximum(ev, 0)
        evec = la.eigh(rho)[1]
        rho = evec @ np.diag(ev_pos) @ evec.conj().T
    rho /= np.trace(rho)
    return rho

def entropy(rho):
    ev = la.eigvalsh(rho); ev = ev[ev > 1e-15]
    return float(-np.sum(ev * np.log2(ev)))

def mi_cv(n, q, g_J, chi_max=20):
    t0 = time.time()
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g_J, 'bc_MPS': 'finite'})
    psi = MPS.from_product_state(model.lat.mps_sites(), [0]*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-5,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-8},
        'max_sweeps': 12,
    })
    E0, _ = eng.run()
    t_dmrg = time.time() - t0

    t1 = time.time()
    tensors = get_mps_tensors(psi, n)
    left_envs, right_envs = compute_environments(tensors, n)
    d = q; mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            rho_ij = compute_rho_ij(tensors, left_envs, right_envs, i, j)
            rho_d = rho_ij.reshape(d, d, d, d)
            rho_i = np.trace(rho_d, axis1=1, axis2=3)
            rho_j = np.trace(rho_d, axis1=0, axis2=2)
            mi = entropy(rho_i) + entropy(rho_j) - entropy(rho_ij)
            mi_vals.append(max(float(np.real(mi)), 0))
    t_mi = time.time() - t1

    mi_pos = [m for m in mi_vals if m > 1e-10]
    cv = float(np.std(mi_pos) / np.mean(mi_pos)) if len(mi_pos) >= 2 else 0.0
    return cv, E0, mi_vals, t_dmrg, t_mi

q = 4
chi = 20
g_values = [0.20, 0.30, 0.40, 0.50, 0.60]

results = {}
print(f"=== q={q} Potts low-g extension (direct MPS, chi={chi}) ===\n")

for n in [8, 12]:
    results[n] = []
    print(f"--- n={n} ---")
    for g in g_values:
        t0 = time.time()
        try:
            cv, E0, mi_vals, t_dmrg, t_mi = mi_cv(n, q, g, chi_max=chi)
            print(f"  g={g:.2f}: CV={cv:.4f}, E0={E0:.4f}, DMRG={t_dmrg:.1f}s")
            results[n].append({'g': g, 'cv': cv, 'E0': float(E0), 'status': 'ok',
                               'time_dmrg': t_dmrg, 'time_mi': t_mi})
        except Exception as e:
            print(f"  g={g:.2f}: FAILED ({e})")
            results[n].append({'g': g, 'status': 'failed', 'error': str(e)[:100]})
    print()

# Compare with n=16 data at same g values
print("=== CROSSING ANALYSIS ===")
print(f"{'g':>6s} {'n=8':>8s} {'n=12':>8s} {'n=16':>8s} {'n12-n8':>8s} {'n16-n12':>8s}")
# n=16 low-g data from exp_046c
n16_data = {0.40: 0.1634, 0.50: 0.2085, 0.60: 0.3515}
for g in g_values:
    cv8 = next((r['cv'] for r in results[8] if r['g'] == g and r['status'] == 'ok'), None)
    cv12 = next((r['cv'] for r in results[12] if r['g'] == g and r['status'] == 'ok'), None)
    cv16 = n16_data.get(g, None)
    if cv8 is not None and cv12 is not None:
        d12_8 = cv12 - cv8
        d16_12 = (cv16 - cv12) if cv16 is not None else None
        cv16_s = f"{cv16:.4f}" if cv16 else "   --  "
        d16_s = f"{d16_12:+.4f}" if d16_12 is not None else "   --  "
        print(f"  {g:.2f}  {cv8:.4f}  {cv12:.4f}  {cv16_s}  {d12_8:+.4f}  {d16_s}")

# Save
with open('results/sprint_046d_q4_lowg.json', 'w') as f:
    json.dump({'sprint': '046d', 'q': q, 'chi_max': chi,
               'n8': results[8], 'n12': results[12]}, f, indent=2)

print("\nResults saved to results/sprint_046d_q4_lowg.json")
