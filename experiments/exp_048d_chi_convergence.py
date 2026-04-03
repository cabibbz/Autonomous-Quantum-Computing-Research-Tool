#!/usr/bin/env python3
"""Sprint 048d: Chi convergence test for q=10 MI-CV.

Test whether absence of crossings is real or chi=20 artifact.
Run single g value at chi=20,30,40 for n=8 and n=12.
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
        return PottsSite(model_params.get('q', 10))
    def init_terms(self, model_params):
        J, g, q = model_params.get('J', 1.0), model_params.get('g', 1.0), model_params.get('q', 10)
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
        'max_sweeps': 20,
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

q = 10
# Test at 2 g values: one in ordered phase, one in disordered
# g=0.20 (ordered) and g=0.30 (disordered)
results = []

for g in [0.20, 0.30]:
    for n in [8, 12]:
        for chi in [20, 30, 40]:
            t0 = time.time()
            try:
                cv, E0, mi_vals, t_dmrg, t_mi = mi_cv(n, q, g, chi_max=chi)
                dt = time.time() - t0
                entry = {'g': g, 'n': n, 'chi': chi, 'cv': cv, 'E0': float(E0),
                         'mi_mean': float(np.mean(mi_vals)), 'time': dt, 'status': 'ok'}
                print(f"  g={g:.2f} n={n} chi={chi}: CV={cv:.4f}, E0={E0:.6f}, time={dt:.0f}s")
            except Exception as e:
                dt = time.time() - t0
                entry = {'g': g, 'n': n, 'chi': chi, 'status': 'failed', 'error': str(e)[:200], 'time': dt}
                print(f"  g={g:.2f} n={n} chi={chi}: FAILED ({str(e)[:60]}), time={dt:.0f}s")
            results.append(entry)

            with open('results/sprint_048d_chi_conv.json', 'w') as f:
                json.dump({'sprint': '048d', 'q': q, 'data': results}, f, indent=2)

# Summary
print("\n=== Chi convergence summary ===")
for g in [0.20, 0.30]:
    print(f"\ng={g:.2f}:")
    for chi in [20, 30, 40]:
        cv8 = [r['cv'] for r in results if r['g']==g and r['n']==8 and r['chi']==chi and r['status']=='ok']
        cv12 = [r['cv'] for r in results if r['g']==g and r['n']==12 and r['chi']==chi and r['status']=='ok']
        if cv8 and cv12:
            diff = cv12[0] - cv8[0]
            print(f"  chi={chi}: n8={cv8[0]:.4f}, n12={cv12[0]:.4f}, diff={diff:+.4f} {'CROSSING' if diff>0 else 'n12<n8'}")

print("\nSaved to results/sprint_048d_chi_conv.json")
