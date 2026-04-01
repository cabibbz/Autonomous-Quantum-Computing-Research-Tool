#!/usr/bin/env python3
"""Sprint 043c: q=10 Potts n=12 crossing test.
Compare CV at n=8 (from exp_043b) vs n=12 at same g values.
Crossing = second-order. No crossing (step) = first-order.
Key g values: 0.10, 0.15, 0.20, 0.25, 0.30 (spanning transition).
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

def mi_cv(n, q, g_J, chi_max=10):
    t0 = time.time()
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g_J, 'bc_MPS': 'finite'})
    psi = MPS.from_product_state(model.lat.mps_sites(), [0]*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-5,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-8},
        'max_sweeps': 8,
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

# ── n=12 crossing test ──
q = 10
n = 12
g_values = [0.10, 0.15, 0.20, 0.25, 0.30]

# n=8 reference data from exp_043b
n8_cv = {0.10: 0.4924, 0.15: 0.6122, 0.20: 0.6412, 0.25: 0.6584, 0.30: 0.6630}

results_all = []
print(f"=== q={q} Potts n={n} crossing test ===")
print(f"n=8 reference: {n8_cv}")
t_start = time.time()

for g in g_values:
    t0 = time.time()
    try:
        cv, E0, mi_vals, t_dmrg, t_mi = mi_cv(n, q, g, chi_max=10)
        dt = time.time() - t0
        n8 = n8_cv[g]
        delta = cv - n8
        direction = "↑" if delta > 0 else "↓"
        print(f"  g={g:.2f}: n12_CV={cv:.4f}, n8_CV={n8:.4f}, delta={delta:+.4f} {direction}, "
              f"DMRG={t_dmrg:.1f}s, MI={t_mi:.1f}s")
        results_all.append({
            'g': g, 'n12_cv': cv, 'n8_cv': n8, 'delta': delta,
            'E0': E0, 'mi_mean': float(np.mean(mi_vals)),
            'time_dmrg': t_dmrg, 'time_mi': t_mi, 'status': 'ok'
        })
    except Exception as e:
        dt = time.time() - t0
        print(f"  g={g:.2f}: FAILED ({e.__class__.__name__}), time={dt:.1f}s")
        results_all.append({'g': g, 'status': 'failed', 'error': str(e)[:100], 'time_s': dt})

    with open('results/sprint_043c_q10_n12_cross.json', 'w') as f:
        json.dump({'sprint': '043c', 'q': q, 'n': n, 'chi_max': 10,
                   'n8_reference': n8_cv, 'data': results_all,
                   'total_time': time.time()-t_start}, f, indent=2)

t_total = time.time() - t_start
print(f"\nTotal: {t_total:.1f}s")

# Crossing analysis
print("\n=== CROSSING ANALYSIS ===")
ok_results = [r for r in results_all if r.get('status') == 'ok']
if len(ok_results) >= 2:
    has_decrease = any(r['delta'] < -0.01 for r in ok_results)
    has_increase = any(r['delta'] > 0.01 for r in ok_results)
    if has_decrease and has_increase:
        print("CROSSING DETECTED → Second-order transition")
        # Find crossing point by interpolation
        for i in range(len(ok_results) - 1):
            r1, r2 = ok_results[i], ok_results[i+1]
            if r1['delta'] * r2['delta'] < 0:
                g1, g2 = r1['g'], r2['g']
                d1, d2 = r1['delta'], r2['delta']
                g_cross = g1 + (g2 - g1) * (-d1) / (d2 - d1)
                print(f"  Crossing at g_c ≈ {g_cross:.3f} (between g={g1} and g={g2})")
    elif has_decrease and not has_increase:
        print("ALL DECREASE → Ordered phase only (need higher g)")
    elif has_increase and not has_decrease:
        print("ALL INCREASE → Possible first-order (no crossing)")
    else:
        print("FLAT → Insufficient resolution")
