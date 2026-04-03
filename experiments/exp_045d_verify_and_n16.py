#!/usr/bin/env python3
"""Sprint 045d: Verify suspicious n=12 points and get n=16 data.

1. Re-run g=0.40, 0.50 at n=12 to verify Sprint 042d data
2. If time, get n=16 at key g values for 3-size collapse
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
        return PottsSite(model_params.get('q', 5))
    def init_terms(self, model_params):
        J, g, q = model_params.get('J', 1.0), model_params.get('g', 1.0), model_params.get('q', 5)
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
        'max_sweeps': 15,
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

q = 5
chi = 20
results = {}

# Part 1: Verify n=12 suspicious points
print("=== PART 1: Verify n=12 g=0.40, 0.50 ===")
print(f"Sprint 042d values: g=0.40 CV=0.516, g=0.50 CV=1.360")
for g in [0.40, 0.50]:
    t0 = time.time()
    cv, E0, mi_vals, t_dmrg, t_mi = mi_cv(12, q, g, chi_max=chi)
    print(f"  g={g:.2f}: CV={cv:.4f} (042d: {0.516 if g==0.40 else 1.360:.3f}), "
          f"E0={E0:.4f}, DMRG={t_dmrg:.1f}s")
    results[f'n12_g{g:.2f}'] = {'g': g, 'n': 12, 'cv': cv, 'E0': E0,
                                 'mi_mean': float(np.mean(mi_vals)),
                                 'time': time.time()-t0}

# Save incrementally
with open('results/sprint_045d_verify_n16.json', 'w') as f:
    json.dump({'sprint': '045d', 'q': q, 'chi': chi, 'results': results}, f, indent=2)

# Part 2: n=16 key points
print("\n=== PART 2: n=16 at key g values ===")
g_values_16 = [0.30, 0.35, 0.40, 0.45]
for g in g_values_16:
    t0 = time.time()
    try:
        cv, E0, mi_vals, t_dmrg, t_mi = mi_cv(16, q, g, chi_max=chi)
        print(f"  n=16 g={g:.2f}: CV={cv:.4f}, E0={E0:.4f}, DMRG={t_dmrg:.1f}s, MI={t_mi:.1f}s")
        results[f'n16_g{g:.2f}'] = {'g': g, 'n': 16, 'cv': cv, 'E0': E0,
                                     'mi_mean': float(np.mean(mi_vals)),
                                     'time': time.time()-t0}
    except Exception as e:
        print(f"  n=16 g={g:.2f}: FAILED ({str(e)[:80]})")
        results[f'n16_g{g:.2f}'] = {'g': g, 'n': 16, 'status': 'failed', 'error': str(e)[:100]}

    # Save after each point
    with open('results/sprint_045d_verify_n16.json', 'w') as f:
        json.dump({'sprint': '045d', 'q': q, 'chi': chi, 'results': results}, f, indent=2)

print("\n=== SUMMARY ===")
# Consolidated n=8,12,16 data at chi=20
n8 = {0.30: 0.3732, 0.35: 0.3554, 0.38: 0.3141, 0.40: 0.2755, 0.42: 0.2325,
      0.45: 0.1700, 0.50: 0.1072, 0.55: 0.1324}
n12_new = {0.30: 0.1342, 0.35: 0.3671, 0.38: 0.3399, 0.42: 0.2456, 0.45: 0.1677}
# Add verified g=0.40, 0.50
for key, r in results.items():
    if key.startswith('n12_') and 'cv' in r:
        n12_new[r['g']] = r['cv']

n16 = {}
for key, r in results.items():
    if key.startswith('n16_') and 'cv' in r:
        n16[r['g']] = r['cv']

print(f"n=8:  {dict(sorted(n8.items()))}")
print(f"n=12: {dict(sorted(n12_new.items()))}")
print(f"n=16: {dict(sorted(n16.items()))}")

# Crossing analysis
print("\n=== CROSSINGS ===")
for g in sorted(set(n8.keys()) & set(n12_new.keys())):
    d = n12_new[g] - n8[g]
    print(f"  g={g:.2f}: n8={n8[g]:.4f}, n12={n12_new[g]:.4f}, delta={d:+.4f}")

with open('results/sprint_045d_verify_n16.json', 'w') as f:
    json.dump({'sprint': '045d', 'q': q, 'chi': chi, 'results': results,
               'consolidated': {'n8': n8, 'n12': {str(k): v for k, v in n12_new.items()},
                                'n16': {str(k): v for k, v in n16.items()}}}, f, indent=2)
print("\nResults saved to results/sprint_045d_verify_n16.json")
