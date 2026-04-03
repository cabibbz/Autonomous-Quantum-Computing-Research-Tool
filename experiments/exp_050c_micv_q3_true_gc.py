#!/usr/bin/env python3
"""Sprint 050c: MI-CV at the TRUE q=3 Potts critical point g_c=1/3.

Key question: do MI-CV crossings appear at the actual critical point?
Prior crossings (Sprint 038-039 at g≈0.9) were in the disordered phase.

Self-duality gives g_c = J/3 ≈ 0.333 exactly.
Scan g near 1/3 at n=8,12 to see if crossing curves appear.
Use direct MPS contraction for all-pairs MI (reliable for d=3).
"""
import numpy as np
import numpy.linalg as la
import json, time, warnings
warnings.filterwarnings('ignore')

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


class PottsChain(CouplingMPOModel, NearestNeighborModel):
    def init_sites(self, model_params):
        return PottsSite(model_params.get('q', 3))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 1.0)
        q = model_params.get('q', 3)
        for a in range(q):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'Xphc')


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


def mi_cv(n, q, g, chi_max=40):
    t0 = time.time()
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g, 'bc_MPS': 'finite'})
    np.random.seed(42 + n + int(g * 1000))
    init = [np.random.randint(q) for _ in range(n)]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max, 'svd_min': 1e-12},
        'max_sweeps': 40,
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
    mean_mi = float(np.mean(mi_pos)) if mi_pos else 0.0
    return cv, float(E0), mean_mi, mi_vals, t_dmrg, t_mi


# Timing test
q = 3
print("=== Timing test: q=3, n=8, chi=40, g=0.333 ===", flush=True)
cv, E0, mean_mi, mi_vals, t_dmrg, t_mi = mi_cv(8, q, 1/3, chi_max=40)
print(f"  CV={cv:.4f}, E0={E0:.6f}, mean_MI={mean_mi:.4f}, DMRG={t_dmrg:.1f}s, MI={t_mi:.1f}s", flush=True)

# g scan near g_c = 1/3 ≈ 0.333
# Ordered phase < 0.333, disordered > 0.333
g_values = [0.15, 0.20, 0.24, 0.28, 0.30, 0.32, 0.333, 0.35, 0.38, 0.42, 0.50]

results = {}
for n in [8, 12]:
    chi = 40
    print(f"\n=== q=3 Potts MI-CV: n={n}, chi={chi} ===", flush=True)
    results[n] = []
    t_start = time.time()
    for g in g_values:
        if time.time() - t_start > 230:
            print(f"  Time limit at n={n}", flush=True)
            break
        t0 = time.time()
        cv, E0, mean_mi, mi_vals, t_dmrg, t_mi = mi_cv(n, q, g, chi_max=chi)
        dt = time.time() - t0
        print(f"  g={g:.3f}: CV={cv:.4f}, mean_MI={mean_mi:.4f}, DMRG={t_dmrg:.1f}s, MI={t_mi:.1f}s", flush=True)
        results[n].append({
            'g': float(g), 'cv': cv, 'E0': E0, 'mean_mi': mean_mi,
            'mi_vals': mi_vals, 'time': dt
        })

    # Save incrementally
    with open('results/sprint_050c_micv_q3.json', 'w') as f:
        json.dump({'sprint': '050c', 'q': 3, 'chi_max': chi,
                   'g_c_selfdual': 1/3,
                   'data': {str(k): v for k, v in results.items()}}, f, indent=2)

# === Analysis ===
print("\n=== MI-CV Comparison: n=8 vs n=12 ===", flush=True)
print(f"{'g':>6s} | {'CV(8)':>7s} | {'CV(12)':>7s} | {'Δ=CV12-CV8':>10s} | Crossing?", flush=True)
print("-" * 60, flush=True)

if 8 in results and 12 in results:
    n8_cv = {round(d['g'], 3): d['cv'] for d in results[8]}
    n12_cv = {round(d['g'], 3): d['cv'] for d in results[12]}
    crossings = []
    prev_delta = None
    for g in sorted(set(n8_cv.keys()) & set(n12_cv.keys())):
        delta = n12_cv[g] - n8_cv[g]
        cross = ""
        if prev_delta is not None and delta * prev_delta < 0:
            cross = " ← CROSSING"
            crossings.append(g)
        print(f"  {g:.3f} | {n8_cv[g]:>7.4f} | {n12_cv[g]:>7.4f} | {delta:>+10.4f} | {cross}", flush=True)
        prev_delta = delta

    if crossings:
        print(f"\n*** MI-CV CROSSINGS found near g = {crossings} ***", flush=True)
    else:
        print(f"\n*** NO MI-CV crossings found ***", flush=True)

print("\nDone!", flush=True)
