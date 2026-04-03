"""
Sprint 042b: True q=5 Potts n=8 sweep — low g region.
CV=1.49 at g=0.8, so transition is likely g<0.5. Sweep g=0.1,0.3,0.5.
~17s/point budget, 3 points in 60s.
"""

import numpy as np
from scipy import linalg as la
import json, time, warnings, sys
warnings.filterwarnings('ignore')
import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)

from tenpy.models.model import CouplingMPOModel, NearestNeighborModel
from tenpy.networks.site import Site
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS
import tenpy.linalg.np_conserved as npc

t0 = time.time()
q = 5

class PottsSite(Site):
    def __init__(self, q):
        leg = npc.LegCharge.from_trivial(q)
        names = [str(a) for a in range(q)]
        Site.__init__(self, leg, names, sort_charge=False)
        for a in range(q):
            P = np.zeros((q, q), dtype=complex); P[a, a] = 1.0
            self.add_op(f'P{a}', P)
        X = np.zeros((q, q), dtype=complex)
        for a in range(q): X[(a+1)%q, a] = 1.0
        self.add_op('X', X, hc='Xhc')
        self.add_op('Xhc', X.conj().T, hc='X')
        self.add_op('Xphc', X + X.conj().T, hc='Xphc')

class PottsModel(CouplingMPOModel):
    def init_sites(self, model_params):
        return PottsSite(model_params.get('q', 5))
    def init_terms(self, model_params):
        J = model_params.get('J', 1.0)
        g = model_params.get('g', 1.0)
        qq = model_params.get('q', 5)
        for a in range(qq):
            self.add_coupling(-J, 0, f'P{a}', 0, f'P{a}', 1)
        self.add_onsite(-g, 0, 'Xphc')

class PottsChain(PottsModel, NearestNeighborModel):
    pass

# Gell-Mann basis
gm_mats = []
for j in range(q):
    for k in range(j+1, q):
        m = np.zeros((q, q), dtype=complex); m[j,k]=1; m[k,j]=1
        gm_mats.append((f'sym_{j}{k}', m.copy()))
for j in range(q):
    for k in range(j+1, q):
        m = np.zeros((q, q), dtype=complex); m[j,k]=-1j; m[k,j]=1j
        gm_mats.append((f'asym_{j}{k}', m.copy()))
for l in range(1, q):
    m = np.zeros((q, q), dtype=complex)
    norm = np.sqrt(2.0/(l*(l+1)))
    for j in range(l): m[j,j]=norm
    m[l,l] = -l*norm
    gm_mats.append((f'diag_{l}', m.copy()))
op_names = [name for name, _ in gm_mats]
op_mats_all = [np.eye(q, dtype=complex)] + [mat for _, mat in gm_mats]

def entropy(rho):
    ev = la.eigvalsh(rho); ev = ev[ev>1e-15]
    return float(-np.sum(ev*np.log2(ev)))

def mi_cv_potts(n, g_J, chi_max=15):
    model = PottsChain({'L': n, 'q': q, 'J': 1.0, 'g': g_J,
                         'bc_MPS': 'finite', 'conserve': None})
    site = model.lat.site(0)
    for name, mat in gm_mats:
        if name not in site.opnames:
            site.add_op(name, mat)
    psi = MPS.from_product_state(model.lat.mps_sites(), [0]*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-6,
        'trunc_params': {'chi_max': chi_max}, 'max_sweeps': 6,
    })
    E0, _ = eng.run()
    exp_vals = {name: psi.expectation_value(name) for name in op_names}
    corr = {}
    for a in op_names:
        for b in op_names:
            corr[(a,b)] = psi.correlation_function(a, b)
    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            rho_ij = np.eye(q*q, dtype=complex)/(q*q)
            for idx_a, a_name in enumerate(op_names):
                a_mat = op_mats_all[idx_a+1]
                ev_i = complex(exp_vals[a_name][i])
                ev_j = complex(exp_vals[a_name][j])
                rho_ij += (ev_i/(2*q))*np.kron(a_mat, np.eye(q))
                rho_ij += (ev_j/(2*q))*np.kron(np.eye(q), a_mat)
                for idx_b, b_name in enumerate(op_names):
                    b_mat = op_mats_all[idx_b+1]
                    c_val = complex(corr[(a_name,b_name)][i,j])
                    rho_ij += (c_val/4.0)*np.kron(a_mat, b_mat)
            rho_ij = (rho_ij + rho_ij.conj().T)/2
            ev = la.eigvalsh(rho_ij)
            if np.min(ev) < -1e-10:
                ev_pos = np.maximum(ev, 0)
                evec = la.eigh(rho_ij)[1]
                rho_ij = evec @ np.diag(ev_pos) @ evec.T
                rho_ij /= np.trace(rho_ij)
            rho_d = rho_ij.reshape(q,q,q,q)
            rho_i = np.trace(rho_d, axis1=1, axis2=3)
            rho_j = np.trace(rho_d, axis1=0, axis2=2)
            mi = entropy(rho_i) + entropy(rho_j) - entropy(rho_ij)
            mi_vals.append(max(float(np.real(mi)), 0))
    mi_pos = [m for m in mi_vals if m > 1e-10]
    if len(mi_pos) < 2: return 0.0, E0
    return float(np.std(mi_pos)/np.mean(mi_pos)), E0

print("=== Sprint 042b: True q=5 Potts n=8 — Low g Sweep ===\n")
sys.stdout.flush()

results = {'experiment': '042b', 'model': 'q=5 Potts', 'n': 8, 'chi_max': 15, 'data': {}}

for g_J in [0.1, 0.3, 0.5]:
    elapsed = time.time() - t0
    if elapsed > 50:
        print(f"  TIMEOUT at {elapsed:.0f}s"); break
    t1 = time.time()
    try:
        cv, E0 = mi_cv_potts(8, g_J)
        dt = time.time()-t1
        print(f"  g/J={g_J:.2f}: CV={cv:.4f}, E0={E0:.6f} ({dt:.1f}s)")
        results['data'][str(g_J)] = {'cv': float(cv), 'E0': float(E0), 'time': float(dt)}
    except Exception as e:
        print(f"  g/J={g_J:.2f}: FAILED: {e}")
        results['data'][str(g_J)] = {'cv': None, 'error': str(e)}
    sys.stdout.flush()
    results['total_runtime'] = time.time()-t0
    with open('results/sprint_042b_potts_q5_n8_low.json', 'w') as f:
        json.dump(results, f, indent=2)

print(f"\nTotal: {time.time()-t0:.1f}s")
# Reference: g=0.8 gave CV=1.4922 (from 042a)
print("Reference: g=0.8 CV=1.4922 (042a)")
