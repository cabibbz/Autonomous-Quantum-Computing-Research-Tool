"""
Sprint 039a2: Potts n=16 continuation — g=1.0, 1.05, 1.1
"""

import numpy as np
from scipy import linalg as la
import json, time, warnings, sys
warnings.filterwarnings('ignore')
import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)
from tenpy.models.clock import ClockChain
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS

t0 = time.time()
d = 3

gm_mats = []
m = np.zeros((3,3), dtype=complex); m[0,1] = 1; m[1,0] = 1; gm_mats.append(('gm1', m.copy()))
m = np.zeros((3,3), dtype=complex); m[0,1] = -1j; m[1,0] = 1j; gm_mats.append(('gm2', m.copy()))
m = np.zeros((3,3), dtype=complex); m[0,0] = 1; m[1,1] = -1; gm_mats.append(('gm3', m.copy()))
m = np.zeros((3,3), dtype=complex); m[0,2] = 1; m[2,0] = 1; gm_mats.append(('gm4', m.copy()))
m = np.zeros((3,3), dtype=complex); m[0,2] = -1j; m[2,0] = 1j; gm_mats.append(('gm5', m.copy()))
m = np.zeros((3,3), dtype=complex); m[1,2] = 1; m[2,1] = 1; gm_mats.append(('gm6', m.copy()))
m = np.zeros((3,3), dtype=complex); m[1,2] = -1j; m[2,1] = 1j; gm_mats.append(('gm7', m.copy()))
m = np.zeros((3,3), dtype=complex); m[0,0] = 1; m[1,1] = 1; m[2,2] = -2; m /= np.sqrt(3); gm_mats.append(('gm8', m.copy()))

op_names_gm = [name for name, _ in gm_mats]
op_mats_all = [np.eye(3, dtype=complex)] + [mat for _, mat in gm_mats]

def entropy(rho):
    ev = la.eigvalsh(rho); ev = ev[ev > 1e-15]
    return float(-np.sum(ev * np.log2(ev)))

def mi_cv_potts(n, g_J, chi_max=30):
    model = ClockChain({'L': n, 'q': 3, 'J': 1.0, 'g': g_J,
                         'bc_MPS': 'finite', 'conserve': None})
    site = model.lat.site(0)
    for name, mat in gm_mats:
        if name not in site.opnames:
            site.add_op(name, mat)
    psi = MPS.from_product_state(model.lat.mps_sites(), [0]*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-8,
        'trunc_params': {'chi_max': chi_max}, 'max_sweeps': 8,
    })
    eng.run()
    exp_vals = {name: psi.expectation_value(name) for name in op_names_gm}
    corr = {}
    for a in op_names_gm:
        for b in op_names_gm:
            corr[(a, b)] = psi.correlation_function(a, b)
    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            rho_ij = np.eye(d*d, dtype=complex) / (d*d)
            for idx_a, a_name in enumerate(op_names_gm):
                a_mat = op_mats_all[idx_a + 1]
                ev_i = complex(exp_vals[a_name][i])
                ev_j = complex(exp_vals[a_name][j])
                rho_ij += (ev_i / (2*d)) * np.kron(a_mat, np.eye(d))
                rho_ij += (ev_j / (2*d)) * np.kron(np.eye(d), a_mat)
                for idx_b, b_name in enumerate(op_names_gm):
                    b_mat = op_mats_all[idx_b + 1]
                    c_val = complex(corr[(a_name, b_name)][i, j])
                    rho_ij += (c_val / 4.0) * np.kron(a_mat, b_mat)
            rho_ij = (rho_ij + rho_ij.conj().T) / 2
            ev = la.eigvalsh(rho_ij)
            if np.min(ev) < -1e-10:
                ev_pos = np.maximum(ev, 0)
                evec = la.eigh(rho_ij)[1]
                rho_ij = evec @ np.diag(ev_pos) @ evec.T
                rho_ij /= np.trace(rho_ij)
            rho_d = rho_ij.reshape(d, d, d, d)
            rho_i = np.trace(rho_d, axis1=1, axis2=3)
            rho_j = np.trace(rho_d, axis1=0, axis2=2)
            mi = entropy(rho_i) + entropy(rho_j) - entropy(rho_ij)
            mi_vals.append(max(float(np.real(mi)), 0))
    mi_pos = [m for m in mi_vals if m > 1e-10]
    if len(mi_pos) < 2: return 0.0
    return float(np.std(mi_pos) / np.mean(mi_pos))

print("=== Sprint 039a2: Potts n=16 continuation ===\n")
sys.stdout.flush()

# Load prior n=16 data
with open('results/sprint_039a_potts_n16.json') as f:
    prior = json.load(f)

results = prior.copy()
results['experiment'] = '039a2'

g_vals = [1.0, 1.05, 1.1]
for g_J in g_vals:
    elapsed = time.time() - t0
    if elapsed > 48:
        print(f"  TIMEOUT at {elapsed:.0f}s")
        break
    t1 = time.time()
    cv = mi_cv_potts(16, g_J, chi_max=30)
    dt = time.time() - t1
    print(f"  g/J={g_J:.2f}: CV={cv:.4f} ({dt:.1f}s)")
    sys.stdout.flush()
    results['data'][str(g_J)] = {'cv': float(cv), 'time': float(dt)}
    results['total_runtime'] = time.time() - t0
    with open('results/sprint_039a_potts_n16.json', 'w') as f:
        json.dump(results, f, indent=2)

print(f"\nTotal runtime: {time.time()-t0:.1f}s")
print(f"Total n=16 points: {len(results['data'])}")
