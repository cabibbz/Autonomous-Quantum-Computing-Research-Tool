"""
Sprint 036a2: TFIM MI-CV at n=32 and n=50

Focused sweep to complete the size scaling from 036a.
Uses correlation function approach for all-pairs MI.
"""

import numpy as np
from scipy import linalg as la
import json, time, warnings
warnings.filterwarnings('ignore')

t0 = time.time()

I2 = np.eye(2)
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
sz = np.array([[1,0],[0,-1]])

def entropy(rho):
    ev = la.eigvalsh(rho); ev = ev[ev > 1e-15]
    return float(-np.sum(ev * np.log2(ev)))

def mi_from_rho2(rho_ij):
    rho_i = np.trace(rho_ij.reshape(2,2,2,2), axis1=1, axis2=3)
    rho_j = np.trace(rho_ij.reshape(2,2,2,2), axis1=0, axis2=2)
    return entropy(rho_i) + entropy(rho_j) - entropy(rho_ij)

def dmrg_all_mi_tfim(n, h_J, chi_max=100):
    import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)
    from tenpy.models.tf_ising import TFIChain
    from tenpy.algorithms import dmrg
    from tenpy.networks.mps import MPS

    model = TFIChain({'L': n, 'J': 1.0, 'g': h_J, 'bc_MPS': 'finite'})
    psi = MPS.from_product_state(model.lat.mps_sites(), ['up']*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max}, 'max_sweeps': 20,
    })
    E0, _ = eng.run()

    exp_x = psi.expectation_value('Sigmax')
    exp_y = psi.expectation_value('Sigmay')
    exp_z = psi.expectation_value('Sigmaz')

    ops = ['Sigmax', 'Sigmay', 'Sigmaz']
    corr = {}
    for a in ops:
        for b in ops:
            corr[(a,b)] = psi.correlation_function(a, b)

    pauli = {'Sigmax': sx, 'Sigmay': sy, 'Sigmaz': sz}
    exp_map = {'Sigmax': exp_x, 'Sigmay': exp_y, 'Sigmaz': exp_z}

    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            rho_ij = np.eye(4, dtype=complex) / 4.0
            for a in ops:
                rho_ij += float(np.real(exp_map[a][i])) * np.kron(pauli[a], I2) / 4.0
                rho_ij += float(np.real(exp_map[a][j])) * np.kron(I2, pauli[a]) / 4.0
            for a in ops:
                for b in ops:
                    rho_ij += float(np.real(corr[(a,b)][i,j])) * np.kron(pauli[a], pauli[b]) / 4.0
            rho_ij = np.real(rho_ij)
            ev = la.eigvalsh(rho_ij)
            if np.min(ev) < -1e-10:
                ev_pos = np.maximum(ev, 0)
                evec = la.eigh(rho_ij)[1]
                rho_ij = evec @ np.diag(ev_pos) @ evec.T
                rho_ij /= np.trace(rho_ij)
            mi_vals.append(max(mi_from_rho2(rho_ij), 0))

    return mi_vals, E0

def mi_cv(mi_values):
    mi_pos = [m for m in mi_values if m > 1e-10]
    if len(mi_pos) < 2: return 0.0
    return float(np.std(mi_pos) / np.mean(mi_pos))

print("=== Sprint 036a2: TFIM n=32,50 MI-CV ===\n")

results = {'experiment': '036a2', 'data': {}}

# n=32 sweep: key points around criticality
h_vals_32 = [0.5, 0.8, 0.9, 0.95, 1.0, 1.1, 1.3, 2.0]

print("--- n=32 (chi=100) ---")
results['data']['32'] = {}
for h_J in h_vals_32:
    elapsed = time.time() - t0
    if elapsed > 48:
        print(f"  TIMEOUT at {elapsed:.0f}s")
        break
    t1 = time.time()
    try:
        mi, _ = dmrg_all_mi_tfim(32, h_J, chi_max=100)
        cv = mi_cv(mi)
        mean_mi = float(np.mean([m for m in mi if m > 1e-10])) if any(m > 1e-10 for m in mi) else 0
        dt = time.time() - t1
        print(f"  h/J={h_J:.2f}: CV={cv:.4f}, mean={mean_mi:.6f} ({dt:.1f}s)")
        results['data']['32'][str(h_J)] = {
            'h_J': h_J, 'mi_cv': float(cv), 'mean_mi': mean_mi,
            'n_pairs': len(mi), 'time': float(dt)
        }
    except Exception as e:
        print(f"  h/J={h_J:.2f}: FAILED: {e}")

    with open('results/sprint_036a2_tfim_n32_50.json', 'w') as f:
        json.dump(results, f, indent=2)

# n=50 if time permits
elapsed = time.time() - t0
if elapsed < 35:
    print(f"\n--- n=50 (chi=64) ---")
    results['data']['50'] = {}
    h_vals_50 = [0.5, 0.9, 1.0, 1.1, 2.0]
    for h_J in h_vals_50:
        elapsed = time.time() - t0
        if elapsed > 48:
            break
        t1 = time.time()
        try:
            mi, _ = dmrg_all_mi_tfim(50, h_J, chi_max=64)
            cv = mi_cv(mi)
            mean_mi = float(np.mean([m for m in mi if m > 1e-10])) if any(m > 1e-10 for m in mi) else 0
            dt = time.time() - t1
            print(f"  h/J={h_J:.2f}: CV={cv:.4f}, mean={mean_mi:.6f} ({dt:.1f}s)")
            results['data']['50'][str(h_J)] = {
                'h_J': h_J, 'mi_cv': float(cv), 'mean_mi': mean_mi,
                'n_pairs': len(mi), 'time': float(dt)
            }
        except Exception as e:
            print(f"  h/J={h_J:.2f}: FAILED: {e}")

        with open('results/sprint_036a2_tfim_n32_50.json', 'w') as f:
            json.dump(results, f, indent=2)

# Summary
print(f"\n{'='*50}")
for n_str in sorted(results['data'].keys(), key=int):
    data = results['data'][n_str]
    if not data: continue
    print(f"\nn={n_str}:")
    for h_str in sorted(data.keys(), key=float):
        d = data[h_str]
        print(f"  h/J={float(h_str):5.2f}: CV={d['mi_cv']:.4f}")

results['total_runtime'] = time.time() - t0
with open('results/sprint_036a2_tfim_n32_50.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nRuntime: {time.time()-t0:.1f}s")
