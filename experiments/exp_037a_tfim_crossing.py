"""
Sprint 037a: TFIM MI-CV crossing points for critical exponent extraction.

Dense h/J sweep near criticality for multiple sizes.
Find crossing points h_c(n1,n2), then fit h_c = h_c(inf) + a*n^{-1/nu}.
Ising universality predicts nu=1.

Uses DMRG for all sizes (faster and consistent).
"""

import numpy as np
from scipy import linalg as la
import json, time, warnings
warnings.filterwarnings('ignore')
import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)
from tenpy.models.tf_ising import TFIChain
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS

t0 = time.time()

I2 = np.eye(2)
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
sz = np.array([[1,0],[0,-1]])
pauli = {'Sigmax': sx, 'Sigmay': sy, 'Sigmaz': sz}
ops = ['Sigmax', 'Sigmay', 'Sigmaz']

def entropy(rho):
    ev = la.eigvalsh(rho); ev = ev[ev > 1e-15]
    return float(-np.sum(ev * np.log2(ev)))

def mi_from_rho2(rho_ij):
    rho_i = np.trace(rho_ij.reshape(2,2,2,2), axis1=1, axis2=3)
    rho_j = np.trace(rho_ij.reshape(2,2,2,2), axis1=0, axis2=2)
    return entropy(rho_i) + entropy(rho_j) - entropy(rho_ij)

def mi_cv(mi_values):
    mi_pos = [m for m in mi_values if m > 1e-10]
    if len(mi_pos) < 2: return 0.0
    return float(np.std(mi_pos) / np.mean(mi_pos))

def dmrg_mi_cv_tfim(n, h_J, chi_max=60):
    model = TFIChain({'L': n, 'J': 1.0, 'g': h_J, 'bc_MPS': 'finite'})
    psi = MPS.from_product_state(model.lat.mps_sites(), ['up']*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max}, 'max_sweeps': 15,
    })
    eng.run()

    exp_v = {op: psi.expectation_value(op) for op in ops}
    corr = {(a,b): psi.correlation_function(a, b) for a in ops for b in ops}

    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            rho_ij = np.eye(4, dtype=complex) / 4.0
            for a in ops:
                rho_ij += float(np.real(exp_v[a][i])) * np.kron(pauli[a], I2) / 4.0
                rho_ij += float(np.real(exp_v[a][j])) * np.kron(I2, pauli[a]) / 4.0
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

    return mi_cv(mi_vals)

print("=== Sprint 037a: TFIM MI-CV Crossing Points ===\n")

# Focused sweep: 12 points around criticality
h_vals = [0.80, 0.85, 0.90, 0.92, 0.94, 0.96, 0.98, 1.00, 1.02, 1.05, 1.10, 1.20]

sizes = [8, 12, 16, 24, 32]
chi_map = {8: 30, 12: 50, 16: 60, 24: 80, 32: 80}

results = {'experiment': '037a', 'sizes': {}, 'crossings': {}}

# Time calibration: one point for n=32
print("Calibrating n=32 timing...")
t_cal = time.time()
_ = dmrg_mi_cv_tfim(32, 1.0, chi_max=80)
dt_32 = time.time() - t_cal
print(f"  n=32 single point: {dt_32:.1f}s")
n_32_budget = int(40 / dt_32)  # how many n=32 points fit in 40s
print(f"  Budget: ~{n_32_budget} points for n=32\n")

for n in sizes:
    elapsed = time.time() - t0
    remaining = 55 - elapsed
    if remaining < 3:
        print(f"TIMEOUT at {elapsed:.0f}s, stopping")
        break

    # For larger sizes, use fewer h values
    if n <= 16:
        h_sweep = h_vals
    elif n == 24:
        h_sweep = [0.85, 0.90, 0.94, 0.98, 1.00, 1.05, 1.10]
    else:  # n=32
        h_sweep = [0.90, 0.95, 1.00, 1.05, 1.10]

    chi = chi_map[n]
    print(f"--- n={n} (chi={chi}, {len(h_sweep)} points) ---")
    results['sizes'][str(n)] = {}

    for h_J in h_sweep:
        elapsed = time.time() - t0
        if elapsed > 52:
            print(f"  TIMEOUT at {elapsed:.0f}s")
            break

        t1 = time.time()
        try:
            cv = dmrg_mi_cv_tfim(n, h_J, chi_max=chi)
            dt = time.time() - t1
            print(f"  h/J={h_J:.2f}: CV={cv:.4f} ({dt:.1f}s)")
            results['sizes'][str(n)][str(h_J)] = {'cv': float(cv), 'time': float(dt)}
        except Exception as e:
            print(f"  h/J={h_J:.2f}: FAILED: {e}")

    # Save after each size
    with open('results/sprint_037a_tfim_crossing.json', 'w') as f:
        json.dump(results, f, indent=2)

# Find crossing points
print("\n=== Crossing Point Analysis ===")
size_keys = sorted([int(k) for k in results['sizes'].keys()])

for i in range(len(size_keys) - 1):
    n1, n2 = size_keys[i], size_keys[i+1]
    d1 = results['sizes'][str(n1)]
    d2 = results['sizes'][str(n2)]

    common_h = sorted(set(d1.keys()) & set(d2.keys()), key=float)
    if len(common_h) < 3:
        print(f"  n=({n1},{n2}): insufficient common h values ({len(common_h)})")
        continue

    h_arr = np.array([float(h) for h in common_h])
    cv1 = np.array([d1[h]['cv'] for h in common_h])
    cv2 = np.array([d2[h]['cv'] for h in common_h])
    diff = cv1 - cv2

    crossings = []
    for j in range(len(diff) - 1):
        if diff[j] * diff[j+1] < 0:
            h_cross = h_arr[j] + (h_arr[j+1] - h_arr[j]) * abs(diff[j]) / (abs(diff[j]) + abs(diff[j+1]))
            crossings.append(float(h_cross))

    key = f"{n1}_{n2}"
    n_eff = np.sqrt(n1 * n2)  # geometric mean
    results['crossings'][key] = {'crossings': crossings, 'n_eff': float(n_eff)}
    print(f"  n=({n1},{n2}), n_eff={n_eff:.1f}: h_c = {crossings}")
    print(f"    diff = {dict(zip([f'{h:.2f}' for h in h_arr], [f'{d:.4f}' for d in diff]))}")

# Extract nu
crossing_pts = []
for key, data in results['crossings'].items():
    n1, n2 = map(int, key.split('_'))
    n_eff = data['n_eff']
    for h_c in data['crossings']:
        if 0.8 < h_c < 1.15:
            crossing_pts.append((n_eff, h_c))

if len(crossing_pts) >= 2:
    print("\n=== Critical Exponent Estimates ===")
    crossing_pts.sort()
    for ne, hc in crossing_pts:
        print(f"  n_eff={ne:.1f}: h_c = {hc:.4f}, shift = {abs(hc - 1.0):.4f}")

    # Two-point nu estimates using h_c(∞)=1.0
    for i in range(len(crossing_pts)-1):
        n1, h1 = crossing_pts[i]
        n2, h2 = crossing_pts[i+1]
        if abs(h1 - 1.0) > 1e-4 and abs(h2 - 1.0) > 1e-4:
            inv_nu = np.log(abs(h1 - 1.0) / abs(h2 - 1.0)) / np.log(n2 / n1)
            print(f"  Two-point ({n1:.0f},{n2:.0f}): 1/nu = {inv_nu:.3f}, nu = {1/inv_nu:.3f}")
            results['nu_estimates'] = results.get('nu_estimates', [])
            results['nu_estimates'].append({
                'n1': float(n1), 'n2': float(n2),
                'h1': float(h1), 'h2': float(h2),
                'inv_nu': float(inv_nu), 'nu': float(1/inv_nu)
            })

    if len(crossing_pts) >= 3:
        from scipy.optimize import curve_fit
        def model(n, a, inv_nu):
            return 1.0 + a * n**(-inv_nu)
        try:
            ns = np.array([x[0] for x in crossing_pts])
            hs = np.array([x[1] for x in crossing_pts])
            popt, pcov = curve_fit(model, ns, hs, p0=[-1.0, 1.0])
            perr = np.sqrt(np.diag(pcov))
            nu = 1.0/popt[1]
            print(f"\n  Fit (h_c_inf=1.0 fixed): 1/nu = {popt[1]:.3f} ± {perr[1]:.3f}")
            print(f"  => nu = {nu:.3f} (Ising exact: 1.0)")
            results['fit'] = {'inv_nu': float(popt[1]), 'inv_nu_err': float(perr[1]),
                              'nu': float(nu), 'a': float(popt[0])}
        except Exception as e:
            print(f"  Fit failed: {e}")

results['total_runtime'] = time.time() - t0
with open('results/sprint_037a_tfim_crossing.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nTotal runtime: {time.time()-t0:.1f}s")
