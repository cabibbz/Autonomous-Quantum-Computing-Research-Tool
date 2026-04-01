"""
Sprint 037a3: Fill gaps for crossing analysis.
Need overlap between (n=16,n=24) and (n=24,n=32).

Missing: n=16 at h=0.93,0.95,0.97,1.00,1.03
         n=24 at 0.92,0.94,0.96
         n=32 at 0.98,1.00
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

print("=== Sprint 037a3: Fill-in Points ===\n")

# Priority ordered: most critical gaps first
tasks = [
    (16, 0.93, 60), (16, 0.95, 60), (16, 0.97, 60), (16, 1.00, 60),
    (24, 0.92, 80), (24, 0.94, 80), (24, 0.96, 80),
    (32, 0.98, 80), (32, 1.00, 80),
]

results = {'experiment': '037a3', 'new_points': {}}

for n, h_J, chi in tasks:
    elapsed = time.time() - t0
    if elapsed > 52:
        print(f"TIMEOUT at {elapsed:.0f}s")
        break

    t1 = time.time()
    try:
        cv = dmrg_mi_cv_tfim(n, h_J, chi_max=chi)
        dt = time.time() - t1
        print(f"n={n}, h/J={h_J:.2f}: CV={cv:.4f} ({dt:.1f}s)")
        key = f"{n}_{h_J}"
        results['new_points'][key] = {'n': n, 'h_J': h_J, 'cv': float(cv), 'time': float(dt)}
    except Exception as e:
        print(f"n={n}, h/J={h_J:.2f}: FAILED: {e}")

    with open('results/sprint_037a3_tfim_fill.json', 'w') as f:
        json.dump(results, f, indent=2)

# Now do crossing analysis with ALL data combined
print("\n=== Combined Analysis ===")

# Load all previous data
all_sizes = {}

with open('results/sprint_037a_tfim_crossing.json') as f:
    d = json.load(f)
    for n_str, data in d['sizes'].items():
        n = int(n_str)
        all_sizes[n] = {float(k): v['cv'] for k, v in data.items()}

with open('results/sprint_037a2_tfim_large.json') as f:
    d = json.load(f)
    for n_str, data in d['sizes'].items():
        n = int(n_str)
        if n not in all_sizes:
            all_sizes[n] = {}
        for h_str, v in data.items():
            all_sizes[n][float(h_str)] = v['cv']

# Add new fill-in points
for key, data in results['new_points'].items():
    n = data['n']
    h = data['h_J']
    if n not in all_sizes:
        all_sizes[n] = {}
    all_sizes[n][h] = data['cv']

# Print all data
for n in sorted(all_sizes.keys()):
    d = all_sizes[n]
    print(f"\nn={n}:")
    for h in sorted(d.keys()):
        print(f"  h/J={h:.2f}: CV={d[h]:.4f}")

# Find crossings
print("\n=== Crossings ===")
size_list = sorted(all_sizes.keys())
crossing_pts = []

for i in range(len(size_list) - 1):
    n1, n2 = size_list[i], size_list[i+1]
    d1, d2 = all_sizes[n1], all_sizes[n2]
    common_h = sorted(set(d1.keys()) & set(d2.keys()))

    h_arr = np.array(common_h)
    diff = np.array([d1[h] - d2[h] for h in common_h])

    crosses = []
    for j in range(len(diff) - 1):
        if diff[j] * diff[j+1] < 0:
            h_cross = h_arr[j] + (h_arr[j+1] - h_arr[j]) * abs(diff[j]) / (abs(diff[j]) + abs(diff[j+1]))
            crosses.append(float(h_cross))

    n_eff = np.sqrt(n1 * n2)
    print(f"  n=({n1},{n2}), n_eff={n_eff:.1f}: h_c = {[f'{x:.4f}' for x in crosses]}")
    for h_c in crosses:
        if 0.8 < h_c < 1.1:
            crossing_pts.append((n_eff, h_c))

results['all_crossings'] = crossing_pts

if len(crossing_pts) >= 2:
    crossing_pts.sort()
    print("\n=== Critical Exponent ===")
    for ne, hc in crossing_pts:
        print(f"  n_eff={ne:.1f}: h_c={hc:.4f}, |shift|={abs(hc-1):.4f}")

    for i in range(len(crossing_pts)-1):
        n1, h1 = crossing_pts[i]
        n2, h2 = crossing_pts[i+1]
        if abs(h1-1) > 1e-4 and abs(h2-1) > 1e-4:
            inv_nu = np.log(abs(h1-1)/abs(h2-1)) / np.log(n2/n1)
            print(f"  Two-point ({n1:.0f},{n2:.0f}): 1/nu={inv_nu:.3f}, nu={1/inv_nu:.3f}")

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
            print(f"\n  Fit: a={popt[0]:.3f}, 1/nu={popt[1]:.3f}±{perr[1]:.3f}")
            print(f"  => nu = {nu:.3f} (Ising: 1.0)")
            results['fit'] = {'a': float(popt[0]), 'inv_nu': float(popt[1]),
                              'inv_nu_err': float(perr[1]), 'nu': float(nu)}
        except Exception as e:
            print(f"  Fit failed: {e}")

results['total_runtime'] = time.time() - t0
with open('results/sprint_037a3_tfim_fill.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nTotal runtime: {time.time()-t0:.1f}s")
