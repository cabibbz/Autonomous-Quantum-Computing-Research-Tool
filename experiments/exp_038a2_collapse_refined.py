"""
Sprint 038a2: Refined data collapse analysis.

1. Collapse quality for LARGE sizes only (n=16,24,32) — less finite-size contamination
2. Corrections-to-scaling: try CV(h,n) = F(x) + n^{-omega}*G(x) form
3. New DMRG data at n=50 for additional collapse test
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar, minimize
import json, time, warnings
warnings.filterwarnings('ignore')
import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)
from tenpy.models.tf_ising import TFIChain
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS
from scipy import linalg as la

t0 = time.time()

# Load consolidated data from 038a
with open('results/sprint_038a_data_collapse.json') as f:
    prior = json.load(f)

all_data = {}
for n_str, points in prior['sizes'].items():
    n = int(n_str)
    all_data[n] = {}
    for h_str, cv in points.items():
        all_data[n][float(h_str)] = cv

# ========== Generate n=50 data ==========
print("=== Sprint 038a2: Refined Collapse + n=50 Data ===\n")

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

def dmrg_mi_cv_tfim(n, h_J, chi_max=100):
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

# n=50: select h values that will help collapse test
h_vals_50 = [0.90, 0.95, 0.97, 1.0, 1.05]
print("Computing n=50 MI-CV...")
all_data[50] = {}
results_50 = {}

for h_J in h_vals_50:
    elapsed = time.time() - t0
    if elapsed > 45:
        print(f"  TIMEOUT at {elapsed:.0f}s")
        break
    t1 = time.time()
    cv = dmrg_mi_cv_tfim(50, h_J, chi_max=100)
    dt = time.time() - t1
    print(f"  h/J={h_J:.2f}: CV={cv:.4f} ({dt:.1f}s)")
    all_data[50][h_J] = cv
    results_50[str(h_J)] = {'cv': float(cv), 'time': float(dt)}

# ========== Collapse analysis with all sizes including n=50 ==========
print("\n=== Collapse Quality: All sizes vs Large sizes only ===")

def collapse_quality(nu, h_c, sizes):
    total_msd = 0.0
    n_pairs = 0
    for i in range(len(sizes)):
        for j in range(i+1, len(sizes)):
            n1, n2 = sizes[i], sizes[j]
            if n1 not in all_data or n2 not in all_data:
                continue
            h_arr1 = np.array(sorted(all_data[n1].keys()))
            cv_arr1 = np.array([all_data[n1][h] for h in sorted(all_data[n1].keys())])
            x_arr1 = (h_arr1 - h_c) * n1**(1.0/nu)

            h_arr2 = np.array(sorted(all_data[n2].keys()))
            cv_arr2 = np.array([all_data[n2][h] for h in sorted(all_data[n2].keys())])
            x_arr2 = (h_arr2 - h_c) * n2**(1.0/nu)

            x_min = max(x_arr1.min(), x_arr2.min())
            x_max = min(x_arr1.max(), x_arr2.max())
            if x_max <= x_min:
                continue

            f1 = interp1d(x_arr1, cv_arr1, kind='linear')
            f2 = interp1d(x_arr2, cv_arr2, kind='linear')
            x_common = np.linspace(x_min, x_max, 50)
            msd = np.mean((f1(x_common) - f2(x_common))**2)
            total_msd += msd
            n_pairs += 1
    return total_msd / n_pairs if n_pairs > 0 else 1e10

# Test different size subsets
size_sets = {
    'all (8-50)': [8, 12, 16, 24, 32, 50],
    'large (16-50)': [16, 24, 32, 50],
    'largest (24-50)': [24, 32, 50],
}

results = {
    'experiment': '038a2',
    'n50_data': results_50,
    'collapse_by_subset': {},
}

for label, sizes in size_sets.items():
    # Optimize nu with h_c=1.0 fixed
    opt = minimize_scalar(lambda nu: collapse_quality(nu, 1.0, sizes),
                          bounds=(0.3, 3.0), method='bounded')
    q_opt = opt.fun
    nu_opt = opt.x
    q_ising = collapse_quality(1.0, 1.0, sizes)

    # Joint optimize
    res = minimize(lambda p: collapse_quality(p[0], p[1], sizes),
                   [nu_opt, 1.0], method='Nelder-Mead',
                   options={'xatol': 0.001, 'fatol': 1e-8})

    print(f"\n{label}:")
    print(f"  nu_opt(hc=1) = {nu_opt:.3f}, quality = {q_opt:.6f}")
    print(f"  Ising (nu=1,hc=1): quality = {q_ising:.6f}, ratio = {q_ising/q_opt:.2f}x")
    print(f"  Joint opt: nu={res.x[0]:.3f}, hc={res.x[1]:.4f}, quality={res.fun:.6f}")

    results['collapse_by_subset'][label] = {
        'sizes': sizes,
        'nu_opt_hc1': float(nu_opt),
        'quality_opt_hc1': float(q_opt),
        'quality_ising': float(q_ising),
        'ratio_ising_to_opt': float(q_ising / q_opt) if q_opt > 0 else float('inf'),
        'joint_nu': float(res.x[0]),
        'joint_hc': float(res.x[1]),
        'joint_quality': float(res.fun),
    }

# ========== Pairwise collapse for ALL pairs including n=50 ==========
print("\n=== Pairwise Collapse Quality (nu=1.0, hc=1.0) ===")
sizes_all = [8, 12, 16, 24, 32, 50]
pair_q = {}
for i in range(len(sizes_all)):
    for j in range(i+1, len(sizes_all)):
        if sizes_all[i] not in all_data or sizes_all[j] not in all_data:
            continue
        q = collapse_quality(1.0, 1.0, [sizes_all[i], sizes_all[j]])
        label = f"{sizes_all[i]}_{sizes_all[j]}"
        pair_q[label] = float(q)
        print(f"  ({sizes_all[i]:2d},{sizes_all[j]:2d}): {q:.6f}")
results['pairwise_ising'] = pair_q

# ========== Slope at transition ==========
print("\n=== Transition slope at h=1.0 for each size ===")
# Compute dCV/dh at h=1.0 by finite difference from nearest points
for n in sorted(all_data.keys()):
    h_arr = np.array(sorted(all_data[n].keys()))
    cv_arr = np.array([all_data[n][h] for h in sorted(all_data[n].keys())])
    # Find points bracketing h=1.0
    below = h_arr[h_arr < 1.0]
    above = h_arr[h_arr > 1.0]
    if len(below) > 0 and len(above) > 0:
        h_lo = below[-1]
        h_hi = above[0]
        cv_lo = all_data[n][h_lo]
        cv_hi = all_data[n][h_hi]
        slope = (cv_hi - cv_lo) / (h_hi - h_lo)
        print(f"  n={n:3d}: slope = {slope:.2f} (h=[{h_lo:.2f},{h_hi:.2f}])")

results['total_runtime'] = time.time() - t0
with open('results/sprint_038a2_collapse_refined.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nTotal runtime: {time.time()-t0:.1f}s")
