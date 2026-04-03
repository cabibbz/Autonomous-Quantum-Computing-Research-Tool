"""
Sprint 037b2: XXZ n=32 near Δ=-1 first-order transition.
From 037b: FM phase has CV=0, jump at Δ≈-1. Need n=32 near transition.
"""

import numpy as np
from scipy import linalg as la
import json, time, warnings
warnings.filterwarnings('ignore')
import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)
from tenpy.models.xxz_chain import XXZChain
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS

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

def mi_cv(mi_values):
    mi_pos = [m for m in mi_values if m > 1e-10]
    if len(mi_pos) < 2: return 0.0
    return float(np.std(mi_pos) / np.mean(mi_pos))

def dmrg_mi_cv_xxz(n, delta, chi_max=80):
    model = XXZChain({'L': n, 'Jxx': 1.0, 'Jz': delta, 'hz': 0.0, 'bc_MPS': 'finite'})
    if delta < -1.0:
        init = ['up'] * n
    else:
        init = ['up', 'down'] * (n // 2)
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max}, 'max_sweeps': 20,
    })
    eng.run()

    spin_ops = ['Sp', 'Sm', 'Sz']
    exp_sp = psi.expectation_value('Sp')
    exp_sm = psi.expectation_value('Sm')
    exp_sz = psi.expectation_value('Sz')
    corr = {(a,b): psi.correlation_function(a, b) for a in spin_ops for b in spin_ops}

    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            sx_i = float(np.real(exp_sp[i] + exp_sm[i]))
            sy_i = float(np.real(-1j * (exp_sp[i] - exp_sm[i])))
            sz_i = float(np.real(2 * exp_sz[i]))
            sx_j = float(np.real(exp_sp[j] + exp_sm[j]))
            sy_j = float(np.real(-1j * (exp_sp[j] - exp_sm[j])))
            sz_j = float(np.real(2 * exp_sz[j]))

            SpSp = corr[('Sp','Sp')][i,j]; SpSm = corr[('Sp','Sm')][i,j]
            SmSp = corr[('Sm','Sp')][i,j]; SmSm = corr[('Sm','Sm')][i,j]
            SzSz = corr[('Sz','Sz')][i,j]
            SpSz = corr[('Sp','Sz')][i,j]; SmSz = corr[('Sm','Sz')][i,j]
            SzSp = corr[('Sz','Sp')][i,j]; SzSm = corr[('Sz','Sm')][i,j]

            sxsx = float(np.real(SpSp + SpSm + SmSp + SmSm))
            sysy = float(np.real(-(SpSp - SpSm - SmSp + SmSm)))
            szsz = float(np.real(4 * SzSz))
            sxsy = float(np.real(-1j * (SpSp - SpSm + SmSp - SmSm)))
            sysx = float(np.real(-1j * (SpSp + SpSm - SmSp - SmSm)))
            sxsz = float(np.real(2 * (SpSz + SmSz)))
            szsx = float(np.real(2 * (SzSp + SzSm)))
            sysz = float(np.real(-2j * (SpSz - SmSz)))
            szsy = float(np.real(-2j * (SzSp - SzSm)))

            rho_ij = np.eye(4, dtype=complex) / 4.0
            rho_ij += sx_i * np.kron(sx, I2) / 4.0 + sy_i * np.kron(sy, I2) / 4.0 + sz_i * np.kron(sz, I2) / 4.0
            rho_ij += sx_j * np.kron(I2, sx) / 4.0 + sy_j * np.kron(I2, sy) / 4.0 + sz_j * np.kron(I2, sz) / 4.0
            rho_ij += sxsx * np.kron(sx, sx) / 4.0 + sxsy * np.kron(sx, sy) / 4.0 + sxsz * np.kron(sx, sz) / 4.0
            rho_ij += sysx * np.kron(sy, sx) / 4.0 + sysy * np.kron(sy, sy) / 4.0 + sysz * np.kron(sy, sz) / 4.0
            rho_ij += szsx * np.kron(sz, sx) / 4.0 + szsy * np.kron(sz, sy) / 4.0 + szsz * np.kron(sz, sz) / 4.0

            rho_ij = np.real(rho_ij)
            ev = la.eigvalsh(rho_ij)
            if np.min(ev) < -1e-10:
                ev_pos = np.maximum(ev, 0)
                evec = la.eigh(rho_ij)[1]
                rho_ij = evec @ np.diag(ev_pos) @ evec.T
                rho_ij /= np.trace(rho_ij)
            mi_vals.append(max(mi_from_rho2(rho_ij), 0))

    return mi_cv(mi_vals)

print("=== Sprint 037b2: XXZ n=32 Near Δ=-1 ===\n")

# Only points near transition + a couple in XY phase
delta_vals = [-1.10, -1.05, -1.02, -1.00, -0.98, -0.95, -0.90, -0.50]

results = {'experiment': '037b2', 'n32': {}}

for delta in delta_vals:
    elapsed = time.time() - t0
    if elapsed > 50:
        print(f"TIMEOUT")
        break
    t1 = time.time()
    try:
        cv = dmrg_mi_cv_xxz(32, delta, chi_max=80)
        dt = time.time() - t1
        print(f"Δ={delta:6.2f}: CV={cv:.4f} ({dt:.1f}s)")
        results['n32'][str(delta)] = {'cv': float(cv), 'time': float(dt)}
    except Exception as e:
        print(f"Δ={delta:6.2f}: FAILED: {e}")

    with open('results/sprint_037b2_xxz_n32.json', 'w') as f:
        json.dump(results, f, indent=2)

# Combine with previous 037b data for full analysis
print("\n=== Combined Analysis (all sizes) ===")

try:
    with open('results/sprint_037b_xxz_firstorder.json') as f:
        prev = json.load(f)
except:
    prev = {'sizes': {}}

all_data = {}
for n_str, data in prev.get('sizes', {}).items():
    n = int(n_str)
    all_data[n] = {}
    for d_str, v in data.items():
        if isinstance(v, dict) and 'cv' in v:
            all_data[n][float(d_str)] = v['cv']

all_data[32] = {}
for d_str, v in results['n32'].items():
    all_data[32][float(d_str)] = v['cv']

# Gradient analysis
print("\n  Gradient at transition:")
for n in sorted(all_data.keys()):
    data = all_data[n]
    deltas = sorted(data.keys())
    cvs = [data[d] for d in deltas]

    max_grad = 0
    max_loc = None
    for i in range(len(deltas)-1):
        if abs(deltas[i+1] - deltas[i]) > 1e-6:
            grad = abs(cvs[i+1] - cvs[i]) / abs(deltas[i+1] - deltas[i])
            if grad > max_grad:
                max_grad = grad
                max_loc = (deltas[i] + deltas[i+1]) / 2

    if max_loc is not None:
        print(f"  n={n:3d}: max|dCV/dΔ| = {max_grad:.3f} at Δ ≈ {max_loc:.3f}")

    # Also check: CV at Δ=-0.98 vs Δ=-1.02
    cv_above = data.get(-0.98, None)
    cv_below = data.get(-1.02, None)
    if cv_above is not None and cv_below is not None:
        jump = cv_above - cv_below
        print(f"         CV jump (-1.02→-0.98) = {jump:.4f}")

# Compare with TFIM slope at h_c=1.0
print("\n  For comparison, TFIM slopes at h_c=1.0 (Sprint 037a):")
print("  n=8: ~1.7, n=16: ~3.9")
print("  XXZ first-order should scale FASTER than Ising (volume vs correlation length)")

results['total_runtime'] = time.time() - t0
with open('results/sprint_037b2_xxz_n32.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nTotal runtime: {time.time()-t0:.1f}s")
