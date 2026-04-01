"""
Sprint 037b: XXZ MI-CV at first-order FM transition (Δ=-1).

XXZ has first-order transition at Δ=-1 (FM→XY).
MI-CV should show discontinuous jump, not smooth inflection (Ising) or dome (BKT).
Test at n=8,16,32 to see how the jump sharpens with size.

XXZChain H = Jxx*(Sx Sx + Sy Sy) + Jz * Sz Sz
Sp/Sm conversion: sigma_x = Sp+Sm, sigma_y = -i(Sp-Sm), sigma_z = 2*Sz
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
    # FM initial state for Δ < -1, Néel for Δ > 0, half-filling otherwise
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
    corr = {}
    for a in spin_ops:
        for b in spin_ops:
            corr[(a,b)] = psi.correlation_function(a, b)

    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            sx_i = float(np.real(exp_sp[i] + exp_sm[i]))
            sy_i = float(np.real(-1j * (exp_sp[i] - exp_sm[i])))
            sz_i = float(np.real(2 * exp_sz[i]))
            sx_j = float(np.real(exp_sp[j] + exp_sm[j]))
            sy_j = float(np.real(-1j * (exp_sp[j] - exp_sm[j])))
            sz_j = float(np.real(2 * exp_sz[j]))

            SpSp = corr[('Sp','Sp')][i,j]
            SpSm = corr[('Sp','Sm')][i,j]
            SmSp = corr[('Sm','Sp')][i,j]
            SmSm = corr[('Sm','Sm')][i,j]
            SzSz = corr[('Sz','Sz')][i,j]
            SpSz = corr[('Sp','Sz')][i,j]
            SmSz = corr[('Sm','Sz')][i,j]
            SzSp = corr[('Sz','Sp')][i,j]
            SzSm = corr[('Sz','Sm')][i,j]

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
            rho_ij += sx_i * np.kron(sx, I2) / 4.0
            rho_ij += sy_i * np.kron(sy, I2) / 4.0
            rho_ij += sz_i * np.kron(sz, I2) / 4.0
            rho_ij += sx_j * np.kron(I2, sx) / 4.0
            rho_ij += sy_j * np.kron(I2, sy) / 4.0
            rho_ij += sz_j * np.kron(I2, sz) / 4.0
            rho_ij += sxsx * np.kron(sx, sx) / 4.0
            rho_ij += sxsy * np.kron(sx, sy) / 4.0
            rho_ij += sxsz * np.kron(sx, sz) / 4.0
            rho_ij += sysx * np.kron(sy, sx) / 4.0
            rho_ij += sysy * np.kron(sy, sy) / 4.0
            rho_ij += sysz * np.kron(sy, sz) / 4.0
            rho_ij += szsx * np.kron(sz, sx) / 4.0
            rho_ij += szsy * np.kron(sz, sy) / 4.0
            rho_ij += szsz * np.kron(sz, sz) / 4.0

            rho_ij = np.real(rho_ij)
            ev = la.eigvalsh(rho_ij)
            if np.min(ev) < -1e-10:
                ev_pos = np.maximum(ev, 0)
                evec = la.eigh(rho_ij)[1]
                rho_ij = evec @ np.diag(ev_pos) @ evec.T
                rho_ij /= np.trace(rho_ij)
            mi_vals.append(max(mi_from_rho2(rho_ij), 0))

    return mi_cv(mi_vals)

print("=== Sprint 037b: XXZ First-Order Transition (Δ=-1) ===\n")

# Sweep Δ from -2 to 0, dense near Δ=-1
delta_vals = [-2.0, -1.5, -1.3, -1.2, -1.15, -1.10, -1.05, -1.02,
              -1.00, -0.98, -0.95, -0.90, -0.80, -0.50, 0.0]

results = {'experiment': '037b', 'sizes': {}}

for n in [8, 16, 32]:
    elapsed = time.time() - t0
    if elapsed > 50:
        print(f"TIMEOUT")
        break

    chi = {8: 30, 16: 60, 32: 80}[n]
    print(f"--- n={n} (chi={chi}) ---")
    results['sizes'][str(n)] = {}

    for delta in delta_vals:
        elapsed = time.time() - t0
        if elapsed > 50:
            print(f"  TIMEOUT")
            break

        t1 = time.time()
        try:
            cv = dmrg_mi_cv_xxz(n, delta, chi_max=chi)
            dt = time.time() - t1
            print(f"  Δ={delta:6.2f}: CV={cv:.4f} ({dt:.1f}s)")
            results['sizes'][str(n)][str(delta)] = {'cv': float(cv), 'time': float(dt)}
        except Exception as e:
            print(f"  Δ={delta:6.2f}: FAILED: {e}")

    with open('results/sprint_037b_xxz_firstorder.json', 'w') as f:
        json.dump(results, f, indent=2)

# Analysis
print("\n=== Analysis ===")
for n_str in sorted(results['sizes'].keys(), key=int):
    data = results['sizes'][n_str]
    if not data:
        continue
    deltas = sorted(data.keys(), key=float)
    cvs = [data[d]['cv'] for d in deltas]

    # Find max CV gradient (steepest change)
    max_grad = 0
    max_grad_loc = None
    for i in range(len(deltas)-1):
        d1, d2 = float(deltas[i]), float(deltas[i+1])
        grad = abs(cvs[i+1] - cvs[i]) / abs(d2 - d1)
        if grad > max_grad:
            max_grad = grad
            max_grad_loc = (d1 + d2) / 2

    print(f"\nn={n_str}:")
    print(f"  Max |dCV/dΔ| = {max_grad:.3f} at Δ ≈ {max_grad_loc:.2f}")
    print(f"  CV range: [{min(cvs):.4f}, {max(cvs):.4f}]")

    # Check for discontinuity: is the gradient growing with n?
    results['sizes'][n_str]['max_gradient'] = {'value': max_grad, 'location': max_grad_loc}

# Compare gradient scaling
grads = []
for n_str in sorted(results['sizes'].keys(), key=int):
    data = results['sizes'][n_str]
    if 'max_gradient' in data:
        grads.append((int(n_str), data['max_gradient']['value'], data['max_gradient']['location']))

if len(grads) >= 2:
    print(f"\n=== Gradient Scaling (first-order test) ===")
    for n, g, loc in grads:
        print(f"  n={n}: max|dCV/dΔ| = {g:.3f} at Δ≈{loc:.2f}")
    # For first-order: gradient should grow as ~n (volume scaling)
    # For second-order: gradient grows as ~n^(1/ν)
    if len(grads) >= 2:
        for i in range(len(grads)-1):
            n1, g1, _ = grads[i]
            n2, g2, _ = grads[i+1]
            if g1 > 0:
                alpha = np.log(g2/g1) / np.log(n2/n1)
                print(f"  Scaling exponent ({n1},{n2}): α = {alpha:.3f}")

results['total_runtime'] = time.time() - t0
with open('results/sprint_037b_xxz_firstorder.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nTotal runtime: {time.time()-t0:.1f}s")
