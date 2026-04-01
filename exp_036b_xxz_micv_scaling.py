"""
Sprint 036b: XXZ MI-CV Finite-Size Scaling (n=8 to 32)

XXZChain uses Sp, Sm, Sz operators. Convert to Pauli correlations:
sigma_x = Sp + Sm, sigma_y = -i(Sp - Sm), sigma_z = 2*Sz
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

def kron_list(ops):
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r

def entropy(rho):
    ev = la.eigvalsh(rho); ev = ev[ev > 1e-15]
    return float(-np.sum(ev * np.log2(ev)))

def mi_from_rho2(rho_ij):
    rho_i = np.trace(rho_ij.reshape(2,2,2,2), axis1=1, axis2=3)
    rho_j = np.trace(rho_ij.reshape(2,2,2,2), axis1=0, axis2=2)
    return entropy(rho_i) + entropy(rho_j) - entropy(rho_ij)

def xxz_ground_exact(n, delta):
    Sp = np.array([[0,1],[0,0]]); Sm = np.array([[0,0],[1,0]])
    dim = 2**n; H = np.zeros((dim, dim), dtype=complex)
    for i in range(n-1):
        ops_pm = [I2]*n; ops_pm[i] = Sp; ops_pm[i+1] = Sm
        ops_mp = [I2]*n; ops_mp[i] = Sm; ops_mp[i+1] = Sp
        ops_zz = [I2]*n; ops_zz[i] = sz/2; ops_zz[i+1] = sz/2
        H += 0.5*(kron_list(ops_pm) + kron_list(ops_mp)) + delta * kron_list(ops_zz)
    H = np.real(H) if np.allclose(H.imag, 0) else H
    evals, evecs = la.eigh(H)
    return evecs[:, 0]

def all_pairs_mi_exact(psi, n):
    rho = np.outer(psi, psi.conj()).reshape([2]*(2*n))
    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            trace_q = [q for q in range(n) if q not in (i, j)]
            rt = rho.copy()
            for q in sorted(trace_q, reverse=True):
                nr = rt.ndim // 2
                rt = np.trace(rt, axis1=q, axis2=q+nr)
            mi_vals.append(max(mi_from_rho2(rt.reshape(4,4)), 0))
    return mi_vals

def dmrg_all_mi_xxz(n, delta, chi_max=100):
    """All-pairs MI using Sp/Sm/Sz correlations from XXZChain."""
    import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)
    from tenpy.models.xxz_chain import XXZChain
    from tenpy.algorithms import dmrg
    from tenpy.networks.mps import MPS

    model = XXZChain({'L': n, 'Jxx': 1.0, 'Jz': delta, 'hz': 0.0, 'bc_MPS': 'finite'})
    psi = MPS.from_product_state(model.lat.mps_sites(), ['up', 'down']*(n//2), bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max}, 'max_sweeps': 20,
    })
    E0, _ = eng.run()

    # Single-site: <Sz_i>, <Sp_i>, <Sm_i>
    exp_sz = psi.expectation_value('Sz')  # real
    exp_sp = psi.expectation_value('Sp')  # complex
    exp_sm = psi.expectation_value('Sm')  # complex

    # Two-point correlations with Sp, Sm, Sz
    # sigma_x = Sp + Sm, sigma_y = -i(Sp - Sm), sigma_z = 2*Sz
    spin_ops = ['Sp', 'Sm', 'Sz']
    corr = {}
    for a in spin_ops:
        for b in spin_ops:
            corr[(a,b)] = psi.correlation_function(a, b)

    # Reconstruct Pauli correlations:
    # <sx_i sx_j> = <(Sp+Sm)_i (Sp+Sm)_j> = <SpSp> + <SpSm> + <SmSp> + <SmSm>
    # <sy_i sy_j> = <(-i(Sp-Sm))_i (-i(Sp-Sm))_j> = -(SpSp - SpSm - SmSp + SmSm)
    # <sz_i sz_j> = <2Sz_i 2Sz_j> = 4*<SzSz>
    # <sx_i sy_j> = <(Sp+Sm)_i (-i)(Sp-Sm)_j> = -i(SpSp - SpSm + SmSp - SmSm)
    # etc.

    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            # Single-site Pauli expectations
            # <sigma_x> = <Sp> + <Sm>, <sigma_y> = -i(<Sp> - <Sm>), <sigma_z> = 2<Sz>
            sx_i = float(np.real(exp_sp[i] + exp_sm[i]))
            sy_i = float(np.real(-1j * (exp_sp[i] - exp_sm[i])))
            sz_i = float(np.real(2 * exp_sz[i]))
            sx_j = float(np.real(exp_sp[j] + exp_sm[j]))
            sy_j = float(np.real(-1j * (exp_sp[j] - exp_sm[j])))
            sz_j = float(np.real(2 * exp_sz[j]))

            # Two-point Pauli correlations
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

            # Build rho_ij = (1/4)(I4 + single-site + two-point)
            rho_ij = np.eye(4, dtype=complex) / 4.0
            # Single-site terms
            rho_ij += sx_i * np.kron(sx, I2) / 4.0
            rho_ij += sy_i * np.kron(sy, I2) / 4.0
            rho_ij += sz_i * np.kron(sz, I2) / 4.0
            rho_ij += sx_j * np.kron(I2, sx) / 4.0
            rho_ij += sy_j * np.kron(I2, sy) / 4.0
            rho_ij += sz_j * np.kron(I2, sz) / 4.0
            # Two-point terms
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

    return mi_vals, E0

def mi_cv(mi_values):
    mi_pos = [m for m in mi_values if m > 1e-10]
    if len(mi_pos) < 2: return 0.0
    return float(np.std(mi_pos) / np.mean(mi_pos))

print("=== Sprint 036b: XXZ MI-CV Size Scaling ===\n")

delta_values = [-0.5, 0.0, 0.5, 0.8, 0.9, 1.0, 1.1, 1.3, 1.5, 2.0]
results = {'experiment': '036b', 'data': {}}

# n=8 exact (fast sanity check)
print("--- n=8 (exact) ---")
results['data']['8'] = {}
for delta in delta_values:
    psi = xxz_ground_exact(8, delta)
    mi = all_pairs_mi_exact(psi, 8)
    cv = mi_cv(mi)
    mean_mi = float(np.mean([m for m in mi if m > 1e-10])) if any(m > 1e-10 for m in mi) else 0
    print(f"  Δ={delta:5.1f}: CV={cv:.4f}, mean={mean_mi:.6f}")
    results['data']['8'][str(delta)] = {'delta': delta, 'mi_cv': float(cv), 'mean_mi': mean_mi}

with open('results/sprint_036b_xxz_micv.json', 'w') as f:
    json.dump(results, f, indent=2)

# n=8 DMRG cross-check at one point
print("\n--- n=8 DMRG cross-check at Δ=1.0 ---")
mi_dmrg, _ = dmrg_all_mi_xxz(8, 1.0, chi_max=100)
cv_dmrg = mi_cv(mi_dmrg)
cv_exact = results['data']['8']['1.0']['mi_cv']
print(f"  Exact CV={cv_exact:.4f}, DMRG CV={cv_dmrg:.4f}, diff={abs(cv_exact-cv_dmrg):.6f}")

# n=16 DMRG
elapsed = time.time() - t0
print(f"\n--- n=16 (chi=100), budget={55-elapsed:.0f}s ---")
results['data']['16'] = {}
for delta in delta_values:
    elapsed = time.time() - t0
    if elapsed > 48: break
    t1 = time.time()
    try:
        mi, _ = dmrg_all_mi_xxz(16, delta, chi_max=100)
        cv = mi_cv(mi)
        mean_mi = float(np.mean([m for m in mi if m > 1e-10])) if any(m > 1e-10 for m in mi) else 0
        dt = time.time() - t1
        print(f"  Δ={delta:5.1f}: CV={cv:.4f}, mean={mean_mi:.6f} ({dt:.1f}s)")
        results['data']['16'][str(delta)] = {'delta': delta, 'mi_cv': float(cv), 'mean_mi': mean_mi, 'time': float(dt)}
    except Exception as e:
        print(f"  Δ={delta:5.1f}: FAILED: {e}")

    with open('results/sprint_036b_xxz_micv.json', 'w') as f:
        json.dump(results, f, indent=2)

# Summary
print(f"\n{'='*60}")
print("XXZ MI-CV vs system size:")
avail = [n for n in ['8', '16'] if n in results['data'] and results['data'][n]]
print(f"{'Δ':>6}", end='')
for n in avail: print(f"  n={n:>3}", end='')
print()
for delta in delta_values:
    print(f"{delta:6.1f}", end='')
    for n in avail:
        d = results['data'][n].get(str(delta), {})
        cv = d.get('mi_cv', float('nan'))
        print(f"  {cv:6.3f}" if not np.isnan(cv) else f"  {'---':>6}", end='')
    print()

# BKT dome analysis
print(f"\nBKT dome (peak CV location and width):")
for n_str in avail:
    data = results['data'][n_str]
    deltas = sorted([float(d) for d in data.keys()])
    cvs = [data[str(d)]['mi_cv'] for d in deltas]
    peak_idx = np.argmax(cvs)
    peak_d = deltas[peak_idx]
    peak_cv = cvs[peak_idx]
    half_max = peak_cv / 2
    above = [d for d, c in zip(deltas, cvs) if c > half_max]
    width = above[-1] - above[0] if above else 0
    print(f"  n={n_str}: peak Δ={peak_d:.1f}, CV={peak_cv:.4f}, FWHM≈{width:.1f}")

results['total_runtime'] = time.time() - t0
with open('results/sprint_036b_xxz_micv.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nRuntime: {time.time()-t0:.1f}s")
