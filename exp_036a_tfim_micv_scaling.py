"""
Sprint 036a: TFIM MI-CV Finite-Size Scaling (n=8 to 50)

Use correlation functions from MPS to reconstruct ALL pairwise 2-site RDMs.
rho_ij = (1/4) sum_{a,b} <sigma_a(i) sigma_b(j)> sigma_a x sigma_b
MI from eigenvalues of rho_ij.

This is O(n^2 * chi^2) via correlation_function(), much faster than
get_rho_segment for each pair.
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
    """MI from 2-qubit density matrix."""
    rho_i = np.trace(rho_ij.reshape(2,2,2,2), axis1=1, axis2=3)
    rho_j = np.trace(rho_ij.reshape(2,2,2,2), axis1=0, axis2=2)
    return entropy(rho_i) + entropy(rho_j) - entropy(rho_ij)

def tfim_ground_exact(n, h_J):
    dim = 2**n; H = np.zeros((dim, dim))
    for i in range(n-1):
        ops = [I2]*n; ops[i] = sz; ops[i+1] = sz; H -= kron_list(ops)
    for i in range(n):
        ops = [I2]*n; ops[i] = sx; H -= h_J * kron_list(ops)
    evals, evecs = la.eigh(H)
    return evecs[:, 0]

def all_pairs_mi_exact(psi, n):
    """All-pairs MI from exact state vector."""
    rho = np.outer(psi, psi.conj()).reshape([2]*(2*n))
    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            keep = sorted([i, j])
            trace_q = [q for q in range(n) if q not in keep]
            rt = rho.copy()
            for q in sorted(trace_q, reverse=True):
                nr = rt.ndim // 2
                rt = np.trace(rt, axis1=q, axis2=q+nr)
            mi = mi_from_rho2(rt.reshape(4,4))
            mi_vals.append(max(mi, 0))
    return mi_vals

def dmrg_all_mi_tfim(n, h_J, chi_max=100):
    """All-pairs MI from DMRG using correlation functions."""
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

    # Single-site expectations
    exp_x = psi.expectation_value('Sigmax')  # <sigma_x(i)> for each site
    exp_y = psi.expectation_value('Sigmay')
    exp_z = psi.expectation_value('Sigmaz')

    # Two-point correlations: <sigma_a(i) sigma_b(j)> matrices
    # correlation_function returns n x n matrix
    ops = ['Sigmax', 'Sigmay', 'Sigmaz']
    corr = {}
    for a in ops:
        for b in ops:
            corr[(a,b)] = psi.correlation_function(a, b)

    # Reconstruct 2-site RDMs and compute MI
    pauli = {'Sigmax': sx, 'Sigmay': sy, 'Sigmaz': sz}
    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            # rho_ij = (1/4)(I4 + sum_a <a_i> a x I + sum_b <b_j> I x b + sum_{ab} <a_i b_j> a x b)
            rho_ij = np.eye(4, dtype=complex) / 4.0

            # Single-site terms
            for a_name, a_val, exp_a in [('x', sx, exp_x), ('y', sy, exp_y), ('z', sz, exp_z)]:
                rho_ij += float(np.real(exp_a[i])) * np.kron(pauli[f'Sigma{a_name}'], I2) / 4.0
                rho_ij += float(np.real(exp_a[j])) * np.kron(I2, pauli[f'Sigma{a_name}']) / 4.0

            # Two-point terms
            for a in ops:
                for b in ops:
                    c_val = float(np.real(corr[(a,b)][i, j]))
                    rho_ij += c_val * np.kron(pauli[a], pauli[b]) / 4.0

            rho_ij = np.real(rho_ij)  # Should be real for physical states

            # Ensure positive semi-definite (numerical noise)
            ev = la.eigvalsh(rho_ij)
            if np.min(ev) < -1e-10:
                # Project to PSD
                ev_pos = np.maximum(ev, 0)
                evec = la.eigh(rho_ij)[1]
                rho_ij = evec @ np.diag(ev_pos) @ evec.T
                rho_ij /= np.trace(rho_ij)

            mi = mi_from_rho2(rho_ij)
            mi_vals.append(max(mi, 0))

    return mi_vals, E0

def mi_cv(mi_values):
    mi_pos = [m for m in mi_values if m > 1e-10]
    if len(mi_pos) < 2: return 0.0
    return float(np.std(mi_pos) / np.mean(mi_pos))

print("=== Sprint 036a: TFIM MI-CV All-Pairs Scaling ===\n")

# Timing test at n=32
print("--- Timing: n=32, chi=100, h/J=1.0 ---")
t1 = time.time()
mi_t, _ = dmrg_all_mi_tfim(32, 1.0, chi_max=100)
dt32 = time.time() - t1
print(f"  n=32: {len(mi_t)} pairs, CV={mi_cv(mi_t):.4f}, time={dt32:.1f}s")

results = {'experiment': '036a', 'description': 'TFIM MI-CV all-pairs via correlations', 'data': {}}
results['timing_n32'] = dt32

with open('results/sprint_036a_tfim_micv.json', 'w') as f:
    json.dump(results, f, indent=2)

# Budget and size selection
budget = 55 - (time.time() - t0)
print(f"  Budget remaining: {budget:.0f}s")

if dt32 < 5:
    sizes = [8, 16, 32, 50]
    chi_map = {8: 100, 16: 100, 32: 100, 50: 64}
elif dt32 < 8:
    sizes = [8, 16, 32]
    chi_map = {8: 100, 16: 100, 32: 100}
else:
    sizes = [8, 16]
    chi_map = {8: 100, 16: 100}

# Dense near criticality
h_J_values = [0.5, 0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.3, 1.5, 2.0]

for n in sizes:
    elapsed = time.time() - t0
    if elapsed > 50:
        print(f"\n  n={n}: TIMEOUT ({elapsed:.0f}s)")
        break

    chi = chi_map.get(n, 100)
    print(f"\n--- n={n} (chi={chi}) ---")
    results['data'][str(n)] = {}

    for h_J in h_J_values:
        elapsed = time.time() - t0
        if elapsed > 50:
            print(f"  TIMEOUT at {elapsed:.0f}s")
            break

        t1 = time.time()
        try:
            if n <= 10:
                psi = tfim_ground_exact(n, h_J)
                mi = all_pairs_mi_exact(psi, n)
            else:
                mi, _ = dmrg_all_mi_tfim(n, h_J, chi)

            cv = mi_cv(mi)
            mean_mi = float(np.mean([m for m in mi if m > 1e-10])) if any(m > 1e-10 for m in mi) else 0
            dt = time.time() - t1
            print(f"  h/J={h_J:.2f}: CV={cv:.4f}, mean={mean_mi:.6f} ({dt:.1f}s)")

            results['data'][str(n)][str(h_J)] = {
                'h_J': h_J, 'mi_cv': float(cv), 'mean_mi': mean_mi,
                'n_pairs': len(mi), 'time': float(dt)
            }
        except Exception as e:
            print(f"  h/J={h_J:.2f}: FAILED: {e}")

    with open('results/sprint_036a_tfim_micv.json', 'w') as f:
        json.dump(results, f, indent=2)

# Summary
print(f"\n{'='*70}")
avail = [n for n in sizes if str(n) in results['data'] and results['data'][str(n)]]
print("TFIM MI-CV (all pairs) vs system size:")
print(f"{'h/J':>6}", end='')
for n in avail: print(f"  n={n:>3}", end='')
print()

all_h = sorted(set(float(h) for ns in results['data'].values() for h in ns.keys()))
for h_J in all_h:
    print(f"{h_J:6.2f}", end='')
    for n in avail:
        d = results['data'].get(str(n), {}).get(str(h_J), {})
        cv = d.get('mi_cv', float('nan'))
        print(f"  {cv:6.3f}" if not np.isnan(cv) else f"  {'---':>6}", end='')
    print()

# Slopes near criticality
print(f"\nTransition slope (dCV/d(h/J)) near h/J=1.0:")
slopes = {}
for n_str, data in sorted(results['data'].items(), key=lambda x: int(x[0])):
    n = int(n_str)
    h_vals = sorted([float(h) for h in data.keys()])
    cv_vals = [data[str(h)]['mi_cv'] for h in h_vals]
    if len(h_vals) >= 3:
        idx = min(range(len(h_vals)), key=lambda i: abs(h_vals[i] - 1.0))
        if 0 < idx < len(h_vals) - 1:
            slope = (cv_vals[idx+1] - cv_vals[idx-1]) / (h_vals[idx+1] - h_vals[idx-1])
            slopes[n] = slope
            print(f"  n={n:>3}: slope = {slope:.4f}")

if len(slopes) >= 2:
    ns = sorted(slopes.keys())
    print(f"\n  Slope scaling:")
    for i in range(1, len(ns)):
        ratio = slopes[ns[i]] / slopes[ns[i-1]] if slopes[ns[i-1]] != 0 else float('inf')
        print(f"    n={ns[i]}/n={ns[i-1]} = {ratio:.2f}")

results['slopes'] = {str(k): v for k, v in slopes.items()}
results['total_runtime'] = time.time() - t0
with open('results/sprint_036a_tfim_micv.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nTotal runtime: {time.time()-t0:.1f}s")
