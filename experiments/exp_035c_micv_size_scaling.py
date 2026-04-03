"""
Sprint 035c: MI-CV Size Scaling for TFIM and XXZ

Sprint 030 found MI-CV (coefficient of variation of pairwise MI)
is a phase transition order parameter. All results were at n=6-8.

Test at n=8,12,16,20: does MI-CV sharpen with system size?
Use DMRG to get 2-site RDMs efficiently.

TFIM: MI-CV jump at h/J≈1 (Ising transition)
XXZ: MI-CV dome at Δ≈1 (BKT transition)
"""

import numpy as np
from scipy import linalg as la
import json, time, warnings
warnings.filterwarnings('ignore')

t0 = time.time()

I2 = np.eye(2)
X = np.array([[0,1],[1,0]])
Z = np.array([[1,0],[0,-1]])

def kron_list(ops):
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r

# ---- Exact diag helpers ----
def tfim_ground_state(n, h_J):
    dim = 2**n; H = np.zeros((dim, dim))
    for i in range(n-1):
        ops = [I2]*n; ops[i] = Z; ops[i+1] = Z; H -= kron_list(ops)
    for i in range(n):
        ops = [I2]*n; ops[i] = X; H -= h_J * kron_list(ops)
    evals, evecs = la.eigh(H)
    return evecs[:, 0], evals[0]

def xxz_ground_state(n, delta):
    Sp = np.array([[0,1],[0,0]]); Sm = np.array([[0,0],[1,0]])
    dim = 2**n; H = np.zeros((dim, dim), dtype=complex)
    for i in range(n-1):
        ops_pm = [I2]*n; ops_pm[i] = Sp; ops_pm[i+1] = Sm
        ops_mp = [I2]*n; ops_mp[i] = Sm; ops_mp[i+1] = Sp
        ops_zz = [I2]*n; ops_zz[i] = Z/2; ops_zz[i+1] = Z/2
        H += 0.5*(kron_list(ops_pm) + kron_list(ops_mp)) + delta * kron_list(ops_zz)
    H = np.real(H) if np.allclose(H.imag, 0) else H
    evals, evecs = la.eigh(H)
    return evecs[:, 0], evals[0]

def two_qubit_rdm_exact(psi, n, i, j):
    """Get 2-qubit RDM for qubits i,j from state vector."""
    rho = np.outer(psi, psi.conj())
    rho_t = rho.reshape([2]*(2*n))
    keep = sorted([i, j])
    trace_q = [q for q in range(n) if q not in keep]
    for q in sorted(trace_q, reverse=True):
        nr = rho_t.ndim // 2
        rho_t = np.trace(rho_t, axis1=q, axis2=q+nr)
    return rho_t.reshape(4, 4)

def single_qubit_rdm_exact(psi, n, i):
    rho = np.outer(psi, psi.conj())
    rho_t = rho.reshape([2]*(2*n))
    trace_q = [q for q in range(n) if q != i]
    for q in sorted(trace_q, reverse=True):
        nr = rho_t.ndim // 2
        rho_t = np.trace(rho_t, axis1=q, axis2=q+nr)
    return rho_t.reshape(2, 2)

# ---- DMRG helpers ----
def dmrg_mi_tfim(n, h_J):
    """Get pairwise MI from DMRG for TFIM."""
    import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)
    from tenpy.models.tf_ising import TFIChain
    from tenpy.algorithms import dmrg
    from tenpy.networks.mps import MPS

    model = TFIChain({'L': n, 'J': 1.0, 'g': h_J, 'bc_MPS': 'finite'})
    psi = MPS.from_product_state(model.lat.mps_sites(), ['up']*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12, 'trunc_params': {'chi_max': 200}, 'max_sweeps': 30,
    })
    E0, _ = eng.run()

    # Compute MI from single and two-site reduced density matrices
    # For MPS, we can get expectation values efficiently
    # Single-site entropies
    S1 = []
    for i in range(n):
        rho_i = psi.get_rho_segment([i]).to_ndarray().reshape(2, 2)
        ev = la.eigvalsh(rho_i)
        ev = ev[ev > 1e-15]
        S1.append(-np.sum(ev * np.log2(ev)))

    # Two-site entropies (only nearest and next-nearest for speed)
    mi_values = []
    for i in range(n):
        for j in range(i+1, min(i+4, n)):  # up to 3 apart
            rho_ij = psi.get_rho_segment([i, j]).to_ndarray().reshape(4, 4)
            ev = la.eigvalsh(rho_ij)
            ev = ev[ev > 1e-15]
            S2 = -np.sum(ev * np.log2(ev))
            mi = S1[i] + S1[j] - S2
            mi_values.append(max(mi, 0))

    return mi_values, E0

def dmrg_mi_xxz(n, delta):
    import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)
    from tenpy.models.xxz_chain import XXZChain
    from tenpy.algorithms import dmrg
    from tenpy.networks.mps import MPS

    model = XXZChain({'L': n, 'Jxx': 1.0, 'Jz': delta, 'hz': 0.0, 'bc_MPS': 'finite'})
    psi = MPS.from_product_state(model.lat.mps_sites(), ['up', 'down']*(n//2), bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12, 'trunc_params': {'chi_max': 200}, 'max_sweeps': 30,
    })
    E0, _ = eng.run()

    S1 = []
    for i in range(n):
        rho_i = psi.get_rho_segment([i]).to_ndarray().reshape(2, 2)
        ev = la.eigvalsh(rho_i); ev = ev[ev > 1e-15]
        S1.append(-np.sum(ev * np.log2(ev)))

    mi_values = []
    for i in range(n):
        for j in range(i+1, min(i+4, n)):
            rho_ij = psi.get_rho_segment([i, j]).to_ndarray().reshape(4, 4)
            ev = la.eigvalsh(rho_ij); ev = ev[ev > 1e-15]
            S2 = -np.sum(ev * np.log2(ev))
            mi_values.append(max(S1[i] + S1[j] - S2, 0))

    return mi_values, E0

def exact_mi(psi, n):
    """All pairwise MI from exact state vector (for n<=10)."""
    S1 = []
    for i in range(n):
        rho_i = single_qubit_rdm_exact(psi, n, i)
        ev = la.eigvalsh(rho_i); ev = ev[ev > 1e-15]
        S1.append(-np.sum(ev * np.log2(ev)))

    mi_values = []
    for i in range(n):
        for j in range(i+1, n):
            rho_ij = two_qubit_rdm_exact(psi, n, i, j)
            ev = la.eigvalsh(rho_ij); ev = ev[ev > 1e-15]
            S2 = -np.sum(ev * np.log2(ev))
            mi_values.append(max(S1[i] + S1[j] - S2, 0))
    return mi_values

def mi_cv(mi_values):
    """MI coefficient of variation (CV = std/mean)."""
    mi_pos = [m for m in mi_values if m > 1e-10]
    if len(mi_pos) < 2:
        return 0.0
    return float(np.std(mi_pos) / np.mean(mi_pos))

# ---- Main: TFIM sweep ----
results = {'experiment': '035c', 'data': {}}

print("=== Sprint 035c: MI-CV Size Scaling ===\n")

# TFIM at several h/J values, multiple sizes
h_J_values = [0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 2.0]

print("--- TFIM ---")
results['data']['tfim'] = {}

for n in [8, 12, 16]:
    elapsed = time.time() - t0
    if elapsed > 45:
        print(f"  n={n}: SKIPPED ({elapsed:.0f}s)")
        continue

    print(f"\n  n={n}:")
    results['data']['tfim'][str(n)] = {}

    for h_J in h_J_values:
        elapsed = time.time() - t0
        if elapsed > 45:
            break
        t1 = time.time()
        try:
            if n <= 10:
                psi, E0 = tfim_ground_state(n, h_J)
                mi = exact_mi(psi, n)
            else:
                mi, E0 = dmrg_mi_tfim(n, h_J)
            cv = mi_cv(mi)
            mean_mi = np.mean([m for m in mi if m > 1e-10]) if any(m > 1e-10 for m in mi) else 0
            dt = time.time() - t1
            print(f"    h/J={h_J:.1f}: MI-CV={cv:.4f}, mean_MI={mean_mi:.4f} ({dt:.1f}s)")
            results['data']['tfim'][str(n)][str(h_J)] = {
                'h_J': h_J, 'mi_cv': float(cv), 'mean_mi': float(mean_mi),
                'n_pairs': len(mi), 'time': float(dt)
            }
        except Exception as e:
            print(f"    h/J={h_J:.1f}: FAILED: {e}")

    with open('results/sprint_035c_micv_scaling.json', 'w') as f:
        json.dump(results, f, indent=2)

# XXZ at several Δ values
delta_values = [-0.5, 0.0, 0.5, 0.8, 0.95, 1.0, 1.05, 1.2, 1.5, 2.0]

print("\n--- XXZ ---")
results['data']['xxz'] = {}

for n in [8, 12]:
    elapsed = time.time() - t0
    if elapsed > 48:
        print(f"  n={n}: SKIPPED ({elapsed:.0f}s)")
        continue

    print(f"\n  n={n}:")
    results['data']['xxz'][str(n)] = {}

    for delta in delta_values:
        elapsed = time.time() - t0
        if elapsed > 48:
            break
        t1 = time.time()
        try:
            if n <= 10:
                psi, E0 = xxz_ground_state(n, delta)
                mi = exact_mi(psi, n)
            else:
                mi, E0 = dmrg_mi_xxz(n, delta)
            cv = mi_cv(mi)
            mean_mi = np.mean([m for m in mi if m > 1e-10]) if any(m > 1e-10 for m in mi) else 0
            dt = time.time() - t1
            print(f"    Δ={delta:.2f}: MI-CV={cv:.4f}, mean_MI={mean_mi:.4f} ({dt:.1f}s)")
            results['data']['xxz'][str(n)][str(delta)] = {
                'delta': delta, 'mi_cv': float(cv), 'mean_mi': float(mean_mi),
                'n_pairs': len(mi), 'time': float(dt)
            }
        except Exception as e:
            print(f"    Δ={delta:.2f}: FAILED: {e}")

    with open('results/sprint_035c_micv_scaling.json', 'w') as f:
        json.dump(results, f, indent=2)

# Summary
print(f"\n{'='*60}")
print("TFIM MI-CV vs size:")
print(f"{'h/J':>5}", end='')
for n in sorted(results['data']['tfim'].keys(), key=int):
    print(f"  n={n:>3}", end='')
print()
for h_J in h_J_values:
    h_str = str(h_J)
    print(f"{h_J:5.1f}", end='')
    for n in sorted(results['data']['tfim'].keys(), key=int):
        d = results['data']['tfim'][n].get(h_str, {})
        cv = d.get('mi_cv', float('nan'))
        print(f"  {cv:6.3f}", end='')
    print()

print(f"\nXXZ MI-CV vs size:")
print(f"{'Δ':>5}", end='')
for n in sorted(results['data']['xxz'].keys(), key=int):
    print(f"  n={n:>3}", end='')
print()
for delta in delta_values:
    d_str = str(delta)
    print(f"{delta:5.2f}", end='')
    for n in sorted(results['data']['xxz'].keys(), key=int):
        d = results['data']['xxz'][n].get(d_str, {})
        cv = d.get('mi_cv', float('nan'))
        print(f"  {cv:6.3f}", end='')
    print()

results['total_runtime'] = time.time() - t0
with open('results/sprint_035c_micv_scaling.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nRuntime: {time.time()-t0:.1f}s")
