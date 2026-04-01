"""
Sprint 038b: q=3 Potts MI-CV transition signature.

The 1D q=3 Potts (clock) model has a quantum phase transition.
Does MI-CV show crossing curves (like Ising) or a different signature?

Uses TeNPy ClockChain with q=3. Reconstructs all-pairs MI from
qutrit correlation functions (generalized Gell-Mann basis).
"""

import numpy as np
from scipy import linalg as la
import json, time, warnings
warnings.filterwarnings('ignore')
import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)
from tenpy.models.clock import ClockChain
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS
from tenpy.linalg.np_conserved import Array as npc_Array

t0 = time.time()

# ========== Qutrit operator basis (Gell-Mann + identity) ==========
d = 3
# Define all 9 basis operators for d=3 single-site Hilbert space
# Standard Gell-Mann matrices (λ_1 through λ_8) plus identity
gm = []  # List of (name, matrix) tuples

# Identity
gm.append(('I3', np.eye(3, dtype=complex)))

# λ_1: |0><1| + |1><0|
m = np.zeros((3,3), dtype=complex); m[0,1] = 1; m[1,0] = 1
gm.append(('gm1', m.copy()))

# λ_2: -i|0><1| + i|1><0|
m = np.zeros((3,3), dtype=complex); m[0,1] = -1j; m[1,0] = 1j
gm.append(('gm2', m.copy()))

# λ_3: |0><0| - |1><1|
m = np.zeros((3,3), dtype=complex); m[0,0] = 1; m[1,1] = -1
gm.append(('gm3', m.copy()))

# λ_4: |0><2| + |2><0|
m = np.zeros((3,3), dtype=complex); m[0,2] = 1; m[2,0] = 1
gm.append(('gm4', m.copy()))

# λ_5: -i|0><2| + i|2><0|
m = np.zeros((3,3), dtype=complex); m[0,2] = -1j; m[2,0] = 1j
gm.append(('gm5', m.copy()))

# λ_6: |1><2| + |2><1|
m = np.zeros((3,3), dtype=complex); m[1,2] = 1; m[2,1] = 1
gm.append(('gm6', m.copy()))

# λ_7: -i|1><2| + i|2><1|
m = np.zeros((3,3), dtype=complex); m[1,2] = -1j; m[2,1] = 1j
gm.append(('gm7', m.copy()))

# λ_8: (|0><0| + |1><1| - 2|2><2|) / sqrt(3)
m = np.zeros((3,3), dtype=complex); m[0,0] = 1; m[1,1] = 1; m[2,2] = -2
m /= np.sqrt(3)
gm.append(('gm8', m.copy()))

# Verify orthogonality: Tr(λ_a λ_b) = 2δ_{ab} for a,b >= 1; Tr(I·I) = 3
for i, (ni, mi) in enumerate(gm):
    for j, (nj, mj) in enumerate(gm):
        tr = np.trace(mi.conj().T @ mj)
        if i == j:
            expected = 3 if i == 0 else 2
            assert abs(tr - expected) < 1e-10, f"Tr({ni}·{nj}) = {tr}, expected {expected}"
        else:
            assert abs(tr) < 1e-10, f"Tr({ni}·{nj}) = {tr}, expected 0"
print("Gell-Mann basis verified orthogonal.\n")

def entropy(rho):
    ev = la.eigvalsh(rho)
    ev = ev[ev > 1e-15]
    return float(-np.sum(ev * np.log2(ev)))

def mi_from_rho2(rho_ij, d=3):
    """Mutual information from two-qutrit density matrix."""
    rho_ij = rho_ij.reshape(d, d, d, d)
    rho_i = np.trace(rho_ij, axis1=1, axis2=3)
    rho_j = np.trace(rho_ij, axis1=0, axis2=2)
    return entropy(rho_i) + entropy(rho_j) - entropy(rho_ij.reshape(d*d, d*d))

def mi_cv(mi_values):
    mi_pos = [m for m in mi_values if m > 1e-10]
    if len(mi_pos) < 2: return 0.0
    return float(np.std(mi_pos) / np.mean(mi_pos))

def potts_mi_cv(n, g_J, chi_max=80):
    """Compute MI-CV for q=3 Potts/clock chain.

    H = -J Σ X_i X†_j + h.c.  - g Σ Z_i + h.c.
    J=1 fixed, vary g. Transition at g/J ≈ 1 (self-dual point).
    """
    model = ClockChain({'L': n, 'q': 3, 'J': 1.0, 'g': g_J,
                         'bc_MPS': 'finite', 'conserve': None})
    site = model.lat.site(0)

    # Register Gell-Mann operators on the site
    for name, mat in gm[1:]:  # skip identity, already have 'Id'
        if name not in site.opnames:
            site.add_op(name, mat)

    psi = MPS.from_product_state(model.lat.mps_sites(), [0]*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-10,
        'trunc_params': {'chi_max': chi_max}, 'max_sweeps': 20,
    })
    eng.run()

    # Compute all single-site expectations and two-site correlators
    op_names = ['Id'] + [name for name, _ in gm[1:]]
    op_mats = [np.eye(3, dtype=complex)] + [mat for _, mat in gm[1:]]

    # Single-site expectations
    exp_vals = {}
    for name in op_names[1:]:  # skip identity (always 1)
        exp_vals[name] = psi.expectation_value(name)

    # Two-site correlators: <O_a(i) O_b(j)>
    corr = {}
    for a_name in op_names[1:]:
        for b_name in op_names[1:]:
            corr[(a_name, b_name)] = psi.correlation_function(a_name, b_name)

    # Reconstruct ρ_ij for all pairs and compute MI
    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            # ρ_ij = (1/d²)[I⊗I + (d/2)Σ_a <λ_a>_i (λ_a⊗I) + (d/2)Σ_a <λ_a>_j (I⊗λ_a)
            #         + (d²/4)Σ_{a,b} <λ_a(i)λ_b(j)> (λ_a⊗λ_b)]
            # For Gell-Mann with Tr(λ_a λ_b) = 2δ_{ab}:
            rho_ij = np.eye(d*d, dtype=complex) / (d*d)

            for idx_a, a_name in enumerate(op_names[1:]):
                a_mat = op_mats[idx_a + 1]
                # Single-site terms
                ev_i = complex(exp_vals[a_name][i])
                ev_j = complex(exp_vals[a_name][j])
                rho_ij += (ev_i / (2*d)) * np.kron(a_mat, np.eye(d))
                rho_ij += (ev_j / (2*d)) * np.kron(np.eye(d), a_mat)

                # Two-site correlator terms
                for idx_b, b_name in enumerate(op_names[1:]):
                    b_mat = op_mats[idx_b + 1]
                    c_val = complex(corr[(a_name, b_name)][i, j])
                    rho_ij += (c_val / 4.0) * np.kron(a_mat, b_mat)

            # Ensure Hermitian and positive
            rho_ij = (rho_ij + rho_ij.conj().T) / 2
            ev = la.eigvalsh(rho_ij)
            if np.min(ev) < -1e-10:
                ev_pos = np.maximum(ev, 0)
                evec = la.eigh(rho_ij)[1]
                rho_ij = evec @ np.diag(ev_pos) @ evec.T
                rho_ij /= np.trace(rho_ij)

            mi = mi_from_rho2(rho_ij, d=d)
            mi_vals.append(max(float(np.real(mi)), 0))

    return mi_cv(mi_vals)


# ========== Main experiment ==========
print("=== Sprint 038b: q=3 Potts MI-CV ===\n")

# First: timing calibration
print("Calibrating n=8...")
t1 = time.time()
cv_test = potts_mi_cv(8, 1.0, chi_max=40)
dt = time.time() - t1
print(f"  n=8, g=1.0: CV={cv_test:.4f} ({dt:.1f}s)\n")

# Sweep g/J for multiple sizes
# Potts self-dual point is at g/J = 1
g_vals = [0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.3, 1.5, 2.0]
sizes = [8, 12, 16]
chi_map = {8: 40, 12: 60, 16: 80}

results = {'experiment': '038b', 'model': 'q=3 Potts (clock)',
           'sizes': {}, 'transition_point': 'g/J ~ 1 (self-dual)'}

for n in sizes:
    elapsed = time.time() - t0
    if elapsed > 50:
        print(f"TIMEOUT at {elapsed:.0f}s, stopping sizes")
        break

    chi = chi_map[n]
    print(f"--- n={n} (chi={chi}) ---")
    results['sizes'][str(n)] = {}

    for g_J in g_vals:
        elapsed = time.time() - t0
        if elapsed > 52:
            print(f"  TIMEOUT at {elapsed:.0f}s")
            break

        t1 = time.time()
        try:
            cv = potts_mi_cv(n, g_J, chi_max=chi)
            dt = time.time() - t1
            print(f"  g/J={g_J:.2f}: CV={cv:.4f} ({dt:.1f}s)")
            results['sizes'][str(n)][str(g_J)] = {'cv': float(cv), 'time': float(dt)}
        except Exception as e:
            print(f"  g/J={g_J:.2f}: FAILED: {e}")
            results['sizes'][str(n)][str(g_J)] = {'cv': None, 'error': str(e)}

    # Save after each size
    results['total_runtime'] = time.time() - t0
    with open('results/sprint_038b_potts_micv.json', 'w') as f:
        json.dump(results, f, indent=2)

print(f"\nTotal runtime: {time.time()-t0:.1f}s")
print("Results saved to results/sprint_038b_potts_micv.json")
