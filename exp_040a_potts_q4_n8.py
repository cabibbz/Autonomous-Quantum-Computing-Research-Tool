"""
Sprint 040a: q=4 Potts MI-CV at n=8.
q=4 is marginal (BKT-like logarithmic corrections at transition).
Does MI-CV show crossing curves or dome signature?

Uses generalized Gell-Mann matrices for d=4 (15 generators of SU(4)).
First: timing test at g=1.0, then sweep g values.
"""

import numpy as np
from scipy import linalg as la
import json, time, warnings, sys
warnings.filterwarnings('ignore')
import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)
from tenpy.models.clock import ClockChain
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS

t0 = time.time()
d = 4  # q=4 Potts

# ========== Generalized Gell-Mann matrices for SU(4) ==========
# d^2 - 1 = 15 generators. Three types:
# 1. Symmetric off-diagonal: (|j><k| + |k><j|) for j<k  -> d(d-1)/2 = 6
# 2. Antisymmetric off-diagonal: -i(|j><k| - |k><j|) for j<k  -> 6
# 3. Diagonal: generalized diagonal matrices -> d-1 = 3

gm_mats = []

# Type 1: Symmetric off-diagonal
for j in range(d):
    for k in range(j+1, d):
        m = np.zeros((d, d), dtype=complex)
        m[j, k] = 1; m[k, j] = 1
        gm_mats.append((f'sym_{j}{k}', m.copy()))

# Type 2: Antisymmetric off-diagonal
for j in range(d):
    for k in range(j+1, d):
        m = np.zeros((d, d), dtype=complex)
        m[j, k] = -1j; m[k, j] = 1j
        gm_mats.append((f'asym_{j}{k}', m.copy()))

# Type 3: Diagonal
for l in range(1, d):
    m = np.zeros((d, d), dtype=complex)
    norm = np.sqrt(2.0 / (l * (l + 1)))
    for j in range(l):
        m[j, j] = norm
    m[l, l] = -l * norm
    gm_mats.append((f'diag_{l}', m.copy()))

assert len(gm_mats) == d*d - 1, f"Expected {d*d-1} generators, got {len(gm_mats)}"

# Verify orthogonality: Tr(λ_a λ_b) = 2δ_{ab}
for i, (ni, mi_mat) in enumerate(gm_mats):
    for j, (nj, mj_mat) in enumerate(gm_mats):
        tr = np.real(np.trace(mi_mat.conj().T @ mj_mat))
        if i == j:
            assert abs(tr - 2.0) < 1e-10, f"Tr({ni}^2) = {tr}, expected 2"
        else:
            assert abs(tr) < 1e-10, f"Tr({ni}·{nj}) = {tr}, expected 0"
print(f"SU({d}) basis verified: {len(gm_mats)} generators, all orthogonal.\n")

op_names = [name for name, _ in gm_mats]
op_mats_all = [np.eye(d, dtype=complex)] + [mat for _, mat in gm_mats]


def entropy(rho):
    ev = la.eigvalsh(rho); ev = ev[ev > 1e-15]
    return float(-np.sum(ev * np.log2(ev)))


def mi_cv_q4(n, g_J, chi_max=40):
    """Compute MI-CV for q=4 Potts/clock chain."""
    model = ClockChain({'L': n, 'q': 4, 'J': 1.0, 'g': g_J,
                         'bc_MPS': 'finite', 'conserve': None})
    site = model.lat.site(0)
    for name, mat in gm_mats:
        if name not in site.opnames:
            site.add_op(name, mat)

    psi = MPS.from_product_state(model.lat.mps_sites(), [0]*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-8,
        'trunc_params': {'chi_max': chi_max}, 'max_sweeps': 10,
    })
    eng.run()

    # Single-site expectations
    exp_vals = {name: psi.expectation_value(name) for name in op_names}
    # Two-site correlators
    corr = {}
    for a in op_names:
        for b in op_names:
            corr[(a, b)] = psi.correlation_function(a, b)

    # Reconstruct ρ_ij and compute MI for all pairs
    mi_vals = []
    for i in range(n):
        for j in range(i+1, n):
            rho_ij = np.eye(d*d, dtype=complex) / (d*d)
            for idx_a, a_name in enumerate(op_names):
                a_mat = op_mats_all[idx_a + 1]
                ev_i = complex(exp_vals[a_name][i])
                ev_j = complex(exp_vals[a_name][j])
                rho_ij += (ev_i / (2*d)) * np.kron(a_mat, np.eye(d))
                rho_ij += (ev_j / (2*d)) * np.kron(np.eye(d), a_mat)
                for idx_b, b_name in enumerate(op_names):
                    b_mat = op_mats_all[idx_b + 1]
                    c_val = complex(corr[(a_name, b_name)][i, j])
                    rho_ij += (c_val / 4.0) * np.kron(a_mat, b_mat)

            rho_ij = (rho_ij + rho_ij.conj().T) / 2
            ev = la.eigvalsh(rho_ij)
            if np.min(ev) < -1e-10:
                ev_pos = np.maximum(ev, 0)
                evec = la.eigh(rho_ij)[1]
                rho_ij = evec @ np.diag(ev_pos) @ evec.T
                rho_ij /= np.trace(rho_ij)

            rho_d = rho_ij.reshape(d, d, d, d)
            rho_i = np.trace(rho_d, axis1=1, axis2=3)
            rho_j = np.trace(rho_d, axis1=0, axis2=2)
            mi = entropy(rho_i) + entropy(rho_j) - entropy(rho_ij)
            mi_vals.append(max(float(np.real(mi)), 0))

    mi_pos = [m for m in mi_vals if m > 1e-10]
    if len(mi_pos) < 2: return 0.0
    return float(np.std(mi_pos) / np.mean(mi_pos))


# ========== Main ==========
print("=== Sprint 040a: q=4 Potts MI-CV at n=8 ===\n")
sys.stdout.flush()

# Timing test first
print("Timing test: n=8, g=1.0, chi=40...")
sys.stdout.flush()
t1 = time.time()
cv_test = mi_cv_q4(8, 1.0, chi_max=40)
dt_test = time.time() - t1
print(f"  CV={cv_test:.4f} ({dt_test:.1f}s)\n")
sys.stdout.flush()

results = {'experiment': '040a', 'model': 'q=4 Potts (clock)', 'n': 8,
           'timing_test': {'g': 1.0, 'cv': cv_test, 'time': dt_test},
           'data': {}}

# Save timing test immediately
with open('results/sprint_040a_potts_q4_n8.json', 'w') as f:
    json.dump(results, f, indent=2)

# Sweep g values — same range as q=3 for comparison
g_vals = [0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.3, 1.5]

for g_J in g_vals:
    elapsed = time.time() - t0
    if elapsed > 52:
        print(f"  TIMEOUT at {elapsed:.0f}s")
        break
    t1 = time.time()
    try:
        cv = mi_cv_q4(8, g_J, chi_max=40)
        dt = time.time() - t1
        print(f"  g/J={g_J:.2f}: CV={cv:.4f} ({dt:.1f}s)")
        results['data'][str(g_J)] = {'cv': float(cv), 'time': float(dt)}
    except Exception as e:
        print(f"  g/J={g_J:.2f}: FAILED: {e}")
        results['data'][str(g_J)] = {'cv': None, 'error': str(e)}
    sys.stdout.flush()
    # Save incrementally
    results['total_runtime'] = time.time() - t0
    with open('results/sprint_040a_potts_q4_n8.json', 'w') as f:
        json.dump(results, f, indent=2)

print(f"\nTotal runtime: {time.time()-t0:.1f}s")
print(f"Points computed: {len(results['data'])}")
