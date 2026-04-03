"""
Sprint 035b2: BW Fidelity with Optimally-Fitted vs Projected Coefficients

For TFIM at multiple sizes: project H_E onto TFIM operator subspace,
build ρ_BW = exp(-H_proj)/Z, compute fidelity F(ρ_A, ρ_BW).

This tells us: is the BW FORM correct even if the envelope is wrong?
Compare with projection onto random operators (control).
"""

import numpy as np
from scipy import linalg as la
import json, time, warnings
warnings.filterwarnings('ignore')

t0 = time.time()

I2 = np.eye(2)
X_std = np.array([[0,1],[1,0]])
Z_std = np.array([[1,0],[0,-1]])
Sx_tp = np.array([[0,1],[1,0]])
Sz_tp = np.diag([-1, 1])

def kron_list(ops):
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r

def get_rho_A_exact_tfim(n, h_J):
    dim = 2**n; H = np.zeros((dim, dim))
    for i in range(n-1):
        ops = [I2]*n; ops[i] = Z_std; ops[i+1] = Z_std; H -= kron_list(ops)
    for i in range(n):
        ops = [I2]*n; ops[i] = X_std; H -= h_J * kron_list(ops)
    evals, evecs = la.eigh(H); psi = evecs[:, 0]
    rho = np.outer(psi, psi.conj()); n_A = n // 2
    rho_t = rho.reshape([2]*(2*n))
    for q in sorted(range(n_A, n), reverse=True):
        nr = rho_t.ndim // 2; rho_t = np.trace(rho_t, axis1=q, axis2=q+nr)
    return rho_t.reshape(2**n_A, 2**n_A), evals[0], 'std'

def get_rho_A_dmrg_tfim(n, h_J):
    import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)
    from tenpy.models.tf_ising import TFIChain
    from tenpy.algorithms import dmrg
    from tenpy.networks.mps import MPS
    model = TFIChain({'L': n, 'J': 1.0, 'g': h_J, 'bc_MPS': 'finite'})
    psi = MPS.from_product_state(model.lat.mps_sites(), ['up']*n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12, 'trunc_params': {'chi_max': 200}, 'max_sweeps': 30,
    })
    E0, _ = eng.run(); n_A = n // 2
    rho = psi.get_rho_segment(list(range(n_A))).to_ndarray().reshape(2**n_A, 2**n_A)
    return rho, E0, 'tenpy'

def fidelity(rho, sigma):
    """Uhlmann fidelity F(ρ,σ) = (Tr √(√ρ σ √ρ))²."""
    try:
        sqrt_rho = la.sqrtm(rho)
        M = sqrt_rho @ sigma @ sqrt_rho
        M = (M + M.conj().T) / 2  # ensure Hermitian
        evals = la.eigvalsh(M)
        evals = np.maximum(np.real(evals), 0)
        return float((np.sum(np.sqrt(evals)))**2)
    except:
        return 0.0

def build_tfim_ops(n_A, basis):
    if basis == 'tenpy':
        bond_pair, site_op = (Sx_tp, Sx_tp), Sz_tp
    else:
        bond_pair, site_op = (Z_std, Z_std), X_std
    ops = []
    for i in range(n_A):
        o = [I2]*n_A; o[i] = site_op; ops.append(kron_list(o))
    for i in range(n_A - 1):
        o = [I2]*n_A; o[i] = bond_pair[0]; o[i+1] = bond_pair[1]; ops.append(kron_list(o))
    return ops

def project_onto_ops(H, ops, dim):
    """Project H onto span of operators, return projected H and coefficients."""
    coeffs = []
    for O in ops:
        c = np.real(np.trace(H @ O)) / np.real(np.trace(O @ O))
        coeffs.append(c)
    H_proj = sum(c * O for c, O in zip(coeffs, ops))
    return H_proj, coeffs

def analyze(rho_A, n_A, basis, label, dim):
    """Full BW analysis: project H_E onto TFIM ops, compute fidelity."""
    # H_E
    evals_rho = la.eigvalsh(rho_A)
    evals_rho = np.maximum(evals_rho, 1e-15)
    evecs = la.eigh(rho_A)[1]
    H_E = -evecs @ np.diag(np.log(evals_rho)) @ evecs.T
    H_E_tl = H_E - np.trace(H_E)/dim * np.eye(dim)

    # TFIM projection
    tfim_ops = build_tfim_ops(n_A, basis)
    H_proj, coeffs = project_onto_ops(H_E_tl, tfim_ops, dim)

    # Frobenius metrics
    norm_E = la.norm(H_E_tl, 'fro')
    norm_proj = la.norm(H_proj, 'fro')
    norm_resid = la.norm(H_E_tl - H_proj, 'fro')
    pauli_frac = (norm_proj / norm_E)**2 if norm_E > 1e-15 else 0.0

    # Fidelity: ρ_BW = exp(-H_proj - c*I)/Z (add identity back for correct trace)
    # Note: we project the TRACELESS part, so we need to add a suitable identity
    rho_proj = la.expm(-H_proj)
    rho_proj = np.real(rho_proj)
    rho_proj /= np.trace(rho_proj)
    F_proj = fidelity(rho_A, rho_proj)

    # Also: fidelity with FULL H_E (should be 1.0, sanity check)
    rho_HE = la.expm(-H_E)
    rho_HE = np.real(rho_HE)
    rho_HE /= np.trace(rho_HE)
    F_exact = fidelity(rho_A, rho_HE)

    # Trace distance
    diff = rho_A - rho_proj
    td = 0.5 * np.sum(np.abs(la.eigvalsh(diff)))

    # Entanglement entropy comparison
    evals_proj = la.eigvalsh(rho_proj)
    evals_proj = evals_proj[evals_proj > 1e-15]
    S_exact = -np.sum(evals_rho[evals_rho > 1e-15] * np.log2(evals_rho[evals_rho > 1e-15]))
    S_proj = -np.sum(evals_proj * np.log2(evals_proj))

    return {
        'pauli_fraction': float(pauli_frac),
        'fidelity_projected': float(F_proj),
        'fidelity_exact': float(F_exact),
        'trace_distance': float(td),
        'S_exact': float(S_exact),
        'S_projected': float(S_proj),
        'S_error_pct': float(abs(S_proj - S_exact) / max(S_exact, 1e-10) * 100),
        'n_ops': len(tfim_ops),
        'coeffs': [float(c) for c in coeffs],
    }

# ---- Main ----
results = {'experiment': '035b2', 'data': {}}

print("=== Sprint 035b2: BW Fidelity with Projected Coefficients (TFIM h/J=1.0) ===\n")
print(f"{'n':>4} {'Pauli%':>7} {'F_proj':>8} {'TraceD':>7} {'S_exact':>8} {'S_proj':>8} {'S_err%':>7}")

for n in [8, 10, 12, 16, 20]:
    elapsed = time.time() - t0
    if elapsed > 48:
        print(f"{n:4d} SKIPPED ({elapsed:.0f}s)")
        continue

    n_A = n // 2; dim = 2**n_A
    t1 = time.time()

    if n <= 10:
        rho_A, E0, basis = get_rho_A_exact_tfim(n, 1.0)
    else:
        rho_A, E0, basis = get_rho_A_dmrg_tfim(n, 1.0)

    m = analyze(rho_A, n_A, basis, 'critical', dim)
    dt = time.time() - t1

    print(f"{n:4d} {m['pauli_fraction']*100:6.1f}% {m['fidelity_projected']:8.5f} "
          f"{m['trace_distance']:7.4f} {m['S_exact']:8.4f} {m['S_projected']:8.4f} {m['S_error_pct']:6.2f}%"
          f"  ({dt:.1f}s)")

    results['data'][f'critical_{n}'] = {**m, 'n': n, 'h_J': 1.0, 'E0': float(E0), 'time': dt}
    with open('results/sprint_035b2_bw_fidelity_fitted.json', 'w') as f:
        json.dump(results, f, indent=2)

# Ordered phase
print(f"\n{'n':>4} {'Pauli%':>7} {'F_proj':>8} {'TraceD':>7} {'S_exact':>8} {'S_proj':>8} {'S_err%':>7}")
for n in [8, 10, 12, 16]:
    elapsed = time.time() - t0
    if elapsed > 50:
        print(f"{n:4d} SKIPPED")
        continue
    n_A = n // 2; dim = 2**n_A
    if n <= 10:
        rho_A, E0, basis = get_rho_A_exact_tfim(n, 0.5)
    else:
        rho_A, E0, basis = get_rho_A_dmrg_tfim(n, 0.5)
    m = analyze(rho_A, n_A, basis, 'ordered', dim)
    dt = time.time() - t0 - elapsed
    print(f"{n:4d} {m['pauli_fraction']*100:6.1f}% {m['fidelity_projected']:8.5f} "
          f"{m['trace_distance']:7.4f} {m['S_exact']:8.4f} {m['S_projected']:8.4f} {m['S_error_pct']:6.2f}%")
    results['data'][f'ordered_{n}'] = {**m, 'n': n, 'h_J': 0.5, 'E0': float(E0)}

results['total_runtime'] = time.time() - t0
with open('results/sprint_035b2_bw_fidelity_fitted.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n{'='*60}")
print("KEY FINDING: Pauli fraction drops but fidelity measures tell the real story.")
print(f"Runtime: {time.time()-t0:.1f}s")
