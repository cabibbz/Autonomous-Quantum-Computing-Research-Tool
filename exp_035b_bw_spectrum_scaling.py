"""
Sprint 035b: BW Entanglement Spectrum Accuracy at Multiple Sizes

Sprint 035a showed Pauli fraction drops with n: 91%→47%→8%.
This is a metric artifact: ||H_E||² grows exponentially from irrelevant eigenvalues.

Better metric: How well does BW predict the IMPORTANT entanglement spectrum?
- Compare top-k Schmidt values (λ_1,...,λ_k where Σλ_i > 0.99)
- Use entanglement energy ε_i = -log(λ_i)
- Fit BW coefficients optimally per term, compute spectrum-level R²

Also test: is the BW envelope shape correct? Fit independent β(i) per term
and compare with predicted sin_inv envelope.

Test TFIM (Z₂) and XXZ (U(1)) at multiple sizes.
"""

import numpy as np
from scipy import linalg as la
from scipy.optimize import minimize
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

# ---- Ground state: exact or DMRG ----
def get_rho_A_exact_tfim(n, h_J):
    dim = 2**n; H = np.zeros((dim, dim))
    for i in range(n-1):
        ops = [I2]*n; ops[i] = Z_std; ops[i+1] = Z_std; H -= kron_list(ops)
    for i in range(n):
        ops = [I2]*n; ops[i] = X_std; H -= h_J * kron_list(ops)
    evals, evecs = la.eigh(H)
    psi = evecs[:, 0]; rho = np.outer(psi, psi.conj())
    n_A = n // 2
    rho_t = rho.reshape([2]*(2*n))
    for q in sorted(range(n_A, n), reverse=True):
        nr = rho_t.ndim // 2
        rho_t = np.trace(rho_t, axis1=q, axis2=q+nr)
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
    E0, _ = eng.run()
    n_A = n // 2
    rho = psi.get_rho_segment(list(range(n_A))).to_ndarray().reshape(2**n_A, 2**n_A)
    return rho, E0, 'tenpy'

def get_rho_A_dmrg_xxz(n, delta):
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
    n_A = n // 2
    rho = psi.get_rho_segment(list(range(n_A))).to_ndarray().reshape(2**n_A, 2**n_A)
    return rho, E0, 'tenpy'

def get_rho_A_exact_xxz(n, delta):
    dim = 2**n; H = np.zeros((dim, dim), dtype=complex)
    Sp = np.array([[0,1],[0,0]]); Sm = np.array([[0,0],[1,0]])
    for i in range(n-1):
        # Jxx (S+S- + S-S+) + Jz SzSz
        ops_pm = [I2]*n; ops_pm[i] = Sp; ops_pm[i+1] = Sm
        ops_mp = [I2]*n; ops_mp[i] = Sm; ops_mp[i+1] = Sp
        ops_zz = [I2]*n; ops_zz[i] = Z_std/2; ops_zz[i+1] = Z_std/2
        H += 0.5*(kron_list(ops_pm) + kron_list(ops_mp)) + delta * kron_list(ops_zz)
    H = np.real(H) if np.allclose(H.imag, 0) else H
    evals, evecs = la.eigh(H)
    psi = evecs[:, 0]; rho = np.outer(psi, psi.conj())
    n_A = n // 2
    rho_t = rho.reshape([2]*(2*n))
    for q in sorted(range(n_A, n), reverse=True):
        nr = rho_t.ndim // 2
        rho_t = np.trace(rho_t, axis1=q, axis2=q+nr)
    return rho_t.reshape(2**n_A, 2**n_A), evals[0], 'std'

# ---- Build TFIM term operators in subsystem ----
def build_tfim_ops(n_A, basis):
    """Return list of (label, operator) for TFIM-type terms."""
    if basis == 'tenpy':
        bond_pair, site_op = (Sx_tp, Sx_tp), Sz_tp
    else:
        bond_pair, site_op = (Z_std, Z_std), X_std

    ops = []
    for i in range(n_A):
        o = [I2]*n_A; o[i] = site_op
        ops.append((f'field_{i}', kron_list(o)))
    for i in range(n_A - 1):
        o = [I2]*n_A; o[i] = bond_pair[0]; o[i+1] = bond_pair[1]
        ops.append((f'bond_{i}', kron_list(o)))
    return ops

def build_xxz_ops(n_A, basis):
    """Return list of (label, operator) for XXZ-type terms."""
    if basis == 'tenpy':
        # TeNPy XXZChain: H = Jxx(SxSx+SySy) + Jz SzSz = Jxx/2(S+S-+S-S+) + Jz SzSz
        # Operators in TeNPy basis: Sx, Sy, Sz are spin-1/2 operators
        Sx = np.array([[0,0.5],[0.5,0]])
        Sy = np.array([[0,-0.5j],[0.5j,0]])
        Sz = np.diag([-0.5, 0.5])
    else:
        Sx = X_std/2
        Sy = np.array([[0,-1j],[1j,0]])/2
        Sz = Z_std/2

    ops = []
    for i in range(n_A - 1):
        # XX + YY
        o_xx = [I2]*n_A; o_xx[i] = 2*Sx; o_xx[i+1] = 2*Sx
        o_yy = [I2]*n_A; o_yy[i] = 2*np.real(Sy @ Sy.conj().T) if np.iscomplexobj(Sy) else 2*Sy; o_yy[i+1] = o_yy[i]
        # Just use SxSx + SySy directly
        o_pm = [I2]*n_A; o_pm[i] = np.array([[0,1],[0,0]]); o_pm[i+1] = np.array([[0,0],[1,0]])
        o_mp = [I2]*n_A; o_mp[i] = np.array([[0,0],[1,0]]); o_mp[i+1] = np.array([[0,1],[0,0]])
        ops.append((f'xy_{i}', 0.5*(kron_list(o_pm) + kron_list(o_mp))))
        # ZZ
        o_zz = [I2]*n_A; o_zz[i] = Sz; o_zz[i+1] = Sz
        ops.append((f'zz_{i}', kron_list(o_zz)))
    return ops

# ---- Spectrum-based BW accuracy ----
def spectrum_bw_accuracy(rho_A, n_A, h_J_or_delta, n_total, basis, model_type='tfim'):
    """
    Fit BW Hamiltonian to match entanglement spectrum.
    Return: R² of entanglement energies, fitted β(i), spectrum comparison.
    """
    dim = 2**n_A

    # Entanglement spectrum
    evals_rho = la.eigvalsh(rho_A)
    evals_rho = np.sort(evals_rho)[::-1]  # descending

    # Keep significant eigenvalues (cumsum > threshold)
    cumsum = np.cumsum(evals_rho)
    k = max(2, int(np.searchsorted(cumsum, 0.9999)) + 1)
    k = min(k, len(evals_rho))

    # Entanglement energies for significant eigenvalues
    eps_exact = -np.log(np.maximum(evals_rho[:k], 1e-15))

    # Build model operators
    if model_type == 'tfim':
        term_ops = build_tfim_ops(n_A, basis)
    else:
        term_ops = build_xxz_ops(n_A, basis)

    n_terms = len(term_ops)

    # Build H_BW with variable coefficients: H_BW = Σ β_i O_i
    # Store operator matrices
    O_matrices = [op for _, op in term_ops]

    def make_H_BW(betas):
        H = np.zeros((dim, dim))
        for b, O in zip(betas, O_matrices):
            H += b * O
        return H

    # Fit: minimize || eps_exact - eps_BW ||² over β
    def objective(betas):
        H_BW = make_H_BW(betas)
        evals_BW = la.eigvalsh(H_BW)
        evals_BW = np.sort(evals_BW)  # ascending (entanglement energies)
        # Map: smallest eigenvalues of H_BW correspond to largest eigenvalues of ρ
        eps_BW = evals_BW[:k]  # lowest k eigenvalues
        # We need a shift (H_BW is defined up to constant)
        shift = np.mean(eps_exact) - np.mean(eps_BW)
        return np.sum((eps_exact - eps_BW - shift)**2)

    # Initialize with BW envelope
    L = n_total
    beta_init = []
    for label, _ in term_ops:
        if label.startswith('field'):
            i = int(label.split('_')[1])
            dist = n_A - 0.5 - i
            beta_init.append(-h_J_or_delta * np.sin(np.pi * dist / L))
        elif label.startswith('bond') or label.startswith('xy') or label.startswith('zz'):
            i = int(label.split('_')[-1])
            dist = n_A - 1 - i
            beta_init.append(-np.sin(np.pi * dist / L))
        else:
            beta_init.append(-1.0)

    result = minimize(objective, beta_init, method='Nelder-Mead',
                     options={'maxiter': 5000, 'xatol': 1e-8})
    betas_opt = result.x

    # Compute R² with optimal betas
    H_BW_opt = make_H_BW(betas_opt)
    evals_BW = np.sort(la.eigvalsh(H_BW_opt))
    eps_BW = evals_BW[:k]
    shift = np.mean(eps_exact) - np.mean(eps_BW)
    eps_BW_shifted = eps_BW + shift

    ss_res = np.sum((eps_exact - eps_BW_shifted)**2)
    ss_tot = np.sum((eps_exact - np.mean(eps_exact))**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 1e-15 else 1.0

    # Also compute R² with BW envelope (not fitted)
    H_BW_env = make_H_BW(beta_init)
    evals_env = np.sort(la.eigvalsh(H_BW_env))
    eps_env = evals_env[:k]
    shift_env = np.mean(eps_exact) - np.mean(eps_env)
    eps_env_shifted = eps_env + shift_env
    ss_res_env = np.sum((eps_exact - eps_env_shifted)**2)
    R2_env = 1 - ss_res_env / ss_tot if ss_tot > 1e-15 else 1.0

    # Extract fitted envelope shape
    fitted_betas = {label: float(b) for (label, _), b in zip(term_ops, betas_opt)}

    return {
        'R2_fitted': float(R2),
        'R2_envelope': float(R2_env),
        'n_significant': int(k),
        'eps_exact_top5': [float(e) for e in eps_exact[:5]],
        'eps_bw_top5': [float(e) for e in eps_BW_shifted[:5]],
        'fitted_betas': fitted_betas,
    }

# ---- Main ----
results = {'experiment': '035b', 'description': 'BW spectrum accuracy scaling', 'data': {}}

print("=== Sprint 035b: BW Entanglement Spectrum Accuracy ===\n")

# TFIM critical
print("--- TFIM h/J=1.0 (critical) ---")
results['data']['tfim_critical'] = {}
for n in [8, 10, 12, 16]:
    elapsed = time.time() - t0
    if elapsed > 45:
        print(f"  n={n}: SKIP ({elapsed:.0f}s)")
        continue
    n_A = n // 2
    t1 = time.time()
    try:
        if n <= 10:
            rho_A, E0, basis = get_rho_A_exact_tfim(n, 1.0)
        else:
            rho_A, E0, basis = get_rho_A_dmrg_tfim(n, 1.0)
        m = spectrum_bw_accuracy(rho_A, n_A, 1.0, n, basis, 'tfim')
        dt = time.time() - t1
        print(f"  n={n}: R²_fitted={m['R2_fitted']:.6f}, R²_env={m['R2_envelope']:.6f}, "
              f"k={m['n_significant']}, time={dt:.1f}s")
        print(f"    ε_exact: {[f'{e:.3f}' for e in m['eps_exact_top5']]}")
        print(f"    ε_BW:    {[f'{e:.3f}' for e in m['eps_bw_top5']]}")
        results['data']['tfim_critical'][str(n)] = {**m, 'n': n, 'E0': float(E0), 'time': float(dt)}
    except Exception as e:
        print(f"  n={n}: FAILED: {e}")
        import traceback; traceback.print_exc()
    with open('results/sprint_035b_bw_spectrum_scaling.json', 'w') as f:
        json.dump(results, f, indent=2)

# TFIM ordered
print("\n--- TFIM h/J=0.5 (ordered) ---")
results['data']['tfim_ordered'] = {}
for n in [8, 10, 12, 16]:
    elapsed = time.time() - t0
    if elapsed > 45:
        print(f"  n={n}: SKIP ({elapsed:.0f}s)")
        continue
    n_A = n // 2
    t1 = time.time()
    try:
        if n <= 10:
            rho_A, E0, basis = get_rho_A_exact_tfim(n, 0.5)
        else:
            rho_A, E0, basis = get_rho_A_dmrg_tfim(n, 0.5)
        m = spectrum_bw_accuracy(rho_A, n_A, 0.5, n, basis, 'tfim')
        dt = time.time() - t1
        print(f"  n={n}: R²_fitted={m['R2_fitted']:.6f}, R²_env={m['R2_envelope']:.6f}, "
              f"k={m['n_significant']}, time={dt:.1f}s")
        results['data']['tfim_ordered'][str(n)] = {**m, 'n': n, 'E0': float(E0), 'time': float(dt)}
    except Exception as e:
        print(f"  n={n}: FAILED: {e}")
    with open('results/sprint_035b_bw_spectrum_scaling.json', 'w') as f:
        json.dump(results, f, indent=2)

# XXZ at Δ=0.5 (XY phase, U(1))
print("\n--- XXZ Δ=0.5 (XY phase, U(1)) ---")
results['data']['xxz_xy'] = {}
for n in [8, 10, 12]:
    elapsed = time.time() - t0
    if elapsed > 50:
        print(f"  n={n}: SKIP ({elapsed:.0f}s)")
        continue
    n_A = n // 2
    t1 = time.time()
    try:
        if n <= 10:
            rho_A, E0, basis = get_rho_A_exact_xxz(n, 0.5)
        else:
            rho_A, E0, basis = get_rho_A_dmrg_xxz(n, 0.5)
        m = spectrum_bw_accuracy(rho_A, n_A, 0.5, n, basis, 'xxz')
        dt = time.time() - t1
        print(f"  n={n}: R²_fitted={m['R2_fitted']:.6f}, R²_env={m['R2_envelope']:.6f}, "
              f"k={m['n_significant']}, time={dt:.1f}s")
        results['data']['xxz_xy'][str(n)] = {**m, 'n': n, 'E0': float(E0), 'time': float(dt)}
    except Exception as e:
        print(f"  n={n}: FAILED: {e}")
    with open('results/sprint_035b_bw_spectrum_scaling.json', 'w') as f:
        json.dump(results, f, indent=2)

# Summary
print(f"\n{'='*60}")
print("SUMMARY: BW Spectrum R² vs System Size")
print(f"{'='*60}")
print(f"{'Model':>15} {'n':>4} {'R²_fit':>8} {'R²_env':>8} {'k':>4}")
for key in results['data']:
    for ns in sorted(results['data'][key].keys(), key=lambda x: int(x)):
        d = results['data'][key][ns]
        if 'R2_fitted' in d:
            print(f"{key:>15} {d['n']:4d} {d['R2_fitted']:8.4f} {d['R2_envelope']:8.4f} {d['n_significant']:4d}")

results['total_runtime'] = time.time() - t0
with open('results/sprint_035b_bw_spectrum_scaling.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nRuntime: {time.time()-t0:.1f}s")
