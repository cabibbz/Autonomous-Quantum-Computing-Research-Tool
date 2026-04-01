"""
Sprint 035a: TFIM BW Locality Size Scaling

BW accuracy vs system size for TFIM (Z₂, d=2).
n=8,10 exact diag, n=12,16,20 DMRG.
h/J=1.0 (critical) and h/J=0.5 (ordered).

Key insight: Pauli fraction is WRONG for scaling — denominator grows exponentially.
Correct metric: fidelity between ρ_A and BW thermal state ρ_BW = exp(-αH_BW)/Z.
Also track entanglement spectrum correlation.

TeNPy TFIChain: H = -J Σ Sigmax_i Sigmax_{i+1} - g Σ Sigmaz_i
Bond terms = XX, site terms = Z (opposite of my convention ZZ+X).
"""

import numpy as np
from scipy import linalg as la
import json, time, warnings
warnings.filterwarnings('ignore')

t0 = time.time()

I2 = np.eye(2)
X_std = np.array([[0,1],[1,0]])
Z_std = np.array([[1,0],[0,-1]])
# TeNPy: Sigmax = [[0,1],[1,0]], Sigmaz = diag(-1,1) (but same form as X)
Sx_tp = np.array([[0,1],[1,0]])
Sz_tp = np.diag([-1, 1])

def kron_list(ops):
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result

# ---- Exact diag ----
def tfim_hamiltonian(n, h_J):
    dim = 2**n
    H = np.zeros((dim, dim))
    for i in range(n-1):
        ops = [I2]*n; ops[i] = Z_std; ops[i+1] = Z_std
        H -= kron_list(ops)
    for i in range(n):
        ops = [I2]*n; ops[i] = X_std
        H -= h_J * kron_list(ops)
    return H

def partial_trace(rho, n_total, keep_qubits):
    n_keep = len(keep_qubits)
    trace_qubits = [q for q in range(n_total) if q not in keep_qubits]
    rho_tensor = rho.reshape([2]*n_total + [2]*n_total)
    for q in sorted(trace_qubits, reverse=True):
        nr = rho_tensor.ndim // 2
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q+nr)
    return rho_tensor.reshape(2**n_keep, 2**n_keep)

def get_rho_A_exact(n, h_J):
    H = tfim_hamiltonian(n, h_J)
    eigvals, eigvecs = la.eigh(H)
    psi = eigvecs[:, 0]
    rho = np.outer(psi, psi.conj())
    n_A = n // 2
    return partial_trace(rho, n, list(range(n_A))), eigvals[0], 'std'

# ---- DMRG ----
def get_rho_A_dmrg(n, h_J):
    import logging; logging.getLogger('tenpy').setLevel(logging.WARNING)
    from tenpy.models.tf_ising import TFIChain
    from tenpy.algorithms import dmrg
    from tenpy.networks.mps import MPS

    model = TFIChain({'L': n, 'J': 1.0, 'g': h_J, 'bc_MPS': 'finite'})
    psi = MPS.from_product_state(model.lat.mps_sites(), ['up'] * n, bc='finite')
    eng = dmrg.TwoSiteDMRGEngine(psi, model, {
        'mixer': True, 'max_E_err': 1e-12, 'max_S_err': 1e-10,
        'trunc_params': {'chi_max': 200, 'svd_min': 1e-14}, 'max_sweeps': 30,
    })
    E0, _ = eng.run()
    n_A = n // 2
    rho_tensor = psi.get_rho_segment(list(range(n_A)))
    rho_np = rho_tensor.to_ndarray().reshape(2**n_A, 2**n_A)
    return rho_np, E0, 'tenpy'

# ---- BW Hamiltonian ----
def build_H_BW(n_A, h_J, n_total, basis='std'):
    """Build BW Hamiltonian with sin_inv envelope."""
    if basis == 'tenpy':
        bond_pair = (Sx_tp, Sx_tp)
        site_op = Sz_tp
    else:
        bond_pair = (Z_std, Z_std)
        site_op = X_std

    dim = 2**n_A
    L = n_total
    H_BW = np.zeros((dim, dim))

    for i in range(n_A - 1):
        dist = n_A - 1 - i
        beta = np.sin(np.pi * dist / L)
        ops = [I2]*n_A; ops[i] = bond_pair[0]; ops[i+1] = bond_pair[1]
        H_BW -= beta * kron_list(ops)
    for i in range(n_A):
        dist = n_A - 0.5 - i
        beta = np.sin(np.pi * dist / L)
        ops = [I2]*n_A; ops[i] = site_op
        H_BW -= h_J * beta * kron_list(ops)

    return H_BW

# ---- Metrics ----
def entanglement_hamiltonian(rho_A):
    eigvals, eigvecs = la.eigh(rho_A)
    eigvals = np.maximum(eigvals, 1e-15)
    return -eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.T

def fidelity(rho, sigma):
    """Quantum fidelity F(ρ,σ) = (Tr √(√ρ σ √ρ))²."""
    sqrt_rho = la.sqrtm(rho)
    M = sqrt_rho @ sigma @ sqrt_rho
    # Eigenvalues of M should be non-negative
    evals = la.eigvalsh(M)
    evals = np.maximum(evals, 0)
    return (np.sum(np.sqrt(evals)))**2

def bw_metrics(rho_A, H_E, n_A, h_J, n_total, basis):
    """Compute BW accuracy via multiple metrics."""
    dim = 2**n_A
    H_BW = build_H_BW(n_A, h_J, n_total, basis)

    # Remove trace from both
    H_E_tl = H_E - np.trace(H_E)/dim * np.eye(dim)
    H_BW_tl = H_BW - np.trace(H_BW)/dim * np.eye(dim)
    norm_E = la.norm(H_E_tl, 'fro')
    norm_BW = la.norm(H_BW_tl, 'fro')

    # 1) Matrix overlap (Frobenius cosine)
    overlap = np.real(np.trace(H_E_tl @ H_BW_tl)) / (norm_E * norm_BW) if norm_E > 1e-15 and norm_BW > 1e-15 else 0.0

    # 2) Optimal rescaling + variance captured
    alpha_opt = np.real(np.trace(H_E_tl @ H_BW_tl)) / np.real(np.trace(H_BW_tl @ H_BW_tl)) if norm_BW > 1e-15 else 0.0
    residual = H_E_tl - alpha_opt * H_BW_tl
    var_captured = 1 - (la.norm(residual, 'fro') / norm_E)**2 if norm_E > 1e-15 else 0.0

    # 3) Density matrix fidelity: compare ρ_A with ρ_BW = exp(-α H_BW)/Z
    rho_BW = la.expm(-alpha_opt * H_BW)
    rho_BW = rho_BW / np.trace(rho_BW)
    F = fidelity(rho_A, rho_BW)

    # 4) Entanglement spectrum correlation
    # Compare sorted eigenvalues of ρ_A and ρ_BW
    evals_A = sorted(la.eigvalsh(rho_A), reverse=True)
    evals_BW = sorted(la.eigvalsh(rho_BW), reverse=True)
    # Keep significant eigenvalues
    k = max(2, np.sum(np.array(evals_A) > 1e-8))
    evals_A_k = evals_A[:k]
    evals_BW_k = evals_BW[:k]
    if len(evals_A_k) > 1:
        spec_corr = np.corrcoef(evals_A_k, evals_BW_k)[0, 1]
    else:
        spec_corr = 1.0

    return {
        'overlap': float(overlap),
        'var_captured': float(var_captured),
        'alpha_opt': float(alpha_opt),
        'fidelity': float(np.real(F)),
        'spectrum_corr': float(spec_corr),
        'n_significant_evals': int(k),
    }

def tfim_pauli_fraction(H_E, n_A, basis):
    """TFIM-type Pauli fraction (for comparison with Sprint 032)."""
    if basis == 'tenpy':
        bond_pair, site_op = (Sx_tp, Sx_tp), Sz_tp
    else:
        bond_pair, site_op = (Z_std, Z_std), X_std

    dim = 2**n_A
    H_tl = H_E - np.trace(H_E)/dim * np.eye(dim)
    total = la.norm(H_tl, 'fro')**2
    if total < 1e-15:
        return 0.0

    tfim_sq = 0.0
    for site in range(n_A):
        ops = [I2]*n_A; ops[site] = site_op
        c = np.real(np.trace(H_tl @ kron_list(ops))) / dim
        tfim_sq += c**2 * dim
    for site in range(n_A - 1):
        ops = [I2]*n_A; ops[site] = bond_pair[0]; ops[site+1] = bond_pair[1]
        c = np.real(np.trace(H_tl @ kron_list(ops))) / dim
        tfim_sq += c**2 * dim

    return tfim_sq / total

# ---- Main ----
system_sizes = [8, 10, 12, 16, 20]
h_J_values = [1.0, 0.5]
h_J_labels = ['critical', 'ordered']

results = {'experiment': '035a', 'description': 'TFIM BW locality size scaling', 'data': {}}

print("=== Sprint 035a: TFIM BW Locality Size Scaling ===")
print("Metrics: fidelity F(ρ_A, ρ_BW), overlap, var_captured, spectrum_corr")

for h_J, label in zip(h_J_values, h_J_labels):
    print(f"\n{'='*60}")
    print(f"h/J = {h_J} ({label})")
    print(f"{'='*60}")
    results['data'][label] = {}

    for n in system_sizes:
        n_A = n // 2
        elapsed = time.time() - t0
        if elapsed > 48:
            print(f"  n={n}: SKIPPED ({elapsed:.0f}s)")
            continue

        print(f"\n  n={n}, n_A={n_A}:")
        t1 = time.time()

        try:
            if n <= 10:
                rho_A, E0, basis = get_rho_A_exact(n, h_J)
            else:
                rho_A, E0, basis = get_rho_A_dmrg(n, h_J)

            evals = la.eigvalsh(rho_A)
            evals_pos = evals[evals > 1e-15]
            S = -np.sum(evals_pos * np.log2(evals_pos))
            H_E = entanglement_hamiltonian(rho_A)

            # BW metrics
            m = bw_metrics(rho_A, H_E, n_A, h_J, n, basis)

            # Pauli fraction (for comparison, only feasible for n_A <= 10)
            pauli_frac = tfim_pauli_fraction(H_E, n_A, basis) if n_A <= 10 else None

            dt = time.time() - t1
            print(f"    E0={E0:.6f}, S={S:.4f} bits")
            print(f"    Fidelity F(ρ,ρ_BW) = {m['fidelity']:.6f}")
            print(f"    Overlap = {m['overlap']:.4f}, Var captured = {m['var_captured']*100:.1f}%")
            print(f"    Spectrum corr = {m['spectrum_corr']:.6f}")
            if pauli_frac is not None:
                print(f"    Pauli fraction = {pauli_frac*100:.1f}% (Sprint 032 metric)")
            print(f"    Time: {dt:.1f}s")

            results['data'][label][str(n)] = {
                'n': n, 'n_A': n_A, 'method': 'exact' if n<=10 else 'DMRG',
                'E0': float(E0), 'entropy_bits': float(S),
                'fidelity': m['fidelity'],
                'overlap': m['overlap'],
                'var_captured': m['var_captured'],
                'alpha_opt': m['alpha_opt'],
                'spectrum_corr': m['spectrum_corr'],
                'n_significant_evals': m['n_significant_evals'],
                'pauli_fraction': float(pauli_frac) if pauli_frac is not None else None,
                'time_seconds': float(dt),
            }
        except Exception as e:
            print(f"    FAILED: {e}")
            import traceback; traceback.print_exc()
            results['data'][label][str(n)] = {'n': n, 'error': str(e)}

        with open('results/sprint_035a_tfim_bw_scaling.json', 'w') as f:
            json.dump(results, f, indent=2)

# Summary
print(f"\n{'='*60}")
print("SUMMARY")
print(f"{'='*60}")
print(f"{'n':>4} {'Phase':>10} {'Fidelity':>9} {'Overlap':>8} {'Var%':>6} {'SpecCorr':>9} {'Pauli%':>7}")
for label in h_J_labels:
    if label not in results['data']:
        continue
    for ns in sorted(results['data'][label].keys(), key=lambda x: int(x)):
        d = results['data'][label][ns]
        if 'error' in d:
            continue
        pf = f"{d['pauli_fraction']*100:.1f}" if d.get('pauli_fraction') is not None else '  -'
        print(f"{d['n']:4d} {label:>10} {d['fidelity']:9.6f} {d['overlap']:8.4f} {d['var_captured']*100:5.1f}% {d['spectrum_corr']:9.6f} {pf:>7}")

results['total_runtime'] = time.time() - t0
with open('results/sprint_035a_tfim_bw_scaling.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nRuntime: {time.time()-t0:.1f}s")
