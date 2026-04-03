"""
Sprint 032a: TFIM Entanglement Hamiltonian at Criticality — Bisognano-Wichmann Test

Compute H_E = -log(ρ_A) for the TFIM ground state at criticality (h/J=1.0).
Decompose H_E into Pauli basis to identify dominant terms.
Compare to BW prediction: H_BW = Σ β(i) h_i where h_i are local TFIM terms in subsystem A
and β(i) follows the CFT envelope.

System: n=8 qubits, subsystem A = left 4 qubits.
"""

import numpy as np
from scipy import linalg as la
import json, time

t0 = time.time()

# ---- Build TFIM Hamiltonian ----
def pauli_matrices():
    I = np.eye(2)
    X = np.array([[0,1],[1,0]])
    Z = np.array([[1,0],[0,-1]])
    return I, X, Z

def kron_list(ops):
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result

def tfim_hamiltonian(n, h_over_J):
    """H = -J Σ Z_i Z_{i+1} - h Σ X_i (open boundary)"""
    I, X, Z = pauli_matrices()
    dim = 2**n
    H = np.zeros((dim, dim))
    # ZZ terms
    for i in range(n-1):
        ops = [I]*n
        ops[i] = Z
        ops[i+1] = Z
        H -= kron_list(ops)
    # X terms
    for i in range(n):
        ops = [I]*n
        ops[i] = X
        H -= h_over_J * kron_list(ops)
    return H

def partial_trace(rho, n_total, keep_qubits):
    """Trace out qubits NOT in keep_qubits."""
    n_keep = len(keep_qubits)
    n_trace = n_total - n_keep
    trace_qubits = [q for q in range(n_total) if q not in keep_qubits]

    # Reshape into tensor
    rho_tensor = rho.reshape([2]*n_total + [2]*n_total)

    # Trace out qubits one at a time (from highest index to avoid reindexing)
    for q in sorted(trace_qubits, reverse=True):
        # Contract axis q with axis q+remaining_total
        n_remaining = rho_tensor.ndim // 2
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q+n_remaining)

    dim_keep = 2**n_keep
    return rho_tensor.reshape(dim_keep, dim_keep)

def entanglement_hamiltonian(rho_A):
    """H_E = -log(ρ_A), regularized for zero eigenvalues."""
    eigvals, eigvecs = la.eigh(rho_A)
    # Regularize: replace zeros/negatives with tiny value
    eigvals = np.maximum(eigvals, 1e-15)
    H_E = -eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.T
    return H_E

# ---- Pauli decomposition of H_E ----
def pauli_basis_4q():
    """Generate all 256 Pauli operators for 4 qubits with labels."""
    I = np.eye(2)
    X = np.array([[0,1],[1,0]])
    Y = np.array([[0,-1j],[1j,0]])
    Z = np.array([[1,0],[0,-1]])
    paulis = {'I': I, 'X': X, 'Y': Y, 'Z': Z}

    basis = {}
    for a in 'IXYZ':
        for b in 'IXYZ':
            for c in 'IXYZ':
                for d in 'IXYZ':
                    label = a+b+c+d
                    op = np.kron(np.kron(np.kron(paulis[a], paulis[b]), paulis[c]), paulis[d])
                    basis[label] = op
    return basis

def decompose_in_pauli(H, basis):
    """Decompose Hermitian matrix H = Σ c_i P_i, c_i = Tr(H P_i) / 2^n."""
    n = int(np.log2(H.shape[0]))
    coeffs = {}
    for label, P in basis.items():
        c = np.real(np.trace(H @ P)) / (2**n)
        if abs(c) > 1e-10:
            coeffs[label] = c
    return coeffs

# ---- BW prediction ----
def bw_hamiltonian_tfim(n_A, h_over_J):
    """
    BW prediction for TFIM subsystem: H_BW = Σ β(i) h_i
    where h_i are local TFIM terms and β(i) is the entanglement temperature envelope.

    For open boundary 1D CFT: β(i) ≈ (π/L) * sqrt(i * (L-i)) for bulk,
    but for half-chain of open-boundary system, the standard form is:
    β(i) = (π v / L_total) * sin(π * x_i / L_total) where x_i is position in full chain.

    We use the lattice BW form: β_bond(i) for bond between sites i and i+1 in subsystem A,
    and β_site(i) for site i. For half-chain cut of open-boundary chain of length L:
    β_bond(i) = (π/L) * (i+1) * (L - i - 1)  ... approximate parabolic
    β_site(i) = average of neighboring bond weights.

    Simplified: use the continuum prediction β(x) = π * x * (1 - x/L_A) / L_A
    evaluated at lattice sites, where x = distance from the entanglement cut.
    """
    I2, X, Z = pauli_matrices()
    dim = 2**n_A

    # Position from the entanglement cut (right boundary of A)
    # Site 0 is furthest from cut, site n_A-1 is adjacent to cut
    # BW envelope: sites far from cut have larger β (cooler), sites near cut have smaller β (hotter)
    # Standard: β(i) ∝ distance from cut for half-infinite, or parabolic for finite

    # For finite chain with open boundaries, half-cut:
    # The BW envelope for site i (0-indexed, in left half of chain of length L):
    # β(i) = (π/L) * sin(π*(i+0.5)/L) where L is total chain length
    L_total = 2 * n_A  # total chain length (we're taking left half)

    H_BW = np.zeros((dim, dim))

    # ZZ terms within subsystem A
    for i in range(n_A - 1):
        beta_bond = np.pi / L_total * np.sin(np.pi * (i + 1) / L_total)
        ops = [I2]*n_A
        ops[i] = Z
        ops[i+1] = Z
        H_BW -= beta_bond * kron_list(ops)

    # X terms within subsystem A
    for i in range(n_A):
        beta_site = np.pi / L_total * np.sin(np.pi * (i + 0.5) / L_total)
        ops = [I2]*n_A
        ops[i] = X
        H_BW -= h_over_J * beta_site * kron_list(ops)

    return H_BW

def matrix_overlap(A, B):
    """Normalized overlap: Tr(A†B) / (||A|| ||B||)"""
    nA = la.norm(A, 'fro')
    nB = la.norm(B, 'fro')
    if nA < 1e-15 or nB < 1e-15:
        return 0.0
    return np.real(np.trace(A.conj().T @ B)) / (nA * nB)

# ---- Main computation ----
n = 8
n_A = 4
h_J_crit = 1.0  # TFIM critical point

print(f"=== Sprint 032a: TFIM Entanglement Hamiltonian (n={n}, n_A={n_A}) ===")
print(f"h/J = {h_J_crit} (critical point)")

# Ground state
H = tfim_hamiltonian(n, h_J_crit)
eigvals, eigvecs = la.eigh(H)
psi_gs = eigvecs[:, 0]
print(f"Ground state energy: {eigvals[0]:.6f}")
print(f"Energy gap: {eigvals[1] - eigvals[0]:.6f}")

# Reduced density matrix
rho = np.outer(psi_gs, psi_gs.conj())
rho_A = partial_trace(rho, n, list(range(n_A)))
print(f"ρ_A shape: {rho_A.shape}")
print(f"Tr(ρ_A) = {np.trace(rho_A).real:.6f}")

# Entanglement entropy
eigvals_rho = la.eigvalsh(rho_A)
eigvals_rho = eigvals_rho[eigvals_rho > 1e-15]
S = -np.sum(eigvals_rho * np.log2(eigvals_rho))
print(f"Entanglement entropy S = {S:.4f} bits")

# Entanglement Hamiltonian
H_E = entanglement_hamiltonian(rho_A)
print(f"\nH_E computed. Frobenius norm: {la.norm(H_E, 'fro'):.4f}")

# Pauli decomposition
print("\n--- Pauli Decomposition of H_E ---")
basis = pauli_basis_4q()
coeffs = decompose_in_pauli(H_E, basis)

# Sort by magnitude
sorted_coeffs = sorted(coeffs.items(), key=lambda x: abs(x[1]), reverse=True)
print(f"Non-zero Pauli terms: {len(sorted_coeffs)}")
print("\nTop 20 terms:")
for label, c in sorted_coeffs[:20]:
    weight = sum(1 for ch in label if ch != 'I')
    print(f"  {label} (weight {weight}): {c:+.6f}")

# Categorize by Pauli weight
weight_norms = {}
for label, c in coeffs.items():
    w = sum(1 for ch in label if ch != 'I')
    weight_norms[w] = weight_norms.get(w, 0) + c**2
total_norm_sq = sum(c**2 for c in coeffs.values())
print("\nNorm² by Pauli weight:")
for w in sorted(weight_norms.keys()):
    frac = weight_norms[w] / total_norm_sq
    print(f"  Weight {w}: {weight_norms[w]:.4f} ({frac*100:.1f}%)")

# Identify TFIM-like terms (ZZ and X)
tfim_terms = {}
non_tfim_terms = {}
for label, c in coeffs.items():
    # Check if this is a ZZ term (exactly two Z's on adjacent sites) or X term (single X)
    is_tfim = False
    if label.count('I') == 3 and label.count('X') == 1:  # single X
        is_tfim = True
    elif label.count('I') == 2 and label.count('Z') == 2:  # check adjacency
        zz_positions = [i for i, ch in enumerate(label) if ch == 'Z']
        if len(zz_positions) == 2 and abs(zz_positions[0] - zz_positions[1]) == 1:
            is_tfim = True
    if label == 'IIII':
        continue  # skip identity
    if is_tfim:
        tfim_terms[label] = c
    else:
        non_tfim_terms[label] = c

tfim_norm_sq = sum(c**2 for c in tfim_terms.values())
non_tfim_norm_sq = sum(c**2 for c in non_tfim_terms.values())
print(f"\nTFIM-type terms: {len(tfim_terms)}, norm² = {tfim_norm_sq:.4f} ({tfim_norm_sq/total_norm_sq*100:.1f}%)")
print(f"Non-TFIM terms: {len(non_tfim_terms)}, norm² = {non_tfim_norm_sq:.4f} ({non_tfim_norm_sq/total_norm_sq*100:.1f}%)")

print("\nTFIM-type terms:")
for label, c in sorted(tfim_terms.items()):
    print(f"  {label}: {c:+.6f}")

print("\nLargest non-TFIM terms:")
sorted_non = sorted(non_tfim_terms.items(), key=lambda x: abs(x[1]), reverse=True)
for label, c in sorted_non[:10]:
    print(f"  {label}: {c:+.6f}")

# ---- BW comparison ----
print("\n--- Bisognano-Wichmann Comparison ---")
H_BW = bw_hamiltonian_tfim(n_A, h_J_crit)

# Remove identity/trace component from both for fair comparison
H_E_traceless = H_E - np.trace(H_E) / H_E.shape[0] * np.eye(H_E.shape[0])
H_BW_traceless = H_BW - np.trace(H_BW) / H_BW.shape[0] * np.eye(H_BW.shape[0])

# Overall BW fidelity
bw_fidelity = matrix_overlap(H_E_traceless, H_BW_traceless)
print(f"BW fidelity (normalized overlap): {bw_fidelity:.6f}")

# Optimal rescaling: find α that minimizes ||H_E - α H_BW||²
# α = Tr(H_E H_BW) / Tr(H_BW H_BW)
alpha = np.real(np.trace(H_E_traceless @ H_BW_traceless)) / np.real(np.trace(H_BW_traceless @ H_BW_traceless))
residual = H_E_traceless - alpha * H_BW_traceless
residual_frac = la.norm(residual, 'fro') / la.norm(H_E_traceless, 'fro')
print(f"Optimal rescaling α = {alpha:.4f}")
print(f"Residual ||H_E - α·H_BW|| / ||H_E|| = {residual_frac:.6f}")
print(f"BW captures {(1 - residual_frac**2)*100:.1f}% of H_E variance")

# Extract effective β(i) from H_E coefficients
print("\n--- Effective Entanglement Temperature ---")
print("BW prediction: β(i) ∝ sin(πi/L)")
I2, X, Z = pauli_matrices()
L_total = 2 * n_A

# ZZ bond weights
print("\nZZ bond coefficients (negative = ferromagnetic coupling in H_E):")
zz_labels = ['ZZII', 'IZZI', 'IIZZ']
bw_zz = []
he_zz = []
for idx, label in enumerate(zz_labels):
    i = idx  # bond between site i and i+1
    beta_bw = np.pi / L_total * np.sin(np.pi * (i + 1) / L_total)
    c_he = coeffs.get(label, 0.0)
    bw_zz.append(beta_bw)
    he_zz.append(-c_he)  # negate because H = -J ZZ, so H_E coeff should be negative
    print(f"  Bond {i}-{i+1} ({label}): H_E coeff = {c_he:+.6f}, BW β = {beta_bw:.6f}")

# X site weights
print("\nX site coefficients:")
x_labels = ['XIII', 'IXII', 'IIXI', 'IIIX']
bw_x = []
he_x = []
for idx, label in enumerate(x_labels):
    i = idx
    beta_bw = np.pi / L_total * np.sin(np.pi * (i + 0.5) / L_total)
    c_he = coeffs.get(label, 0.0)
    bw_x.append(beta_bw)
    he_x.append(-c_he)
    print(f"  Site {i} ({label}): H_E coeff = {c_he:+.6f}, BW β = {beta_bw:.6f}")

# Fit: are the ratios β_HE / β_BW constant?
print("\nRatio test (should be constant if BW holds):")
all_ratios = []
for idx, label in enumerate(zz_labels):
    if bw_zz[idx] > 1e-10:
        ratio = he_zz[idx] / bw_zz[idx]
        all_ratios.append(ratio)
        print(f"  ZZ bond {idx}-{idx+1}: ratio = {ratio:.4f}")
for idx, label in enumerate(x_labels):
    if bw_x[idx] > 1e-10:
        ratio = he_x[idx] / bw_x[idx]
        all_ratios.append(ratio)
        print(f"  X site {idx}: ratio = {ratio:.4f}")

if all_ratios:
    mean_ratio = np.mean(all_ratios)
    std_ratio = np.std(all_ratios)
    cv_ratio = std_ratio / abs(mean_ratio) if abs(mean_ratio) > 1e-10 else float('inf')
    print(f"\n  Mean ratio: {mean_ratio:.4f}")
    print(f"  Std ratio: {std_ratio:.4f}")
    print(f"  CV (should be ~0 if BW perfect): {cv_ratio:.4f}")

# ---- Decompose residual ----
print("\n--- Residual Analysis (what BW misses) ---")
residual_coeffs = decompose_in_pauli(residual, basis)
sorted_res = sorted(residual_coeffs.items(), key=lambda x: abs(x[1]), reverse=True)
print("Top 10 residual Pauli terms:")
for label, c in sorted_res[:10]:
    weight = sum(1 for ch in label if ch != 'I')
    print(f"  {label} (weight {weight}): {c:+.6f}")

# Weight distribution of residual
res_weight_norms = {}
res_total = sum(c**2 for c in residual_coeffs.values())
for label, c in residual_coeffs.items():
    w = sum(1 for ch in label if ch != 'I')
    res_weight_norms[w] = res_weight_norms.get(w, 0) + c**2
print("\nResidual norm² by weight:")
for w in sorted(res_weight_norms.keys()):
    frac = res_weight_norms[w] / res_total if res_total > 0 else 0
    print(f"  Weight {w}: {res_weight_norms[w]:.6f} ({frac*100:.1f}%)")

# ---- Save results ----
results = {
    'experiment': '032a',
    'description': 'TFIM entanglement Hamiltonian at criticality - BW test',
    'parameters': {'n': n, 'n_A': n_A, 'h_over_J': h_J_crit},
    'ground_state_energy': float(eigvals[0]),
    'energy_gap': float(eigvals[1] - eigvals[0]),
    'entanglement_entropy_bits': float(S),
    'H_E_frobenius_norm': float(la.norm(H_E, 'fro')),
    'n_nonzero_pauli_terms': len(sorted_coeffs),
    'top_20_pauli_terms': {label: float(c) for label, c in sorted_coeffs[:20]},
    'weight_distribution': {str(w): float(weight_norms[w]) for w in sorted(weight_norms.keys())},
    'tfim_terms_fraction': float(tfim_norm_sq / total_norm_sq),
    'bw_fidelity': float(bw_fidelity),
    'bw_optimal_rescaling': float(alpha),
    'bw_residual_fraction': float(residual_frac),
    'bw_variance_captured': float(1 - residual_frac**2),
    'zz_bond_coeffs': {f'bond_{i}_{i+1}': float(he_zz[i]) for i in range(len(zz_labels))},
    'x_site_coeffs': {f'site_{i}': float(he_x[i]) for i in range(len(x_labels))},
    'bw_zz_predictions': {f'bond_{i}_{i+1}': float(bw_zz[i]) for i in range(len(zz_labels))},
    'bw_x_predictions': {f'site_{i}': float(bw_x[i]) for i in range(len(x_labels))},
    'ratio_mean': float(mean_ratio) if all_ratios else None,
    'ratio_cv': float(cv_ratio) if all_ratios else None,
    'top_10_residual_terms': {label: float(c) for label, c in sorted_res[:10]},
    'residual_weight_distribution': {str(w): float(res_weight_norms.get(w, 0)) for w in range(5)},
    'runtime_seconds': time.time() - t0
}

with open('results/sprint_032a_tfim_entanglement_hamiltonian.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved. Runtime: {time.time()-t0:.1f}s")
