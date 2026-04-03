"""
Sprint 034a: G-Invariant Operator Counting

For a subsystem of n_A sites with local dimension d, count:
1. Total operator space dimension: d^{2*n_A} (Hermitian operators on subsystem)
2. Hamiltonian-type operator dimension: nearest-neighbor bonds + single-site terms
3. G-invariant operator dimension: operators commuting with all group elements

The ratio dim(H_terms)/dim(G_inv) should predict BW locality.

Groups tested:
- Z₂ (d=2): TFIM, measured BW locality 91%
- U(1) (d=2): XXZ, measured BW locality 100%
- S₃ (d=3): Potts, measured BW locality 76.5%
- Z₃ (d=3): Clock model, predicted < 76.5%

Method: Represent group action on operator space as superoperator,
project onto invariant subspace, count dimensions.
"""

import numpy as np
from scipy import linalg as la
import json, time
from itertools import product

t0 = time.time()

def group_action_on_operator(U, O):
    """Apply symmetry U to operator O: O → U O U†"""
    return U @ O @ U.conj().T

def projector_onto_invariants(group_reps, dim):
    """
    For a group with representations {U_g}, compute the projector onto
    G-invariant operators: P = (1/|G|) Σ_g (U_g ⊗ U_g*)

    An operator O (as a dim×dim matrix, flattened to dim² vector) is
    G-invariant iff U_g O U_g† = O for all g.

    The superoperator acting on vec(O) is: U_g ⊗ conj(U_g)
    The projector onto the invariant subspace is: P = (1/|G|) Σ_g U_g ⊗ conj(U_g)
    """
    dim_sq = dim * dim
    P = np.zeros((dim_sq, dim_sq), dtype=complex)
    for U in group_reps:
        # Superoperator: (U ⊗ U*) acting on vec(O) gives vec(U O U†)
        P += np.kron(U, U.conj())
    P /= len(group_reps)
    return P

def count_invariant_hermitian(group_reps, dim):
    """
    Count dimension of G-invariant Hermitian operators (traceless).

    Returns: dim of invariant subspace of traceless Hermitian operators.
    """
    P = projector_onto_invariants(group_reps, dim)

    # Rank of P = dimension of invariant subspace (all operators, not just Hermitian)
    eigvals = la.eigvalsh(P @ P.conj().T)
    rank_all = np.sum(eigvals > 0.5)  # eigenvalues are 0 or 1

    # But we want Hermitian operators. Since P commutes with the Hermitian constraint
    # (if O is Hermitian and G-invariant, so is U O U†), we can work in Hermitian subspace.
    # For d×d matrices: dim(Hermitian) = d², dim(traceless Hermitian) = d²-1
    # The G-invariant Hermitian subspace has dimension = rank of P restricted to Hermitian subspace

    # Build a real basis for Hermitian matrices (d×d)
    herm_basis = []
    # Diagonal elements
    for i in range(dim):
        M = np.zeros((dim, dim), dtype=complex)
        M[i, i] = 1.0
        herm_basis.append(M)
    # Off-diagonal elements (real and imaginary parts)
    for i in range(dim):
        for j in range(i+1, dim):
            M_re = np.zeros((dim, dim), dtype=complex)
            M_re[i, j] = 1.0
            M_re[j, i] = 1.0
            herm_basis.append(M_re)
            M_im = np.zeros((dim, dim), dtype=complex)
            M_im[i, j] = 1j
            M_im[j, i] = -1j
            herm_basis.append(M_im)

    # Apply projector to each basis element and check linear independence
    projected = []
    for M in herm_basis:
        v = M.ravel()
        Pv = P @ v
        # Pv should be close to a Hermitian matrix
        projected.append(np.real(Pv))  # take real part of vec

    # Stack and find rank
    mat = np.array(projected)  # (d², d²) matrix
    _, s, _ = la.svd(mat)
    rank_herm = np.sum(s > 1e-10)

    # Subtract 1 for traceless constraint (identity is always invariant)
    rank_traceless = rank_herm - 1

    return rank_traceless, rank_herm, rank_all

def count_hamiltonian_terms(n_A, d, model_type):
    """
    Count the number of independent Hamiltonian-type operator degrees of freedom
    on a subsystem of n_A sites with local dimension d.

    For nearest-neighbor models:
    - Bond operators: (n_A - 1) bonds, each with some number of independent terms
    - Site operators: n_A sites, each with some number of independent terms

    Returns: number of independent real parameters in the Hamiltonian-type subspace
    """
    if model_type == 'TFIM':
        # ZZ bonds: 1 parameter per bond (n_A-1 bonds)
        # X sites: 1 parameter per site (n_A sites)
        return (n_A - 1) + n_A
    elif model_type == 'XXZ':
        # XX+YY+ZZ bonds: 2 parameters per bond (XX+YY coupled, ZZ separate)
        # Actually for the H = -J(XX+YY+Δ ZZ) there are 2 independent terms per bond
        # But in BW we fit each independently, so:
        # XX, YY, ZZ per bond: but XX and YY are related by U(1), so 2 independent
        # Site terms: none for XXZ (no field)
        # Actually let's count: each bond has XX+YY (1 term) and ZZ (1 term) = 2 per bond
        return 2 * (n_A - 1)
    elif model_type == 'Potts':
        # δ(σ_i, σ_{i+1}) bonds: 1 parameter per bond (n_A-1 bonds)
        # (P + P†) sites: 1 parameter per site (n_A sites)
        return (n_A - 1) + n_A
    elif model_type == 'Clock':
        # Z_i Z_{i+1}† + h.c. bonds: 1 parameter per bond
        # X + X† sites: 1 parameter per site
        return (n_A - 1) + n_A
    return 0

def multi_site_group_reps(single_site_reps, n_sites):
    """Build tensor product representations for n sites."""
    multi_reps = []
    for U in single_site_reps:
        U_multi = U
        for _ in range(n_sites - 1):
            U_multi = np.kron(U_multi, U)
        multi_reps.append(U_multi)
    return multi_reps

# ---- Define symmetry groups ----

# Z₂ on d=2 (TFIM): X generates Z₂
X2 = np.array([[0, 1], [1, 0]], dtype=complex)
I2 = np.eye(2, dtype=complex)
Z2_reps = [I2, X2]

# U(1) on d=2 (XXZ): Rz(θ) for many θ values (approximate continuous group)
# For counting, use enough angles to capture the full invariant subspace
n_angles = 20
U1_reps = []
for k in range(n_angles):
    theta = 2 * np.pi * k / n_angles
    Rz = np.array([[np.exp(1j*theta/2), 0], [0, np.exp(-1j*theta/2)]], dtype=complex)
    U1_reps.append(Rz)

# S₃ on d=3 (Potts): All permutations of {0,1,2}
from itertools import permutations
S3_reps = []
for perm in permutations(range(3)):
    U = np.zeros((3, 3), dtype=complex)
    for i, j in enumerate(perm):
        U[j, i] = 1.0
    S3_reps.append(U)

# Z₃ on d=3 (Clock): Cyclic permutation X|k⟩ = |k+1 mod 3⟩
I3 = np.eye(3, dtype=complex)
omega = np.exp(2j * np.pi / 3)
X3 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=complex)  # cyclic perm
Z3_reps = [I3, X3, X3 @ X3]  # {I, X, X²}

print("=== Sprint 034a: G-Invariant Operator Counting ===\n")

# ---- Single-site analysis ----
print("--- Single-site (1-body) operator counting ---")
print(f"{'Group':<8} {'d':<4} {'|G|':<5} {'dim(Herm_tl)':<14} {'dim(G-inv_tl)':<15} {'ratio':<8}")

single_site_results = {}

for name, d, reps, order in [
    ('Z₂', 2, Z2_reps, 2),
    ('U(1)', 2, U1_reps, '∞'),
    ('Z₃', 3, Z3_reps, 3),
    ('S₃', 3, S3_reps, 6),
]:
    n_tl, n_h, n_all = count_invariant_hermitian(reps, d)
    total_tl = d*d - 1
    ratio = n_tl / total_tl if total_tl > 0 else 0
    print(f"{name:<8} {d:<4} {str(order):<5} {total_tl:<14} {n_tl:<15} {ratio:<8.4f}")
    single_site_results[name] = {
        'd': d, 'group_order': str(order),
        'total_traceless_hermitian': total_tl,
        'invariant_traceless_hermitian': n_tl,
        'ratio': ratio
    }

# ---- Two-site analysis (2-body operators) ----
print("\n--- Two-site (2-body) operator counting ---")
print(f"{'Group':<8} {'d':<4} {'|G|':<5} {'dim(Herm_tl)':<14} {'dim(G-inv_tl)':<15} {'ratio':<8}")

two_site_results = {}

for name, d, single_reps, order in [
    ('Z₂', 2, Z2_reps, 2),
    ('U(1)', 2, U1_reps, '∞'),
    ('Z₃', 3, Z3_reps, 3),
    ('S₃', 3, S3_reps, 6),
]:
    multi_reps = multi_site_group_reps(single_reps, 2)
    dim2 = d * d
    n_tl, n_h, n_all = count_invariant_hermitian(multi_reps, dim2)
    total_tl = dim2*dim2 - 1
    ratio = n_tl / total_tl if total_tl > 0 else 0
    print(f"{name:<8} {d:<4} {str(order):<5} {total_tl:<14} {n_tl:<15} {ratio:<8.4f}")
    two_site_results[name] = {
        'd': d, 'group_order': str(order),
        'total_traceless_hermitian': total_tl,
        'invariant_traceless_hermitian': n_tl,
        'ratio': ratio
    }

# ---- n_A=4 site analysis (full subsystem as used in BW experiments) ----
print("\n--- 4-site subsystem operator counting (matching BW experiments) ---")
print(f"{'Group':<8} {'d':<4} {'dim_A':<7} {'dim(Herm_tl)':<14} {'dim(G-inv_tl)':<15} {'ratio':<8} {'H_terms':<8}")

four_site_results = {}

for name, d, single_reps, order, model in [
    ('Z₂', 2, Z2_reps, 2, 'TFIM'),
    ('U(1)', 2, U1_reps, '∞', 'XXZ'),
]:
    # d=2 models: 4 sites → dim_A=16, 256 operators
    multi_reps = multi_site_group_reps(single_reps, 4)
    dim_A = d**4
    n_tl, n_h, n_all = count_invariant_hermitian(multi_reps, dim_A)
    total_tl = dim_A*dim_A - 1
    ratio = n_tl / total_tl if total_tl > 0 else 0
    h_terms = count_hamiltonian_terms(4, d, model)
    h_ratio = h_terms / n_tl if n_tl > 0 else 0
    print(f"{name:<8} {d:<4} {dim_A:<7} {total_tl:<14} {n_tl:<15} {ratio:<8.4f} {h_terms:<8} (H/G-inv={h_ratio:.4f})")
    four_site_results[name] = {
        'd': d, 'dim_A': dim_A, 'group_order': str(order), 'model': model,
        'total_traceless_hermitian': total_tl,
        'invariant_traceless_hermitian': n_tl,
        'ratio_inv_to_total': ratio,
        'hamiltonian_terms': h_terms,
        'ratio_H_to_inv': h_ratio
    }

# d=3 models: 4 sites → dim_A=81, 6561 operators — too large for full SVD
# Use 2-site and extrapolate, or compute for n_A=2,3 only
print("\n--- d=3 models: computing for n_A=2 (tractable) ---")
for name, d, single_reps, order, model in [
    ('Z₃', 3, Z3_reps, 3, 'Clock'),
    ('S₃', 3, S3_reps, 6, 'Potts'),
]:
    multi_reps = multi_site_group_reps(single_reps, 2)
    dim_A = d**2  # 9
    n_tl, n_h, n_all = count_invariant_hermitian(multi_reps, dim_A)
    total_tl = dim_A*dim_A - 1  # 80
    ratio = n_tl / total_tl if total_tl > 0 else 0
    h_terms = count_hamiltonian_terms(2, d, model)
    h_ratio = h_terms / n_tl if n_tl > 0 else 0
    print(f"  {name:<8} n_A=2: {total_tl} total, {n_tl} G-inv, ratio={ratio:.4f}, H_terms={h_terms}, H/G-inv={h_ratio:.4f}")
    four_site_results[f"{name}_2site"] = {
        'd': d, 'dim_A': dim_A, 'group_order': str(order), 'model': model,
        'n_A': 2,
        'total_traceless_hermitian': total_tl,
        'invariant_traceless_hermitian': n_tl,
        'ratio_inv_to_total': ratio,
        'hamiltonian_terms': h_terms,
        'ratio_H_to_inv': h_ratio
    }

# Try n_A=3 for d=3 (dim_A=27, 729 operators — should be tractable)
print("\n--- d=3 models: computing for n_A=3 ---")
for name, d, single_reps, order, model in [
    ('Z₃', 3, Z3_reps, 3, 'Clock'),
    ('S₃', 3, S3_reps, 6, 'Potts'),
]:
    t1 = time.time()
    multi_reps = multi_site_group_reps(single_reps, 3)
    dim_A = d**3  # 27
    n_tl, n_h, n_all = count_invariant_hermitian(multi_reps, dim_A)
    total_tl = dim_A*dim_A - 1  # 728
    ratio = n_tl / total_tl if total_tl > 0 else 0
    h_terms = count_hamiltonian_terms(3, d, model)
    h_ratio = h_terms / n_tl if n_tl > 0 else 0
    dt = time.time() - t1
    print(f"  {name:<8} n_A=3: {total_tl} total, {n_tl} G-inv, ratio={ratio:.4f}, H_terms={h_terms}, H/G-inv={h_ratio:.4f} [{dt:.1f}s]")
    four_site_results[f"{name}_3site"] = {
        'd': d, 'dim_A': dim_A, 'group_order': str(order), 'model': model,
        'n_A': 3,
        'total_traceless_hermitian': total_tl,
        'invariant_traceless_hermitian': n_tl,
        'ratio_inv_to_total': ratio,
        'hamiltonian_terms': h_terms,
        'ratio_H_to_inv': h_ratio
    }

# ---- Summary and predictions ----
print("\n=== SUMMARY: Operator Counting vs Measured BW Locality ===")
print(f"{'Model':<12} {'Group':<8} {'d':<4} {'G-inv/total':<13} {'H/G-inv':<10} {'Measured BW':<12}")
print(f"{'TFIM':<12} {'Z₂':<8} {'2':<4} {'see above':<13} {'see above':<10} {'91%':<12}")
print(f"{'XXZ':<12} {'U(1)':<8} {'2':<4} {'see above':<13} {'see above':<10} {'100%':<12}")
print(f"{'Potts':<12} {'S₃':<8} {'3':<4} {'see above':<13} {'see above':<10} {'76.5%':<12}")
print(f"{'Clock':<12} {'Z₃':<8} {'3':<4} {'see above':<13} {'see above':<10} {'???':<12}")

print(f"\nPREDICTION: If H/G-inv correlates with measured BW locality,")
print(f"then Z₃ (smaller group, fewer constraints) should have LOWER BW locality than S₃.")

# ---- Save ----
results = {
    'experiment': '034a',
    'description': 'G-invariant operator counting for BW locality prediction',
    'single_site': single_site_results,
    'two_site': two_site_results,
    'multi_site': four_site_results,
    'runtime_seconds': time.time() - t0
}

with open('results/sprint_034a_operator_counting.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nSaved. Runtime: {time.time()-t0:.1f}s")
