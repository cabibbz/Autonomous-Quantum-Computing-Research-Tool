"""Sprint 101a: Symmetry-Resolved Entanglement Entropy (SREE) for S_q Potts.

Decompose ρ_A into Z_q charge sectors via subsystem symmetry operator G_A.
Measure S(α) for each charge sector α = 0, 1, ..., q-1.
Test equipartition: S(α) ≈ S_total / q at leading order.

Uses periodic BC at g_c = 1/q. Midchain bipartition nA = n/2.
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from gpu_utils import eigsh
import json, time

def potts_hamiltonian_periodic(n, q, g):
    """H = -Σ δ(s_i,s_j) - g Σ (X + X†)_i, periodic BC."""
    dim = q**n
    eye_q = np.eye(q)
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = X + X.T
    delta_2 = np.zeros((q**2, q**2))
    for s in range(q):
        delta_2[s*q+s, s*q+s] = 1.0

    H = csr_matrix((dim, dim))
    for i in range(n):
        j = (i + 1) % n
        # Build delta(s_i, s_j)
        ops_before = q**i if i > 0 else 1
        if j == i + 1:
            op = csr_matrix(delta_2)
            if i > 0:
                op = sp_kron(sp_eye(ops_before), op, format='csr')
            remaining = q**(n - j - 1) if n - j - 1 > 0 else 1
            if n - j - 1 > 0:
                op = sp_kron(op, sp_eye(remaining), format='csr')
            H = H - op
        else:
            # Wrap-around bond: sites n-1 and 0
            for s in range(q):
                proj = np.zeros((q, q))
                proj[s, s] = 1.0
                ops = [eye_q]*n
                ops[0] = proj
                ops[n-1] = proj
                r = csr_matrix(ops[0])
                for op in ops[1:]:
                    r = sp_kron(r, csr_matrix(op), format='csr')
                H = H - r

    for i in range(n):
        ops = [eye_q]*n
        ops[i] = XpXd
        r = csr_matrix(ops[0])
        for op in ops[1:]:
            r = sp_kron(r, csr_matrix(op), format='csr')
        H = H - g * r
    return H

def get_ground_state(H, dim):
    """Get ground state vector."""
    evals, evecs = eigsh(H, k=1, which='SA')
    return evals[0], evecs[:, 0]

def reduced_density_matrix(psi, n, q, nA):
    """Compute ρ_A by tracing out sites nA..n-1."""
    dimA = q**nA
    dimB = q**(n - nA)
    psi_mat = psi.reshape(dimB, dimA).T  # shape (dimA, dimB) — sites 0..nA-1 are subsystem A
    # Actually need to be careful with ordering. State index = s0 + q*s1 + q^2*s2 + ...
    # So first index (fastest) = site 0. Reshape as (q^nA, q^(n-nA)) with A = sites 0..nA-1
    psi_mat = psi.reshape((dimA, dimB), order='C')  # C order: last index fastest...
    # Wait: index = s0 + q*s1 + ... + q^(n-1)*s_{n-1}
    # So psi[s0 + q*s1 + ... + q^(nA-1)*s_{nA-1} + q^nA * (s_nA + q*s_{nA+1} + ...)]
    # = psi[idxA + dimA * idxB]
    # So reshape as (dimA, dimB) with Fortran order, or (dimB, dimA) with C order
    psi_mat = psi.reshape((dimB, dimA)).T  # now psi_mat[idxA, idxB]
    # Wait, let me think again. psi[i] where i = idxA + dimA * idxB
    # So psi.reshape(dimA, dimB, order='C') gives psi_mat[idxA, idxB] — NO
    # C order: last index varies fastest. psi.reshape((dimA, dimB)) means
    # psi_mat[a, b] = psi[a * dimB + b]. But we want psi[idxA + dimA * idxB] = psi_mat[idxA, idxB]
    # So a*dimB + b = idxA + dimA * idxB doesn't match unless we use Fortran order.
    psi_mat = psi.reshape((dimA, dimB), order='F')
    # Check: psi_mat[idxA, idxB] = psi[idxA + dimA*idxB] ✓ (Fortran: first index varies fastest)
    rho_A = psi_mat @ psi_mat.conj().T
    return rho_A

def subsystem_symmetry_projector(nA, q, charge):
    """Project onto charge sector α of subsystem A.

    G_A = X_0 ⊗ X_1 ⊗ ... ⊗ X_{nA-1} has eigenvalues ω^α where ω = e^{2πi/q}.
    Projector P_α = (1/q) Σ_{k=0}^{q-1} ω^{-αk} G_A^k.
    """
    dimA = q**nA
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0

    # Build G_A = X⊗X⊗...⊗X
    GA = X.copy()
    for _ in range(nA - 1):
        GA = np.kron(GA, X)

    # Build projector: P_α = (1/q) Σ_k ω^{-αk} G_A^k
    omega = np.exp(2j * np.pi / q)
    P = np.zeros((dimA, dimA), dtype=complex)
    GA_k = np.eye(dimA)  # G_A^0 = I
    for k in range(q):
        P += omega**(-charge * k) * GA_k
        GA_k = GA_k @ GA  # G_A^{k+1}
    P /= q
    return P

def symmetry_resolved_entropy(rho_A, nA, q):
    """Compute SREE: S(α) for each Z_q charge sector.

    Returns dict with S_total, S(α) for each α, probabilities p(α),
    configurational entropy S_config, and number entropy S_number.
    """
    dimA = q**nA
    results = {}

    # Total entropy
    evals_total = np.linalg.eigvalsh(rho_A)
    evals_total = evals_total[evals_total > 1e-15]
    S_total = -np.sum(evals_total * np.log(evals_total))
    results['S_total'] = float(S_total)

    # Project into each charge sector
    p_alpha = []  # probability in each sector
    S_alpha = []  # entropy in each sector

    for alpha in range(q):
        P = subsystem_symmetry_projector(nA, q, alpha)
        # ρ_α = P_α ρ_A P_α (block of ρ_A in sector α)
        rho_alpha = P @ rho_A @ P.conj().T
        # Probability in this sector
        p = np.real(np.trace(rho_alpha))
        p_alpha.append(float(p))

        # Entropy of normalized block
        if p > 1e-15:
            rho_norm = rho_alpha / p
            evals_sec = np.linalg.eigvalsh(rho_norm)
            evals_sec = evals_sec[evals_sec > 1e-15]
            S_sec = -np.sum(evals_sec * np.log(evals_sec))
            S_alpha.append(float(S_sec))
        else:
            S_alpha.append(0.0)

    results['p_alpha'] = p_alpha
    results['S_alpha'] = S_alpha

    # Configurational (number) entropy: S_n = -Σ p(α) ln p(α)
    p_arr = np.array(p_alpha)
    p_pos = p_arr[p_arr > 1e-15]
    S_number = -np.sum(p_pos * np.log(p_pos))
    results['S_number'] = float(S_number)

    # Check decomposition: S_total = Σ p(α) S(α) + S_number
    S_config = sum(p * s for p, s in zip(p_alpha, S_alpha))
    results['S_config'] = float(S_config)
    results['S_reconstruct'] = float(S_config + S_number)
    results['reconstruction_error'] = float(abs(S_total - S_config - S_number))

    # Equipartition measure: deviation from uniform p(α) = 1/q
    equi_dev = np.std(p_alpha) / np.mean(p_alpha)  # CV of sector probabilities
    results['equipartition_CV'] = float(equi_dev)

    # Entropy equipartition: deviation of S(α) from mean
    S_arr = np.array(S_alpha)
    S_mean = np.mean(S_arr)
    if S_mean > 0:
        entropy_equi_CV = np.std(S_arr) / S_mean
    else:
        entropy_equi_CV = 0.0
    results['entropy_equi_CV'] = float(entropy_equi_CV)

    return results

# ============================================================
# Main experiment: SREE for q=2,3,4,5,7 at g_c = 1/q
# ============================================================

print("=" * 70)
print("Sprint 101a: Symmetry-Resolved Entanglement Entropy")
print("=" * 70)

all_results = {}

configs = [
    # (q, n, nA) — test timing first with small case
    (2, 8, 4),
    (3, 8, 4),
    (4, 6, 3),
    (5, 8, 4),
    (7, 6, 3),
]

for q, n, nA in configs:
    g_c = 1.0 / q
    dim = q**n
    print(f"\n{'='*60}")
    print(f"q={q}, n={n}, nA={nA}, dim={dim}, g_c={g_c:.4f}")
    print(f"{'='*60}")

    t0 = time.time()
    H = potts_hamiltonian_periodic(n, q, g_c)
    t_build = time.time() - t0
    print(f"  H built in {t_build:.1f}s")

    E0, psi = get_ground_state(H, dim)
    t_diag = time.time() - t0
    print(f"  Ground state E0={E0:.8f} in {t_diag:.1f}s")

    rho_A = reduced_density_matrix(psi, n, q, nA)
    t_rho = time.time() - t0
    print(f"  ρ_A ({q**nA}×{q**nA}) in {t_rho:.1f}s")

    sree = symmetry_resolved_entropy(rho_A, nA, q)
    t_sree = time.time() - t0
    print(f"  SREE computed in {t_sree:.1f}s")

    # Print results
    print(f"\n  S_total = {sree['S_total']:.6f}")
    print(f"  S_number (config) = {sree['S_number']:.6f}")
    print(f"  S_config (weighted) = {sree['S_config']:.6f}")
    print(f"  S_reconstruct = {sree['S_reconstruct']:.6f} (error: {sree['reconstruction_error']:.2e})")

    print(f"\n  Charge sector probabilities p(α):")
    for alpha in range(q):
        print(f"    α={alpha}: p={sree['p_alpha'][alpha]:.6f}, S={sree['S_alpha'][alpha]:.6f}")

    print(f"\n  Equipartition CV (probability): {sree['equipartition_CV']:.6f}")
    print(f"  Entropy equipartition CV: {sree['entropy_equi_CV']:.6f}")

    # Check conjugate pairing: S(α) should equal S(q-α) for Z_q
    if q > 2:
        print(f"\n  Conjugate pairing check (S(α) vs S(q-α)):")
        for alpha in range(1, (q+1)//2):
            conj = q - alpha
            diff = abs(sree['S_alpha'][alpha] - sree['S_alpha'][conj])
            print(f"    S({alpha})={sree['S_alpha'][alpha]:.6f}, S({conj})={sree['S_alpha'][conj]:.6f}, diff={diff:.2e}")

    key = f"q{q}_n{n}_nA{nA}"
    all_results[key] = {
        'q': q, 'n': n, 'nA': nA, 'dim': dim, 'g_c': g_c,
        'E0': float(E0), 'time_s': round(time.time() - t0, 2),
        **sree
    }

# Summary table
print("\n\n" + "=" * 70)
print("SUMMARY: Symmetry-Resolved Entanglement Entropy")
print("=" * 70)
print(f"\n  {'q':>2} {'n':>3} {'nA':>3} {'S_total':>8} {'S_number':>9} {'S_num/S_tot':>11} {'p_CV':>8} {'S_CV':>8}")
print(f"  {'-'*60}")
for key, r in sorted(all_results.items()):
    ratio = r['S_number'] / r['S_total'] if r['S_total'] > 0 else 0
    print(f"  {r['q']:>2} {r['n']:>3} {r['nA']:>3} {r['S_total']:>8.4f} {r['S_number']:>9.4f} {ratio:>11.4f} {r['equipartition_CV']:>8.5f} {r['entropy_equi_CV']:>8.5f}")

with open("results/sprint_101a_sree.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print("\nResults saved to results/sprint_101a_sree.json")
