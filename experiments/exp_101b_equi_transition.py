"""Sprint 101b: Equipartition transition — dense q scan + hybrid model comparison.

Map p(α=0) and equipartition CV across q=2-10 for S_q Potts and hybrid model.
Test whether the transition aligns with q=4 (real→complex CFT boundary).
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from gpu_utils import eigsh
import json, time

def potts_hamiltonian_periodic(n, q, g, model='S_q'):
    """Build Hamiltonian. model='S_q' uses Σ X^k field, 'hybrid' uses X+X†."""
    dim = q**n
    eye_q = np.eye(q)

    # Clock operator X
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0

    # Transverse field
    if model == 'hybrid':
        field_op = X + X.T  # X + X†
    else:  # S_q: Σ_{k=1}^{q-1} X^k
        field_op = np.zeros((q, q))
        Xk = np.eye(q)
        for k in range(1, q):
            Xk = Xk @ X
            field_op += Xk

    # Potts coupling δ(s_i, s_j)
    delta_2 = np.zeros((q**2, q**2))
    for s in range(q):
        delta_2[s*q+s, s*q+s] = 1.0

    H = csr_matrix((dim, dim))
    for i in range(n):
        j = (i + 1) % n
        if j == i + 1:
            left = q**i if i > 0 else 1
            right = q**(n-j-1) if n-j-1 > 0 else 1
            op = csr_matrix(delta_2)
            if i > 0:
                op = sp_kron(sp_eye(left), op, format='csr')
            if n-j-1 > 0:
                op = sp_kron(op, sp_eye(right), format='csr')
            H = H - op
        else:
            for s in range(q):
                proj = np.zeros((q, q)); proj[s, s] = 1.0
                ops = [eye_q]*n; ops[0] = proj; ops[n-1] = proj
                r = csr_matrix(ops[0])
                for op in ops[1:]:
                    r = sp_kron(r, csr_matrix(op), format='csr')
                H = H - r

    for i in range(n):
        ops = [eye_q]*n; ops[i] = field_op
        r = csr_matrix(ops[0])
        for op in ops[1:]:
            r = sp_kron(r, csr_matrix(op), format='csr')
        H = H - g * r
    return H

def reduced_density_matrix(psi, n, q, nA):
    dimA = q**nA
    dimB = q**(n - nA)
    psi_mat = psi.reshape((dimA, dimB), order='F')
    return psi_mat @ psi_mat.conj().T

def subsystem_charge_projectors(nA, q):
    """Return list of q projectors P_α for α=0,...,q-1."""
    dimA = q**nA
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    GA = X.copy()
    for _ in range(nA - 1):
        GA = np.kron(GA, X)

    omega = np.exp(2j * np.pi / q)
    projectors = []
    GA_k = np.eye(dimA)
    GA_powers = [GA_k.copy()]
    for k in range(1, q):
        GA_k = GA_k @ GA
        GA_powers.append(GA_k.copy())

    for alpha in range(q):
        P = np.zeros((dimA, dimA), dtype=complex)
        for k in range(q):
            P += omega**(-alpha * k) * GA_powers[k]
        P /= q
        projectors.append(P)
    return projectors

def compute_sree(rho_A, projectors, q):
    """Compute SREE given precomputed projectors."""
    # Total entropy
    evals_total = np.linalg.eigvalsh(rho_A)
    evals_total = evals_total[evals_total > 1e-15]
    S_total = -np.sum(evals_total * np.log(evals_total))

    p_alpha = []
    S_alpha = []
    for alpha in range(q):
        P = projectors[alpha]
        rho_alpha = P @ rho_A @ P.conj().T
        p = np.real(np.trace(rho_alpha))
        p_alpha.append(float(p))
        if p > 1e-15:
            rho_norm = rho_alpha / p
            evals_sec = np.linalg.eigvalsh(rho_norm)
            evals_sec = evals_sec[evals_sec > 1e-15]
            S_alpha.append(float(-np.sum(evals_sec * np.log(evals_sec))))
        else:
            S_alpha.append(0.0)

    p_arr = np.array(p_alpha)
    p_pos = p_arr[p_arr > 1e-15]
    S_number = float(-np.sum(p_pos * np.log(p_pos)))
    S_config = sum(p * s for p, s in zip(p_alpha, S_alpha))

    equi_cv = np.std(p_alpha) / np.mean(p_alpha) if np.mean(p_alpha) > 0 else 0

    return {
        'S_total': float(S_total),
        'S_number': float(S_number),
        'S_config': float(S_config),
        'p_alpha': p_alpha,
        'S_alpha': S_alpha,
        'p0': p_alpha[0],
        'equi_cv': float(equi_cv),
        'sn_ratio': float(S_number / S_total) if S_total > 0 else 0,
    }

# ============================================================
# Dense q scan: S_q Potts at n=6, nA=3
# ============================================================
print("=" * 70)
print("Sprint 101b: Equipartition Transition — Dense q Scan")
print("=" * 70)

all_results = {}

# S_q Potts: q=2-10 at n=6, nA=3 (uniform comparison)
n_sites = 6
nA = 3

for q in range(2, 11):
    g_c = 1.0 / q
    dim = q**n_sites
    if dim > 2_000_000:
        print(f"  q={q}: dim={dim} too large for n={n_sites}, skipping")
        continue

    t0 = time.time()
    H = potts_hamiltonian_periodic(n_sites, q, g_c, model='S_q')
    E0, psi = eigsh(H, k=1, which='SA')
    E0 = E0[0]; psi = psi[:, 0]
    rho_A = reduced_density_matrix(psi, n_sites, q, nA)
    projs = subsystem_charge_projectors(nA, q)
    sree = compute_sree(rho_A, projs, q)
    elapsed = time.time() - t0

    print(f"  q={q:>2} S_q:  p0={sree['p0']:.4f}, CV={sree['equi_cv']:.5f}, "
          f"S_n/S_t={sree['sn_ratio']:.4f}, S_t={sree['S_total']:.4f} ({elapsed:.1f}s)")

    key = f"Sq_q{q}_n{n_sites}"
    all_results[key] = {'q': q, 'n': n_sites, 'nA': nA, 'model': 'S_q', 'g_c': g_c,
                        'time_s': round(elapsed, 2), **sree}

# Hybrid model: q=2-8 at n=6, nA=3
# Need hybrid g_c values (from KNOWLEDGE.md)
hybrid_gc = {2: 0.5, 3: 0.333, 4: 0.392, 5: 0.441, 6: 0.475, 7: 0.504, 8: 0.527}

for q in sorted(hybrid_gc.keys()):
    if q > 8:
        continue
    g_c = hybrid_gc[q]
    dim = q**n_sites
    if dim > 2_000_000:
        print(f"  q={q}: dim={dim} too large for hybrid n={n_sites}, skipping")
        continue

    t0 = time.time()
    H = potts_hamiltonian_periodic(n_sites, q, g_c, model='hybrid')
    E0, psi = eigsh(H, k=1, which='SA')
    E0 = E0[0]; psi = psi[:, 0]
    rho_A = reduced_density_matrix(psi, n_sites, q, nA)
    projs = subsystem_charge_projectors(nA, q)
    sree = compute_sree(rho_A, projs, q)
    elapsed = time.time() - t0

    print(f"  q={q:>2} hyb:  p0={sree['p0']:.4f}, CV={sree['equi_cv']:.5f}, "
          f"S_n/S_t={sree['sn_ratio']:.4f}, S_t={sree['S_total']:.4f} ({elapsed:.1f}s)")

    key = f"hyb_q{q}_n{n_sites}"
    all_results[key] = {'q': q, 'n': n_sites, 'nA': nA, 'model': 'hybrid', 'g_c': g_c,
                        'time_s': round(elapsed, 2), **sree}

# Summary
print("\n\n" + "=" * 70)
print("EQUIPARTITION TRANSITION SUMMARY")
print("=" * 70)
print(f"\n  {'model':>6} {'q':>2} {'p(0)':>8} {'1/q':>8} {'p0*q':>6} {'CV':>8} {'S_n/S_t':>8}")
print(f"  {'-'*55}")
for key in sorted(all_results.keys()):
    r = all_results[key]
    mdl = r['model'][:3]
    p0q = r['p0'] * r['q']
    print(f"  {mdl:>6} {r['q']:>2} {r['p0']:>8.5f} {1/r['q']:>8.5f} {p0q:>6.3f} {r['equi_cv']:>8.5f} {r['sn_ratio']:>8.5f}")

# Check: does p(0)*q approach 1 (perfect equipartition) at large q?
print("\n\n  Equipartition quality: p(0)*q (=1 is perfect)")
for mdl in ['S_q', 'hybrid']:
    print(f"\n  {mdl}:")
    for key in sorted(all_results.keys()):
        r = all_results[key]
        if r['model'] != mdl:
            continue
        p0q = r['p0'] * r['q']
        deviation = abs(p0q - 1)
        bar = '#' * int(50 * min(deviation, 1))
        print(f"    q={r['q']:>2}: p0*q={p0q:.4f} (dev={deviation:.4f}) {bar}")

with open("results/sprint_101b_equi_transition.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print("\nResults saved to results/sprint_101b_equi_transition.json")
