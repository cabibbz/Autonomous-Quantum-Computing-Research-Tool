"""Sprint 060b: OPE coefficients for q=3 Potts at true g_c=1/3.

Uses ratio method validated in 060a: C_{sigma,sigma,epsilon} = |<eps|Z|sigma>| / |<0|Z|sigma>|
where Z is the clock operator at site 0.

q=3 Potts CFT (W_3 minimal model M(5,6)): c=4/5, x_sigma=2/15, x_epsilon=4/5.
Spin fields: (sigma, sigma^2) conjugate pair.
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron
from scipy.sparse.linalg import eigsh
import json, time

def tensor_op(ops_list):
    result = csr_matrix(ops_list[0])
    for op in ops_list[1:]:
        result = sp_kron(result, csr_matrix(op), format='csr')
    return result

def potts_hamiltonian_periodic(n, q, g):
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
        if j == i + 1:
            left_dim = q**i if i > 0 else 1
            right_dim = q**(n-j-1) if n-j-1 > 0 else 1
            op = csr_matrix(delta_2)
            if i > 0:
                op = sp_kron(csr_matrix(np.eye(left_dim)), op, format='csr')
            if n-j-1 > 0:
                op = sp_kron(op, csr_matrix(np.eye(right_dim)), format='csr')
            H = H - op
        else:
            for s in range(q):
                proj = np.zeros((q, q)); proj[s, s] = 1.0
                ops = [eye_q] * n; ops[0] = proj; ops[n-1] = proj
                H = H - tensor_op(ops)
    for i in range(n):
        ops = [eye_q] * n; ops[i] = XpXd
        H = H - g * tensor_op(ops)
    return H

def build_translation_operator(n, q):
    dim = q**n
    rows, cols = [], []
    for idx in range(dim):
        state = []
        tmp = idx
        for _ in range(n):
            state.append(tmp % q)
            tmp //= q
        new_state = state[1:] + [state[0]]
        new_idx = sum(s * q**k for k, s in enumerate(new_state))
        rows.append(new_idx)
        cols.append(idx)
    return csr_matrix((np.ones(dim), (rows, cols)), shape=(dim, dim))

def build_local_clock_op(n, q, site=0):
    """Z = diag(1, omega, omega^2, ...) at given site."""
    omega = np.exp(2j * np.pi / q)
    Z = np.diag([omega**s for s in range(q)])
    eye_q = np.eye(q)
    ops = [eye_q] * n
    ops[site] = Z
    return tensor_op(ops)

def get_momentum_states(evals, evecs, T, n):
    n_states = len(evals)
    momenta = np.zeros(n_states, dtype=int)
    tol = 1e-8 * abs(evals[-1] - evals[0]) if len(evals) > 1 else 1e-10
    i = 0
    while i < n_states:
        group = [i]
        j = i + 1
        while j < n_states and abs(evals[j] - evals[i]) < tol:
            group.append(j)
            j += 1
        sub = evecs[:, group]
        T_sub = (sub.conj().T) @ (T @ sub)
        t_evals = np.linalg.eigvals(T_sub)
        for idx_g, g_idx in enumerate(group):
            phase = np.angle(t_evals[idx_g])
            k = phase / (2 * np.pi / n)
            momenta[g_idx] = int(round(k)) % n
        i = j
    return momenta


def extract_ope_ratio(n_sites, q, g_c, R_eps_target, n_eig=30):
    """Extract OPE C_{sigma,sigma,epsilon} using ratio method.

    C = |<epsilon|Z_0|sigma>| / |<0|Z_0|sigma>|

    Also extract C_{sigma,sigma,sigma2} for higher spin fields.
    """
    t0 = time.time()
    dim = q**n_sites
    print(f"\n  n={n_sites} (dim={dim})...", end=" ", flush=True)

    H = potts_hamiltonian_periodic(n_sites, q, g_c)
    T = build_translation_operator(n_sites, q)
    Z_local = build_local_clock_op(n_sites, q, site=0)

    n_eig = min(n_eig, dim - 2)
    evals, evecs = eigsh(H, k=n_eig, which='SA')
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]

    momenta = get_momentum_states(evals, evecs, T, n_sites)
    gaps = evals - evals[0]
    delta1 = gaps[1] if gaps[1] > 1e-15 else gaps[2]
    ratios = gaps / delta1

    dt = time.time() - t0
    print(f"{dt:.1f}s")

    # Print spectrum
    print(f"  {'i':>3} {'R':>8} {'k':>4} {'gap':>12}")
    for i in range(min(20, len(ratios))):
        k = int(momenta[i])
        s = k if k <= n_sites//2 else k - n_sites
        print(f"  {i:>3} {ratios[i]:>8.4f} {s:>4} {gaps[i]:>12.6f}")

    # Compute ALL matrix elements of Z with low-lying states
    n_me = min(20, n_eig)
    me_Z = np.zeros((n_me, n_me), dtype=complex)
    for a in range(n_me):
        Z_b = Z_local @ evecs[:, :n_me]
        for a2 in range(n_me):
            me_Z[a2, :n_me] = evecs[:, a2].conj() @ Z_b

    # For q>=3, the first excited level is a DEGENERATE PAIR (sigma, sigma*)
    # The ground state is k=0. sigma has k=0 (for symmetric combo of sigma pair).
    # We need to identify the sigma pair.

    # Find sigma: first non-trivial k=0 state(s) with R~1
    sigma_indices = []
    for i in range(1, min(10, len(ratios))):
        if 0.5 < ratios[i] < 1.5:
            sigma_indices.append(i)

    print(f"\n  Sigma candidates: {sigma_indices} with R = {[round(float(ratios[i]),4) for i in sigma_indices]}")
    print(f"  Their momenta: {[int(momenta[i]) for i in sigma_indices]}")

    # For the ratio method, we need the LARGEST |<0|Z|sigma>|
    # among the sigma-level states
    best_sigma = None
    best_me = 0
    for si in sigma_indices:
        me = np.abs(me_Z[0, si])
        if me > best_me:
            best_me = me
            best_sigma = si

    if best_sigma is None:
        print("  ERROR: No sigma state found!")
        return None

    me_0_sigma = best_me
    print(f"\n  Best sigma: index {best_sigma}, R={ratios[best_sigma]:.4f}, |<0|Z|sigma>| = {me_0_sigma:.6f}")

    # Find epsilon: k=0 state near R_eps_target
    eps_idx = None
    best_eps_dist = 999
    for i in range(1, min(20, len(ratios))):
        k = int(momenta[i])
        s = k if k <= n_sites//2 else k - n_sites
        if s == 0 and abs(ratios[i] - R_eps_target) < best_eps_dist and i not in sigma_indices:
            best_eps_dist = abs(ratios[i] - R_eps_target)
            eps_idx = i

    # Extract ALL OPE coefficients via ratio
    # C_{a,sigma,b} = |<a|Z|b>| / |<0|Z|sigma>| when b=sigma
    # But more useful: for each state a, compute <a|Z|sigma_best> / <0|Z|sigma_best>

    print(f"\n  === OPE coefficients C(a, Z, sigma) = |<a|Z|sigma>| / |<0|Z|sigma>| ===")
    print(f"  {'a':>3} {'R_a':>8} {'k_a':>4} {'|<a|Z|sigma>|':>15} {'C_ratio':>10}")

    ope_results = []
    for a in range(min(20, n_me)):
        me = np.abs(me_Z[a, best_sigma])
        if me > 1e-8:
            k = int(momenta[a])
            s = k if k <= n_sites//2 else k - n_sites
            C = me / me_0_sigma
            label = ""
            if a == 0:
                label = " (Identity → normalization)"
            elif a == eps_idx:
                label = " *** EPSILON ***"
            elif a in sigma_indices:
                label = " (sigma self)"
            print(f"  {a:>3} {ratios[a]:>8.4f} {s:>4} {me:>15.6f} {C:>10.4f}{label}")
            ope_results.append({
                "a": int(a), "R_a": round(float(ratios[a]), 6),
                "k_a": int(s), "me_abs": round(float(me), 8),
                "C_ratio": round(float(C), 6),
                "is_epsilon": (a == eps_idx)
            })

    # Also extract using sigma* (Z^dag)
    # <a|Z^dag|sigma> — this probes the conjugate channel
    Z_dag = Z_local.conj().T
    me_Zd = np.zeros(n_me, dtype=complex)
    for a in range(n_me):
        me_Zd[a] = np.dot(evecs[:, a].conj(), Z_dag @ evecs[:, best_sigma])

    me_0_sigma_dag = np.abs(me_Zd[0])
    print(f"\n  |<0|Z^dag|sigma>| = {me_0_sigma_dag:.6f} (should equal |<0|Z|sigma>| = {me_0_sigma:.6f})")

    # For q=3: also check Z^2 operator (second harmonic)
    if q >= 3:
        Z2_local = Z_local @ Z_local  # Z^2 at site 0 (but need single-site Z^2)
        # Actually Z_local is already the full Hilbert space operator.
        # Z^2 at site 0 = diag(1, omega^2, omega^4, ...) at site 0
        omega = np.exp(2j * np.pi / q)
        Z2_single = np.diag([omega**(2*s) for s in range(q)])
        eye_q = np.eye(q)
        ops = [eye_q] * n_sites
        ops[0] = Z2_single
        Z2_local = tensor_op(ops)

        # sigma2 channel: <a|Z^2|sigma>
        print(f"\n  === Z^2 channel (probes sigma^2 field) ===")
        for si in sigma_indices:
            me_0_Z2_s = np.abs(np.dot(evecs[:, 0].conj(), Z2_local @ evecs[:, si]))
            if me_0_Z2_s > 1e-8:
                print(f"  |<0|Z^2|sigma[{si}]>| = {me_0_Z2_s:.6f}")

        # Also <a|Z^2|0> — probes the sigma^2 field from vacuum
        print(f"  OPE from Z^2:")
        for a in range(min(15, n_me)):
            me = np.abs(np.dot(evecs[:, a].conj(), Z2_local @ evecs[:, 0]))
            if me > 1e-6:
                k = int(momenta[a])
                s = k if k <= n_sites//2 else k - n_sites
                print(f"    |<{a}|Z^2|0>| = {me:.6f}, R={ratios[a]:.4f}, k={s}")

    result = {
        "q": q, "n": n_sites, "dim": dim, "g_c": g_c,
        "E0": float(evals[0]), "delta1": float(delta1),
        "sigma_idx": int(best_sigma),
        "R_sigma": float(ratios[best_sigma]),
        "me_0_Z_sigma": float(me_0_sigma),
        "epsilon_idx": int(eps_idx) if eps_idx is not None else None,
        "R_epsilon": float(ratios[eps_idx]) if eps_idx is not None else None,
        "C_sigma_sigma_epsilon": float(np.abs(me_Z[eps_idx, best_sigma]) / me_0_sigma) if eps_idx is not None else None,
        "ope_data": ope_results,
        "time_s": round(time.time() - t0, 2)
    }
    return result


# ===== q=3 Potts at g_c = 1/3 =====
print("=" * 70)
print("SPRINT 060b: OPE COEFFICIENTS — q=3 POTTS (W₃ MINIMAL MODEL)")
print("=" * 70)
print("CFT: c=4/5, x_sigma=2/15, x_epsilon=4/5")
print("Spin fields: 2 (conjugate pair sigma, sigma*)")
print("R_epsilon = x_eps/x_sigma = (4/5)/(2/15) = 6.0")

all_results = {}
# q=3: n up to 10 (3^10 = 59049)
for n in [6, 8, 10]:
    r = extract_ope_ratio(n, q=3, g_c=1/3, R_eps_target=6.0, n_eig=30)
    if r:
        all_results[f"n={n}"] = r

print("\n\n" + "=" * 70)
print("SUMMARY: C_{sigma,sigma,epsilon} for q=3")
print("=" * 70)
for key, res in sorted(all_results.items()):
    C = res.get("C_sigma_sigma_epsilon")
    if C is not None:
        print(f"  {key}: C = {C:.6f}, R_eps = {res['R_epsilon']:.4f}")

# Scaling of <0|Z|sigma>
print("\nScaling check: |<0|Z|sigma>| * N^{x_sigma}:")
for key, res in sorted(all_results.items()):
    n = res["n"]
    me = res["me_0_Z_sigma"]
    x_sigma = 2/15  # q=3
    scaled = me * n**x_sigma
    print(f"  {key}: |<0|Z|sigma>| = {me:.6f}, * N^(2/15) = {scaled:.6f}")

with open("results/sprint_060b_ope_q3.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print("\nResults saved.")
