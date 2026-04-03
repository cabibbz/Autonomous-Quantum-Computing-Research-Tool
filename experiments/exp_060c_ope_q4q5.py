"""Sprint 060c: OPE coefficients for q=4,5 Potts — entering novel CFT territory.

q=4 (Ashkin-Teller): c=1, x_sigma~0.117, R_eps~6.6. Logarithmic corrections expected.
q=5 (novel): c~1.10, x_sigma~0.101, R_eps~7.1. First OPE measurement for Hermitian q=5.

Reuses validated ratio method: C_{sigma,sigma,epsilon} = |<eps|Z|sigma>| / |<0|Z|sigma>|
Also extracts C_{sigma,sigma,sigma2} for higher harmonics.
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
    omega = np.exp(2j * np.pi / q)
    Z = np.diag([omega**s for s in range(q)])
    eye_q = np.eye(q)
    ops = [eye_q] * n
    ops[site] = Z
    return tensor_op(ops)

def build_local_clock_power(n, q, power, site=0):
    """Z^power at given site."""
    omega = np.exp(2j * np.pi / q)
    Zp = np.diag([omega**(power*s) for s in range(q)])
    eye_q = np.eye(q)
    ops = [eye_q] * n
    ops[site] = Zp
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

def extract_ope_multi(n_sites, q, g_c, R_eps_target, n_eig=30):
    """Extract multiple OPE coefficients using ratio method."""
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
    print(f"  {'i':>3} {'R':>8} {'k':>4}")
    for i in range(min(20, len(ratios))):
        k = int(momenta[i]); s = k if k <= n_sites//2 else k - n_sites
        marker = ""
        if 0.8 < ratios[i] < 1.2 and i > 0: marker = " <-- sigma"
        elif abs(ratios[i] - R_eps_target) < 1.5 and s == 0 and i > 0: marker = " <-- epsilon?"
        print(f"  {i:>3} {ratios[i]:>8.4f} {s:>4}{marker}")

    # Compute matrix elements of Z with low-lying states
    n_me = min(20, n_eig)
    Z_evecs = Z_local @ evecs[:, :n_me]
    me_Z = np.zeros((n_me, n_me), dtype=complex)
    for a in range(n_me):
        me_Z[a, :] = evecs[:, a].conj() @ Z_evecs

    # Find sigma: first excited level (possibly degenerate pair)
    sigma_indices = []
    for i in range(1, min(10, len(ratios))):
        if 0.5 < ratios[i] < 1.5:
            sigma_indices.append(i)

    # Pick sigma with largest |<0|Z|sigma>|
    best_sigma = max(sigma_indices, key=lambda i: np.abs(me_Z[0, i]))
    me_0_sigma = np.abs(me_Z[0, best_sigma])
    print(f"\n  Sigma: idx={best_sigma}, R={ratios[best_sigma]:.4f}, |<0|Z|sigma>|={me_0_sigma:.6f}")
    print(f"  Sigma pair: {sigma_indices}, |<0|Z|s>| = {[round(float(np.abs(me_Z[0,i])),6) for i in sigma_indices]}")

    # Find epsilon: k=0 state near R_eps_target, not in sigma
    eps_candidates = []
    for i in range(1, min(20, len(ratios))):
        k = int(momenta[i]); s = k if k <= n_sites//2 else k - n_sites
        if s == 0 and i not in sigma_indices and abs(ratios[i] - R_eps_target) < 2.0:
            eps_candidates.append((i, ratios[i]))
    eps_idx = eps_candidates[0][0] if eps_candidates else None
    print(f"  Epsilon: idx={eps_idx}, R={ratios[eps_idx]:.4f}" if eps_idx else "  Epsilon: NOT FOUND")

    # Extract OPE ratios
    print(f"\n  === OPE ratios C(a, sigma) = |<a|Z|sigma>| / |<0|Z|sigma>| ===")
    ope = {}
    for a in range(min(20, n_me)):
        me = np.abs(me_Z[a, best_sigma])
        if me > 1e-8:
            k = int(momenta[a]); s = k if k <= n_sites//2 else k - n_sites
            C = me / me_0_sigma
            label = "I" if a == 0 else ("eps" if a == eps_idx else ("sigma" if a in sigma_indices else f"L{a}"))
            ope[label] = {"idx": a, "R": float(ratios[a]), "k": s, "C": float(C)}
            if a < 12 or a == eps_idx:
                print(f"  {label:>8}: R={ratios[a]:>8.4f}, k={s:>3}, C={C:.4f}")

    # For q>=4: also compute Z^2 channel OPE (probes sigma^2 field)
    # sigma^2 pair: states at R corresponding to second harmonic
    if q >= 4:
        Z2_local = build_local_clock_power(n_sites, q, 2, site=0)
        Z2_evecs = Z2_local @ evecs[:, :n_me]
        me_Z2 = np.zeros((n_me, n_me), dtype=complex)
        for a in range(n_me):
            me_Z2[a, :] = evecs[:, a].conj() @ Z2_evecs

        # |<0|Z^2|sigma2>| where sigma2 is the second-harmonic field
        # sigma2 states: R near known R(sigma^2)
        R_sigma2_target = {4: 1.68, 5: 2.41, 7: 3.32, 10: 3.66}
        R_s2 = R_sigma2_target.get(q, 2.0)
        sigma2_candidates = []
        for i in range(1, min(20, len(ratios))):
            if abs(ratios[i] - R_s2) < 1.0 and i not in sigma_indices:
                sigma2_candidates.append(i)

        if sigma2_candidates:
            print(f"\n  === Z^2 channel (sigma^2 field, R~{R_s2}) ===")
            best_s2 = max(sigma2_candidates, key=lambda i: np.abs(me_Z2[0, i]))
            me_0_s2 = np.abs(me_Z2[0, best_s2])
            print(f"  sigma^2: idx={best_s2}, R={ratios[best_s2]:.4f}, |<0|Z^2|sigma^2>|={me_0_s2:.6f}")

            # C_{sigma2, sigma2, epsilon} = |<eps|Z^2|sigma2>| / |<0|Z^2|sigma2>|
            if eps_idx is not None and me_0_s2 > 1e-8:
                me_eps_s2 = np.abs(me_Z2[eps_idx, best_s2])
                C_s2_eps = me_eps_s2 / me_0_s2
                print(f"  C_{{sigma2,sigma2,epsilon}} = {C_s2_eps:.4f}")
                ope["sigma2_eps"] = {"C": float(C_s2_eps), "R_sigma2": float(ratios[best_s2])}

            # Cross-channel: C_{sigma,sigma2,...} from <sigma2|Z|sigma>
            me_s2_Z_sigma = np.abs(me_Z[best_s2, best_sigma])
            C_cross = me_s2_Z_sigma / me_0_sigma
            print(f"  C_{{sigma2,Z,sigma}} = {C_cross:.4f} (cross-channel)")
            ope["sigma_sigma2_cross"] = {"C": float(C_cross)}

    C_sse = float(np.abs(me_Z[eps_idx, best_sigma]) / me_0_sigma) if eps_idx is not None else None

    result = {
        "q": q, "n": n_sites, "dim": dim, "g_c": g_c,
        "E0": float(evals[0]), "delta1": float(delta1),
        "me_0_Z_sigma": float(me_0_sigma),
        "R_epsilon": float(ratios[eps_idx]) if eps_idx else None,
        "C_sigma_sigma_epsilon": C_sse,
        "ope": {k: v for k, v in ope.items()},
        "time_s": round(time.time() - t0, 2)
    }
    return result


# ===== q=4 Potts =====
print("=" * 70)
print("q=4 POTTS at g_c=0.392 — Ashkin-Teller point, c=1")
print("=" * 70)
results_q4 = {}
for n in [4, 6, 8]:  # 4^8 = 65536
    r = extract_ope_multi(n, q=4, g_c=0.392, R_eps_target=6.6, n_eig=30)
    if r:
        results_q4[f"n={n}"] = r

# ===== q=5 Potts =====
print("\n" + "=" * 70)
print("q=5 POTTS at g_c=0.441 — novel CFT, c~1.10")
print("=" * 70)
results_q5 = {}
for n in [4, 6, 8]:  # 5^8 = 390625
    r = extract_ope_multi(n, q=5, g_c=0.441, R_eps_target=7.1, n_eig=25)
    if r:
        results_q5[f"n={n}"] = r

# ===== SUMMARY =====
print("\n\n" + "=" * 70)
print("OPE SUMMARY: C_{sigma,sigma,epsilon}")
print("=" * 70)
print(f"  {'q':>3} {'n':>4} {'C_sse':>8} {'R_eps':>8}")
print(f"  {'-'*30}")

# Include q=2 and q=3 known values
print(f"  {'2':>3} {'12':>4} {'0.4979':>8} {'7.97':>8}")
print(f"  {'3':>3} {'10':>4} {'0.5156':>8} {'6.19':>8}")

for label, results in [("q=4", results_q4), ("q=5", results_q5)]:
    q = int(label.split("=")[1])
    for key, res in sorted(results.items()):
        C = res.get("C_sigma_sigma_epsilon")
        R = res.get("R_epsilon")
        if C is not None:
            print(f"  {q:>3} {res['n']:>4} {C:>8.4f} {R:>8.2f}")

all_results = {"q4": results_q4, "q5": results_q5}
with open("results/sprint_060c_ope_q4q5.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print("\nResults saved.")
