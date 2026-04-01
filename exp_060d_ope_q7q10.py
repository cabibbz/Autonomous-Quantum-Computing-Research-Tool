"""Sprint 060d: OPE coefficients for q=7,10 Potts — deep in the novel CFT regime.

q=7: c~1.30, x_sigma~0.086, R_eps~8.1. Sizes n=4,6 (7^6=117649).
q=10: c~1.40, x_sigma~0.083, R_eps~8.3. Sizes n=4,6 (10^6=1M).

Also computes a comprehensive summary across all q values with 1/N extrapolation.
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

def extract_ope(n_sites, q, g_c, R_eps_target, n_eig=25):
    """Extract OPE coefficients."""
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
        print(f"  {i:>3} {ratios[i]:>8.4f} {s:>4}")

    # Matrix elements of Z
    n_me = min(20, n_eig)
    Z_evecs = Z_local @ evecs[:, :n_me]
    me_Z = np.zeros((n_me, n_me), dtype=complex)
    for a in range(n_me):
        me_Z[a, :] = evecs[:, a].conj() @ Z_evecs

    # Find sigma
    sigma_indices = [i for i in range(1, min(10, len(ratios))) if 0.5 < ratios[i] < 1.5]
    best_sigma = max(sigma_indices, key=lambda i: np.abs(me_Z[0, i]))
    me_0_sigma = np.abs(me_Z[0, best_sigma])
    print(f"\n  Sigma: idx={best_sigma}, R={ratios[best_sigma]:.4f}, |<0|Z|sigma>|={me_0_sigma:.6f}")

    # Find epsilon
    eps_candidates = []
    for i in range(1, min(20, len(ratios))):
        k = int(momenta[i]); s = k if k <= n_sites//2 else k - n_sites
        if s == 0 and i not in sigma_indices and abs(ratios[i] - R_eps_target) < 2.0:
            eps_candidates.append((i, ratios[i]))
    eps_idx = eps_candidates[0][0] if eps_candidates else None

    # OPE ratios
    print(f"\n  Key OPE ratios:")
    ope = {}
    for a in range(min(15, n_me)):
        me = np.abs(me_Z[a, best_sigma])
        if me > 1e-8:
            k = int(momenta[a]); s = k if k <= n_sites//2 else k - n_sites
            C = me / me_0_sigma
            if a == 0 or a == eps_idx or a in sigma_indices or (s == 0 and C > 0.1):
                label = "I" if a == 0 else ("eps" if a == eps_idx else ("sigma" if a in sigma_indices else f"L{a}"))
                print(f"  {label:>8}: R={ratios[a]:>8.4f}, k={s:>3}, C={C:.4f}")
                ope[label] = {"R": float(ratios[a]), "k": s, "C": float(C)}

    C_sse = float(np.abs(me_Z[eps_idx, best_sigma]) / me_0_sigma) if eps_idx is not None else None

    # Higher harmonic channels
    for power in range(2, min(q//2 + 1, 4)):
        Zp_local = build_local_clock_power(n_sites, q, power, site=0)
        Zp_evecs = Zp_local @ evecs[:, :n_me]
        me_Zp_0 = np.array([np.dot(evecs[:, 0].conj(), Zp_evecs[:, i]) for i in range(n_me)])

        # Find sigma^power state
        R_target = {
            (4, 2): 1.68, (5, 2): 2.41,
            (7, 2): 3.32, (7, 3): 5.95,
            (10, 2): 3.66, (10, 3): 7.69
        }.get((q, power), power * 1.5)

        sp_candidates = [(i, ratios[i]) for i in range(1, min(20, len(ratios)))
                         if abs(ratios[i] - R_target) < 1.5 and i not in sigma_indices]
        if sp_candidates:
            best_sp = max(sp_candidates, key=lambda x: np.abs(me_Zp_0[x[0]]))
            me_0_sp = np.abs(me_Zp_0[best_sp[0]])
            if me_0_sp > 1e-8:
                # C_{sp, sp, epsilon}
                if eps_idx is not None:
                    me_eps_sp = np.abs(np.dot(evecs[:, eps_idx].conj(), Zp_evecs[:, best_sp[0]]))
                    C_sp_eps = me_eps_sp / me_0_sp
                    print(f"  sigma^{power}: R={best_sp[1]:.4f}, |<0|Z^{power}|sigma^{power}>|={me_0_sp:.4f}, C_{{s{power}s{power}eps}}={C_sp_eps:.4f}")
                    ope[f"sigma{power}_eps"] = {"C": float(C_sp_eps), "R": float(best_sp[1])}

    result = {
        "q": q, "n": n_sites, "dim": dim, "g_c": g_c,
        "me_0_Z_sigma": float(me_0_sigma),
        "C_sigma_sigma_epsilon": C_sse,
        "R_epsilon": float(ratios[eps_idx]) if eps_idx else None,
        "ope": ope,
        "time_s": round(time.time() - t0, 2)
    }
    return result


# ===== q=7 Potts =====
print("=" * 70)
print("q=7 POTTS at g_c=0.535 — novel CFT, c~1.30")
print("=" * 70)
results_q7 = {}
for n in [4, 6]:  # 7^6=117649
    r = extract_ope(n, q=7, g_c=0.535, R_eps_target=8.1, n_eig=25)
    if r: results_q7[f"n={n}"] = r

# ===== q=10 Potts =====
print("\n" + "=" * 70)
print("q=10 POTTS at g_c=0.684 — novel CFT, c~1.40")
print("=" * 70)
results_q10 = {}
for n in [4, 6]:  # 10^6=1M
    ne = 15 if n == 6 else 25
    r = extract_ope(n, q=10, g_c=0.684, R_eps_target=8.3, n_eig=ne)
    if r: results_q10[f"n={n}"] = r

# ===== COMPREHENSIVE SUMMARY =====
print("\n\n" + "=" * 70)
print("COMPREHENSIVE OPE SUMMARY ACROSS ALL q")
print("=" * 70)

# Collect all data
all_data = [
    (2, 12, 0.4979, 7.97),
    (2, 10, 0.4969, 7.95),
    (2, 8, 0.4951, 7.92),
    (3, 10, 0.5156, 6.19),
    (3, 8, 0.5095, 6.22),
    (3, 6, 0.4995, 6.24),
]

for q_label, results in [("q=4", {}), ("q=5", {}), ("q=7", results_q7), ("q=10", results_q10)]:
    q = int(q_label.split("=")[1])
    for key, res in results.items():
        C = res.get("C_sigma_sigma_epsilon")
        R = res.get("R_epsilon")
        if C is not None:
            all_data.append((q, res["n"], C, R))

# Add q=4,5 from 060c
all_data.extend([
    (4, 8, 0.4301, 6.58),
    (4, 6, 0.4195, 6.53),
    (4, 4, 0.3964, 6.44),
    (5, 8, 0.3491, 7.23),
    (5, 6, 0.3411, 7.14),
    (5, 4, 0.3224, 6.98),
])

print(f"\n  {'q':>3} {'n':>4} {'C_sse':>8} {'R_eps':>8}")
print(f"  {'-'*30}")
for q, n, C, R in sorted(all_data, key=lambda x: (x[0], -x[1])):
    print(f"  {q:>3} {n:>4} {C:>8.4f} {R:>8.2f}")

# 1/N extrapolation for each q
print(f"\n\n  === 1/N EXTRAPOLATION of C_{{sigma,sigma,epsilon}} ===")
print(f"  {'q':>3} {'C(largest n)':>14} {'C(N→∞ est.)':>14} {'method':>20}")
by_q = {}
for q, n, C, R in all_data:
    if q not in by_q:
        by_q[q] = []
    by_q[q].append((n, C))

for q in sorted(by_q.keys()):
    pairs = sorted(by_q[q], key=lambda x: x[0])
    n_max, C_max = pairs[-1]
    if len(pairs) >= 2:
        n1, C1 = pairs[-2]
        n2, C2 = pairs[-1]
        # Linear extrapolation in 1/N
        inv_n1, inv_n2 = 1/n1, 1/n2
        slope = (C2 - C1) / (inv_n2 - inv_n1)
        C_inf = C2 - slope * inv_n2
        method = f"linear ({n1},{n2})"
    else:
        C_inf = C_max
        method = "single point"
    print(f"  {q:>3} {C_max:>14.4f} {C_inf:>14.4f} {method:>20}")

all_results = {"q7": results_q7, "q10": results_q10}
with open("results/sprint_060d_ope_q7q10.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print("\nResults saved.")
