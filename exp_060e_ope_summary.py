"""Sprint 060e: Corrected OPE summary with proper epsilon identification.

The bug in 060d: for q=10, sigma^3 pair (R~7.86, Z_q charge 3) was misidentified
as epsilon (Z_q charge 0). Fix: identify epsilon as the k=0 state with NONZERO
|<state|Z|sigma>|, since Z_q selection rules force <sigma^k|Z|sigma*>=0 for k!=0.

This script re-extracts q=10 with corrected logic, and produces the final
comprehensive OPE table with 1/N extrapolation.
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

def build_zq_symmetry(n, q):
    """Global Z_q symmetry operator S = product(X_i).
    Cyclic shift of ALL spins: |s_0,...,s_{n-1}> -> |s_0+1,...,s_{n-1}+1 mod q>."""
    dim = q**n
    rows, cols = [], []
    for idx in range(dim):
        state = []
        tmp = idx
        for _ in range(n):
            state.append(tmp % q)
            tmp //= q
        new_state = [(s + 1) % q for s in state]
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

def get_quantum_numbers(evals, evecs, T, S, n, q):
    """Get both momentum k and Z_q charge for each eigenstate."""
    n_states = len(evals)
    momenta = np.zeros(n_states, dtype=int)
    charges = np.zeros(n_states, dtype=int)
    omega = np.exp(2j * np.pi / q)

    # Momentum from T
    for i in range(n_states):
        t_val = np.dot(evecs[:, i].conj(), T @ evecs[:, i])
        phase = np.angle(t_val)
        k = phase / (2 * np.pi / n)
        momenta[i] = int(round(k)) % n

    # Z_q charge from S
    for i in range(n_states):
        s_val = np.dot(evecs[:, i].conj(), S @ evecs[:, i])
        phase = np.angle(s_val)
        ch = phase / (2 * np.pi / q)
        charges[i] = int(round(ch)) % q

    return momenta, charges


def extract_ope_corrected(n_sites, q, g_c, n_eig=25):
    """Extract OPE with corrected epsilon identification using Z_q charge."""
    t0 = time.time()
    dim = q**n_sites
    print(f"\n  q={q}, n={n_sites} (dim={dim})...", end=" ", flush=True)

    H = potts_hamiltonian_periodic(n_sites, q, g_c)
    T = build_translation_operator(n_sites, q)
    S = build_zq_symmetry(n_sites, q)
    Z_local = build_local_clock_op(n_sites, q, site=0)

    n_eig = min(n_eig, dim - 2)
    evals, evecs = eigsh(H, k=n_eig, which='SA')
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]

    momenta, charges = get_quantum_numbers(evals, evecs, T, S, n_sites, q)
    gaps = evals - evals[0]
    delta1 = gaps[1] if gaps[1] > 1e-15 else gaps[2]
    ratios = gaps / delta1

    dt = time.time() - t0
    print(f"{dt:.1f}s")

    # Print spectrum with charge
    print(f"  {'i':>3} {'R':>8} {'k':>4} {'ch':>4} {'label':>12}")
    for i in range(min(25, len(ratios))):
        k = int(momenta[i]); sk = k if k <= n_sites//2 else k - n_sites
        ch = int(charges[i])
        label = ""
        if i == 0: label = "Identity"
        elif ch == 0 and sk == 0: label = "charge-0"
        elif abs(ratios[i] - 1.0) < 0.2: label = f"sigma(ch={ch})"
        else: label = f"ch={ch}"
        print(f"  {i:>3} {ratios[i]:>8.4f} {sk:>4} {ch:>4} {label:>12}")

    # Matrix elements of Z
    n_me = min(25, n_eig)
    Z_evecs = Z_local @ evecs[:, :n_me]
    me_Z = np.zeros((n_me, n_me), dtype=complex)
    for a in range(n_me):
        me_Z[a, :] = evecs[:, a].conj() @ Z_evecs

    # Sigma: first excited level (charge 1 or q-1)
    sigma_indices = [i for i in range(1, min(10, len(ratios))) if 0.5 < ratios[i] < 1.5]
    # Pick the one with charge q-1 (sigma*) since Z couples 0 → q-1
    sigma_star = [i for i in sigma_indices if int(charges[i]) == q-1]
    sigma_one = [i for i in sigma_indices if int(charges[i]) == 1]

    # For degenerate pair, eigsh may give superposition. Use the one with largest <0|Z|sigma>
    best_sigma = max(sigma_indices, key=lambda i: np.abs(me_Z[0, i]))
    me_0_sigma = np.abs(me_Z[0, best_sigma])
    print(f"\n  Best sigma: idx={best_sigma}, charge={int(charges[best_sigma])}, |<0|Z|sigma>|={me_0_sigma:.6f}")
    print(f"  Sigma pair charges: {[int(charges[i]) for i in sigma_indices]}")

    # Epsilon: FIRST k=0, charge=0 state above ground
    eps_idx = None
    for i in range(1, min(25, len(ratios))):
        k = int(momenta[i]); sk = k if k <= n_sites//2 else k - n_sites
        ch = int(charges[i])
        if sk == 0 and ch == 0 and i not in sigma_indices:
            eps_idx = i
            break

    if eps_idx is not None:
        print(f"  Epsilon: idx={eps_idx}, R={ratios[eps_idx]:.4f}, charge={int(charges[eps_idx])}")
        C_sse = float(np.abs(me_Z[eps_idx, best_sigma]) / me_0_sigma) if me_0_sigma > 1e-10 else 0.0
        print(f"  C_{{sigma,sigma*,epsilon}} = {C_sse:.4f}")
    else:
        C_sse = None
        print("  Epsilon: NOT FOUND")

    # ALL charge-0, k=0 states
    print(f"\n  === All charge-0, k=0 states (potential epsilon/descendants) ===")
    for i in range(min(25, n_me)):
        k = int(momenta[i]); sk = k if k <= n_sites//2 else k - n_sites
        ch = int(charges[i])
        if sk == 0 and ch == 0:
            me = np.abs(me_Z[i, best_sigma])
            C = me / me_0_sigma if me_0_sigma > 1e-10 else 0
            print(f"  i={i:>3}, R={ratios[i]:>8.4f}, |<i|Z|sigma>|={me:.6f}, C={C:.4f}")

    # Cross-check: sigma^k states should have charge k and <sigma^k|Z|sigma*>=0 for k≠0
    print(f"\n  === Z_q selection rule check ===")
    for i in range(min(20, n_me)):
        me = np.abs(me_Z[i, best_sigma])
        ch = int(charges[i])
        if me > 0.01:
            expected_charge = 0  # Z(charge 1) applied to sigma*(charge q-1) → charge 0
            match = "✓" if ch == 0 else "✗ VIOLATION"
            print(f"  i={i}, R={ratios[i]:.4f}, charge={ch}, |<i|Z|sigma>|={me:.4f} {match}")

    result = {
        "q": q, "n": n_sites, "dim": dim, "g_c": g_c,
        "me_0_Z_sigma": float(me_0_sigma),
        "C_sigma_sigma_epsilon": C_sse,
        "R_epsilon": float(ratios[eps_idx]) if eps_idx is not None else None,
        "epsilon_charge": int(charges[eps_idx]) if eps_idx is not None else None,
        "time_s": round(time.time() - t0, 2)
    }
    return result


# ===== Run all q values with charge-resolved epsilon ID =====
all_results = {}
configs = [
    (2, [6, 8, 10, 12], 0.250),
    (3, [6, 8, 10], 1/3),
    (4, [4, 6, 8], 0.392),
    (5, [4, 6, 8], 0.441),
    (7, [4, 6], 0.535),
    (10, [4, 6], 0.684),
]

for q, sizes, g_c in configs:
    print(f"\n{'='*70}")
    print(f"q={q} POTTS at g_c={g_c:.4f}")
    print(f"{'='*70}")
    all_results[f"q={q}"] = {}
    for n in sizes:
        ne = 15 if (q >= 10 and n >= 6) else 25
        r = extract_ope_corrected(n, q, g_c, n_eig=ne)
        if r:
            all_results[f"q={q}"][f"n={n}"] = r


# ===== FINAL COMPREHENSIVE SUMMARY =====
print("\n\n" + "=" * 70)
print("FINAL OPE SUMMARY: C_{sigma,sigma*,epsilon}(q)")
print("=" * 70)
print(f"  {'q':>3} {'n':>4} {'C_sse':>8} {'R_eps':>8}")
print(f"  {'-'*30}")

by_q = {}
for q_label, q_data in sorted(all_results.items()):
    q = int(q_label.split("=")[1])
    for n_label, res in sorted(q_data.items(), key=lambda x: -x[1]["n"]):
        C = res.get("C_sigma_sigma_epsilon")
        R = res.get("R_epsilon")
        if C is not None:
            print(f"  {q:>3} {res['n']:>4} {C:>8.4f} {R:>8.2f}" if R else f"  {q:>3} {res['n']:>4} {C:>8.4f}")
            if q not in by_q: by_q[q] = []
            by_q[q].append((res["n"], C))

# 1/N extrapolation
print(f"\n  === 1/N EXTRAPOLATION ===")
print(f"  {'q':>3} {'C(largest)':>12} {'C(∞) est':>12} {'sizes':>12}")
extrapolated = {}
for q in sorted(by_q.keys()):
    pairs = sorted(by_q[q], key=lambda x: x[0])
    n_max, C_max = pairs[-1]
    if len(pairs) >= 3:
        # Use last 3 points for quadratic extrapolation
        ns = [p[0] for p in pairs[-3:]]
        Cs = [p[1] for p in pairs[-3:]]
        inv_ns = [1/n for n in ns]
        # Linear fit in 1/N
        A = np.array([[1, x] for x in inv_ns])
        coeffs = np.linalg.lstsq(A, Cs, rcond=None)[0]
        C_inf = coeffs[0]
        method = f"({ns[0]},{ns[1]},{ns[2]})"
    elif len(pairs) >= 2:
        n1, C1 = pairs[-2]
        n2, C2 = pairs[-1]
        slope = (C2 - C1) / (1/n2 - 1/n1)
        C_inf = C2 - slope / n2
        method = f"({n1},{n2})"
    else:
        C_inf = C_max
        method = f"({n_max})"
    extrapolated[q] = C_inf
    print(f"  {q:>3} {C_max:>12.4f} {C_inf:>12.4f} {method:>12}")

# Check if C_sse ~ 1/q or another simple formula
print(f"\n  === FUNCTIONAL FORM ===")
print(f"  {'q':>3} {'C(∞)':>8} {'1/sqrt(q)':>10} {'1/q':>8} {'2/(q+1)':>8}")
for q in sorted(extrapolated.keys()):
    C = extrapolated[q]
    print(f"  {q:>3} {C:>8.4f} {1/np.sqrt(q):>10.4f} {1/q:>8.4f} {2/(q+1):>8.4f}")

with open("results/sprint_060e_ope_final.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print("\nResults saved.")
