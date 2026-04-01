"""Sprint 061b: Full CFT data for q=15,20,25,30 at n=4.

Extract: spectrum ratios, operator content, OPE C_sse, c/x_1 ratio,
harmonic ratios x(sigma^k)/x(sigma), descendant structure.

Using formula g_c values: g_c ≈ (1/5)√(q-1) + 1/20.
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import eigsh
import json, time

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
        ops = [eye_q]*n; ops[i] = XpXd
        r = csr_matrix(ops[0])
        for op in ops[1:]:
            r = sp_kron(r, csr_matrix(op), format='csr')
        H = H - g * r
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

def build_clock_op(n, q, power=1, site=0):
    omega = np.exp(2j * np.pi / q)
    Z = np.diag([omega**(power*s) for s in range(q)])
    eye_q = np.eye(q)
    ops = [eye_q]*n; ops[site] = Z
    r = csr_matrix(ops[0])
    for op in ops[1:]:
        r = sp_kron(r, csr_matrix(op), format='csr')
    return r

def analyze_cft(q, g_c, n_sites=4, n_eig=None):
    """Full CFT analysis at a single (q, g_c, n)."""
    dim = q**n_sites
    if n_eig is None:
        n_eig = min(40, dim - 2)
    n_eig = min(n_eig, dim - 2)

    print(f"\n{'='*60}")
    print(f"q={q}, g_c={g_c:.4f}, n={n_sites}, dim={dim}")
    print(f"{'='*60}")

    t0 = time.time()

    # Diagonalize
    H = potts_hamiltonian_periodic(n_sites, q, g_c)
    T = build_translation_operator(n_sites, q)
    evals, evecs = eigsh(H, k=n_eig, which='SA')
    order = np.argsort(evals)
    evals = evals[order]; evecs = evecs[:, order]

    # Gaps and ratios
    gaps = evals - evals[0]
    delta1 = gaps[1]
    ratios = gaps / delta1

    # Momenta
    n_states = len(evals)
    momenta = np.zeros(n_states, dtype=int)
    tol = 1e-8 * abs(evals[-1] - evals[0])
    i = 0
    while i < n_states:
        group = [i]; j = i + 1
        while j < n_states and abs(evals[j] - evals[i]) < tol:
            group.append(j); j += 1
        sub = evecs[:, group]
        T_sub = (sub.conj().T) @ (T @ sub)
        t_evals = np.linalg.eigvals(T_sub)
        for gi, g_idx in enumerate(group):
            phase = np.angle(t_evals[gi])
            k = phase / (2 * np.pi / n_sites)
            momenta[g_idx] = int(round(k)) % n_sites
        i = j

    t_diag = time.time() - t0

    # Print spectrum
    print(f"\n  Spectrum (first 30 levels, {t_diag:.1f}s):")
    print(f"  {'i':>3} {'R':>8} {'k':>4} {'degen':>6}")
    degen_count = {}
    for i in range(min(30, len(ratios))):
        k = int(momenta[i]); s = k if k <= n_sites//2 else k - n_sites
        # Count degeneracy
        r_key = round(ratios[i], 3)
        degen_count[r_key] = degen_count.get(r_key, 0) + 1
        print(f"  {i:>3} {ratios[i]:>8.4f} {s:>4}")

    # === Identify operators ===
    # Sigma: lowest k=0 states above ground (degenerate pair for q>=3)
    sigma_indices = [i for i in range(1, min(2*q, len(ratios)))
                     if ratios[i] < 1.5 and int(momenta[i]) == 0]
    if not sigma_indices:
        sigma_indices = [1]

    # Harmonic fields: sigma^k states
    # For each power, the sigma^k field has R(sigma^k)/R(sigma) ~ k^alpha
    harmonic_ratios = {}

    # === OPE via clock operator ===
    Z_local = build_clock_op(n_sites, q, power=1, site=0)
    n_me = min(30, n_eig)
    Z_evecs = Z_local @ evecs[:, :n_me]
    me_Z = np.zeros((n_me, n_me), dtype=complex)
    for a in range(n_me):
        me_Z[a, :] = evecs[:, a].conj() @ Z_evecs

    # Find best sigma via matrix element
    best_sigma = max(sigma_indices, key=lambda i: abs(me_Z[0, i]) if i < n_me else 0)
    me_0_sigma = abs(me_Z[0, best_sigma])

    print(f"\n  Sigma: idx={best_sigma}, R={ratios[best_sigma]:.4f}, |<0|Z|sigma>|={me_0_sigma:.6f}")

    # Find epsilon: k=0 singlet above all sigma harmonics
    # Strategy: look for k=0 states with large <eps|Z|sigma>/<0|Z|sigma>
    # and R > max sigma harmonic R
    eps_candidates = []
    for i in range(1, n_me):
        k = int(momenta[i]); s = k if k <= n_sites//2 else k - n_sites
        if s == 0 and abs(me_Z[i, best_sigma]) > 1e-8:
            C = abs(me_Z[i, best_sigma]) / me_0_sigma
            if C > 0.05 and ratios[i] > 2.0:  # Must be well above sigma
                eps_candidates.append((i, ratios[i], C))

    # Sort by ratio — epsilon should be at characteristic R
    eps_candidates.sort(key=lambda x: x[1])

    # The true epsilon is the FIRST k=0 state that couples via Z to sigma
    # and is NOT a sigma harmonic (not degenerate with a partner)
    eps_idx = None
    C_sse = None
    if eps_candidates:
        # Pick the one with highest OPE coupling
        best_eps = max(eps_candidates, key=lambda x: x[2])
        eps_idx = best_eps[0]
        C_sse = best_eps[2]
        print(f"  Epsilon: idx={eps_idx}, R={ratios[eps_idx]:.4f}, C_sse={C_sse:.4f}")

    # === Harmonic ratios via Z^power ===
    max_power = min(q//2, 6)
    for power in range(2, max_power + 1):
        Zp = build_clock_op(n_sites, q, power=power, site=0)
        Zp_evecs = Zp @ evecs[:, :n_me]
        me_Zp = np.zeros(n_me, dtype=complex)
        for a in range(n_me):
            me_Zp[a] = np.dot(evecs[:, 0].conj(), Zp_evecs[:, a])

        # Find sigma^power: k=0 state with largest <0|Z^p|state>
        candidates = [(i, ratios[i], abs(me_Zp[i])) for i in range(1, n_me)
                      if int(momenta[i]) == 0 and abs(me_Zp[i]) > 1e-8]
        if candidates:
            best = max(candidates, key=lambda x: x[2])
            R_sp = best[1]
            harmonic_ratios[power] = R_sp / ratios[best_sigma]
            print(f"  sigma^{power}: idx={best[0]}, R={R_sp:.4f}, R/R_sigma={harmonic_ratios[power]:.4f}, "
                  f"k^2={power**2}")

    # === c/x1 from Casimir energy ===
    E0_per_site = evals[0] / n_sites
    delta1_N = delta1 * n_sites
    # At CFT: E0/N = e_inf - pi*v*c/(6*N^2), and delta = 2*pi*v*x1/N
    # So c/x1 = -3*(E0/N - e_inf) * N / delta1 = c/(6*x1) * (2*pi)
    # Actually: ratio method from (E0(N), Delta(N)) at two sizes is better.
    # With one size, we can't separate e_inf from c. Report delta1*N only.

    # === Descendant check ===
    # Look for states with k=±1 near R = R_sigma + 1/x_sigma
    desc_candidates = [(i, ratios[i]) for i in range(1, n_me)
                       if abs(int(momenta[i]) if int(momenta[i]) <= n_sites//2
                              else int(momenta[i]) - n_sites) == 1
                       and ratios[i] > ratios[best_sigma] + 2]
    desc_R = desc_candidates[0][1] if desc_candidates else None
    desc_gap = (desc_R - ratios[best_sigma]) if desc_R else None

    print(f"\n  Descendant L_{{-1}}sigma: R={desc_R:.4f}, gap from sigma={desc_gap:.4f}" if desc_R else "")
    print(f"  Delta_1 * N = {delta1_N:.6f}")

    # Count sigma-type fields
    # At k=0, count how many states have R < R_epsilon and are degenerate pairs
    n_sigma_fields = 0
    if eps_idx:
        R_eps = ratios[eps_idx]
        for i in range(1, n_me):
            if ratios[i] < R_eps * 0.95 and int(momenta[i]) == 0:
                n_sigma_fields += 1

    result = {
        "q": q, "n": n_sites, "dim": dim, "g_c": g_c,
        "E0": float(evals[0]), "E0_per_site": float(E0_per_site),
        "delta1": float(delta1), "delta1_N": float(delta1_N),
        "R_sigma": float(ratios[best_sigma]),
        "R_epsilon": float(ratios[eps_idx]) if eps_idx else None,
        "C_sse": float(C_sse) if C_sse else None,
        "harmonic_ratios": {str(k): float(v) for k, v in harmonic_ratios.items()},
        "free_boson_k2": {str(k): k**2 for k in harmonic_ratios.keys()},
        "n_sigma_fields": n_sigma_fields,
        "expected_sigma_fields": q - 1,
        "desc_gap": float(desc_gap) if desc_gap else None,
        "time_s": round(time.time() - t0, 2),
    }

    return result


# Formula for g_c
def g_c_formula(q):
    return 0.2 * np.sqrt(q - 1) + 0.05

# ===== Run for q = 15, 20, 25, 30 =====
all_results = {}

for q in [15, 20, 25, 30]:
    gc = g_c_formula(q)
    dim4 = q**4
    if dim4 > 1e6:
        n_eig = 25  # Limit eigenvalues for large dim
    else:
        n_eig = 40
    result = analyze_cft(q, gc, n_sites=4, n_eig=n_eig)
    all_results[f"q={q}"] = result

# ===== Summary table =====
print("\n\n" + "=" * 70)
print("COMPREHENSIVE SUMMARY — Large-q CFT Data (n=4)")
print("=" * 70)

print(f"\n  {'q':>3} {'g_c':>6} {'D1N':>8} {'R_eps':>7} {'C_sse':>7} {'#sig':>5} {'desc_gap':>9}")
print(f"  {'-'*55}")
for key, r in sorted(all_results.items(), key=lambda x: x[1]["q"]):
    q = r["q"]
    R_eps = f"{r['R_epsilon']:.3f}" if r['R_epsilon'] else "N/A"
    C_sse = f"{r['C_sse']:.4f}" if r['C_sse'] else "N/A"
    dg = f"{r['desc_gap']:.3f}" if r['desc_gap'] else "N/A"
    print(f"  {q:>3} {r['g_c']:>6.3f} {r['delta1_N']:>8.4f} {R_eps:>7} {C_sse:>7} {r['n_sigma_fields']:>5} {dg:>9}")

# Harmonic ratios vs k^2
print(f"\n  Harmonic ratios R(sigma^k)/R(sigma) vs free boson k²:")
print(f"  {'q':>3}", end="")
for k in range(2, 7):
    print(f"  {'k='+str(k):>8}", end="")
print()
for key, r in sorted(all_results.items(), key=lambda x: x[1]["q"]):
    q = r["q"]
    hr = r["harmonic_ratios"]
    print(f"  {q:>3}", end="")
    for k in range(2, 7):
        val = hr.get(str(k), None)
        if val is not None:
            print(f"  {val:>8.3f}", end="")
        else:
            print(f"  {'---':>8}", end="")
    print()
print(f"  {'k²':>3}", end="")
for k in range(2, 7):
    print(f"  {k**2:>8.1f}", end="")
print()

with open("results/sprint_061b_cft_large_q.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print("\n\nResults saved to results/sprint_061b_cft_large_q.json")
