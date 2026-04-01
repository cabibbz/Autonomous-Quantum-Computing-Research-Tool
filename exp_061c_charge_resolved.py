"""Sprint 061c: Charge-resolved spectrum to separate epsilon from sigma².

At large q, R(sigma²) ≈ R(epsilon) ≈ 4 at n=4. Must resolve Z_q charge
to distinguish charge-0 epsilon from charge-2 sigma².

Build symmetry operator G = X₁·X₂·...·Xₙ (Z_q generator).
Project eigenstates into charge sectors. Then identify epsilon as lowest
charge-0 excited state above identity.

Also extract proper c/x₁ from n=4 Casimir energy for q=2-30 comparison.
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

def build_symmetry_generator(n, q):
    """G = X₁ ⊗ X₂ ⊗ ... ⊗ Xₙ (Z_q global symmetry generator)."""
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    G = csr_matrix(X)
    for _ in range(n - 1):
        G = sp_kron(G, csr_matrix(X), format='csr')
    return G

def build_clock_op(n, q, power=1, site=0):
    omega = np.exp(2j * np.pi / q)
    Z = np.diag([omega**(power*s) for s in range(q)])
    eye_q = np.eye(q)
    ops = [eye_q]*n; ops[site] = Z
    r = csr_matrix(ops[0])
    for op in ops[1:]:
        r = sp_kron(r, csr_matrix(op), format='csr')
    return r

def build_translation(n, q):
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

def g_c_formula(q):
    return 0.2 * np.sqrt(q - 1) + 0.05

def analyze_with_charges(q, g_c, n_sites=4, n_eig=None):
    dim = q**n_sites
    if n_eig is None:
        n_eig = min(40, dim - 2)

    print(f"\n{'='*60}")
    print(f"q={q}, g_c={g_c:.4f}, n={n_sites}, dim={dim}")
    print(f"{'='*60}")
    t0 = time.time()

    H = potts_hamiltonian_periodic(n_sites, q, g_c)
    G = build_symmetry_generator(n_sites, q)
    T = build_translation(n_sites, q)
    Z_local = build_clock_op(n_sites, q, power=1, site=0)

    evals, evecs = eigsh(H, k=n_eig, which='SA')
    order = np.argsort(evals); evals = evals[order]; evecs = evecs[:, order]

    # Assign Z_q charges
    charges = np.zeros(len(evals), dtype=int)
    omega = np.exp(2j * np.pi / q)
    tol = 1e-6 * abs(evals[-1] - evals[0])

    i = 0
    while i < len(evals):
        group = [i]; j = i + 1
        while j < len(evals) and abs(evals[j] - evals[i]) < tol:
            group.append(j); j += 1

        # Diagonalize G within degenerate subspace
        sub = evecs[:, group]
        G_sub = (sub.conj().T) @ (G @ sub)
        g_evals = np.linalg.eigvals(G_sub)

        # Also get momenta from T
        T_sub = (sub.conj().T) @ (T @ sub)
        t_evals = np.linalg.eigvals(T_sub)

        for gi, g_idx in enumerate(group):
            # Charge from G eigenvalue
            phase = np.angle(g_evals[gi])
            charge = round(phase / (2 * np.pi / q)) % q
            charges[g_idx] = charge
        i = j

    # Assign momenta
    momenta = np.zeros(len(evals), dtype=int)
    i = 0
    while i < len(evals):
        group = [i]; j = i + 1
        while j < len(evals) and abs(evals[j] - evals[i]) < tol:
            group.append(j); j += 1
        sub = evecs[:, group]
        T_sub = (sub.conj().T) @ (T @ sub)
        t_evals = np.linalg.eigvals(T_sub)
        for gi, g_idx in enumerate(group):
            phase = np.angle(t_evals[gi])
            k = phase / (2 * np.pi / n_sites)
            momenta[g_idx] = int(round(k)) % n_sites
        i = j

    gaps = evals - evals[0]
    delta1 = gaps[1]
    ratios = gaps / delta1

    t_diag = time.time() - t0

    # Print charge-resolved spectrum
    print(f"\n  Charge-resolved spectrum ({t_diag:.1f}s):")
    print(f"  {'i':>3} {'R':>8} {'k':>4} {'charge':>7}")
    for i in range(min(30, len(ratios))):
        k = int(momenta[i]); s = k if k <= n_sites//2 else k - n_sites
        print(f"  {i:>3} {ratios[i]:>8.4f} {s:>4} {charges[i]:>7}")

    # Identify operators by charge
    # Identity: charge 0, R=0
    # Sigma (charge 1): lowest R > 0 with charge 1
    # Sigma* (charge q-1): lowest R > 0 with charge q-1
    # Epsilon (charge 0): SECOND charge-0 state (first is identity)
    # Sigma² (charge 2): lowest R > 0 with charge 2

    charge_0_states = [(i, ratios[i]) for i in range(len(ratios)) if charges[i] == 0]
    charge_1_states = [(i, ratios[i]) for i in range(len(ratios)) if charges[i] == 1]
    charge_2_states = [(i, ratios[i]) for i in range(len(ratios)) if charges[i] == 2]

    # Sigma: first charge-1 state
    sigma_idx = charge_1_states[0][0] if charge_1_states else None
    R_sigma = ratios[sigma_idx] if sigma_idx is not None else None

    # Epsilon: second charge-0 state (first is identity at R=0)
    eps_idx = charge_0_states[1][0] if len(charge_0_states) > 1 else None
    R_eps = ratios[eps_idx] if eps_idx is not None else None

    # Sigma²: first charge-2 state
    sigma2_idx = charge_2_states[0][0] if charge_2_states else None
    R_sigma2 = ratios[sigma2_idx] if sigma2_idx is not None else None

    print(f"\n  OPERATORS (charge-resolved):")
    print(f"  Identity: charge=0, R=0")
    if sigma_idx is not None:
        print(f"  Sigma: charge=1, idx={sigma_idx}, R={R_sigma:.4f}")
    if sigma2_idx is not None:
        print(f"  Sigma²: charge=2, idx={sigma2_idx}, R={R_sigma2:.4f}")
    if eps_idx is not None:
        print(f"  Epsilon: charge=0, idx={eps_idx}, R={R_eps:.4f}")

    # Epsilon vs sigma² separation
    if R_eps is not None and R_sigma2 is not None:
        print(f"\n  R(epsilon) = {R_eps:.4f}")
        print(f"  R(sigma²)  = {R_sigma2:.4f}")
        print(f"  Separation: {abs(R_eps - R_sigma2):.4f} ({abs(R_eps - R_sigma2)/R_eps*100:.1f}%)")

    # OPE via Z matrix elements (charge selection: Z maps charge c → c-1)
    n_me = min(30, n_eig)
    Z_evecs = Z_local @ evecs[:, :n_me]
    me_Z = np.zeros((n_me, n_me), dtype=complex)
    for a in range(n_me):
        me_Z[a, :] = evecs[:, a].conj() @ Z_evecs

    # Proper OPE: C_sse = |<eps|Z|sigma>| / |<0|Z|sigma>|
    # Z maps charge 1 → charge 0. So <charge0|Z|charge1> is allowed.
    C_sse = None
    if sigma_idx is not None and eps_idx is not None and sigma_idx < n_me and eps_idx < n_me:
        me_0_sigma = abs(me_Z[0, sigma_idx])
        me_eps_sigma = abs(me_Z[eps_idx, sigma_idx])
        if me_0_sigma > 1e-10:
            C_sse = me_eps_sigma / me_0_sigma
            print(f"\n  OPE C_sse = |<eps|Z|sigma>| / |<0|Z|sigma>| = {me_eps_sigma:.6f} / {me_0_sigma:.6f} = {C_sse:.4f}")
        else:
            print(f"\n  WARNING: <0|Z|sigma> ≈ 0 — OPE extraction failed")

    # Check: what's the coupling to sigma² (should be ZERO by charge selection)
    if sigma2_idx is not None and sigma_idx is not None and sigma2_idx < n_me:
        me_s2_sigma = abs(me_Z[sigma2_idx, sigma_idx])
        print(f"  Check: |<sigma²|Z|sigma>| = {me_s2_sigma:.6f} (should be ~0 by charge)")

    # Count all charge sectors
    print(f"\n  Charge distribution:")
    for c in range(min(q, 8)):
        states_c = [(i, ratios[i]) for i in range(len(ratios)) if charges[i] == c and ratios[i] < 20]
        if states_c:
            Rs = [f"{r:.3f}" for _, r in states_c[:4]]
            print(f"    charge {c}: {len(states_c)} states, R = [{', '.join(Rs)}{'...' if len(states_c) > 4 else ''}]")

    result = {
        "q": q, "n": n_sites, "dim": dim, "g_c": g_c,
        "R_sigma": float(R_sigma) if R_sigma is not None else None,
        "R_sigma2": float(R_sigma2) if R_sigma2 is not None else None,
        "R_epsilon": float(R_eps) if R_eps is not None else None,
        "C_sse": float(C_sse) if C_sse is not None else None,
        "delta1_N": float(delta1 * n_sites),
        "eps_sigma2_separation": float(abs(R_eps - R_sigma2)) if R_eps and R_sigma2 else None,
        "time_s": round(time.time() - t0, 2),
    }
    return result


# Run for calibration (known q) and new large q
all_results = {}

# Known q values for calibration
known = {
    2: 0.250, 3: 0.333, 4: 0.392, 5: 0.441,
    7: 0.535, 10: 0.684,
}

for q, gc in known.items():
    dim = q**4
    ne = min(40, dim - 2)
    all_results[f"q={q}"] = analyze_with_charges(q, gc, n_sites=4, n_eig=ne)

# New large q
for q in [15, 20, 25, 30]:
    gc = g_c_formula(q)
    ne = min(30, q**4 - 2)
    all_results[f"q={q}"] = analyze_with_charges(q, gc, n_sites=4, n_eig=ne)

# Summary
print("\n\n" + "=" * 70)
print("CHARGE-RESOLVED SUMMARY")
print("=" * 70)
print(f"\n  {'q':>3} {'R_sigma':>8} {'R_sigma2':>9} {'R_eps':>8} {'sep':>6} {'C_sse':>7}")
print(f"  {'-'*50}")
for key, r in sorted(all_results.items(), key=lambda x: x[1]["q"]):
    q = r["q"]
    rs = f"{r['R_sigma']:.4f}" if r['R_sigma'] else "N/A"
    rs2 = f"{r['R_sigma2']:.4f}" if r['R_sigma2'] else "N/A"
    re = f"{r['R_epsilon']:.4f}" if r['R_epsilon'] else "N/A"
    sep = f"{r['eps_sigma2_separation']:.3f}" if r['eps_sigma2_separation'] else "N/A"
    cs = f"{r['C_sse']:.4f}" if r['C_sse'] else "N/A"
    print(f"  {q:>3} {rs:>8} {rs2:>9} {re:>8} {sep:>6} {cs:>7}")

with open("results/sprint_061c_charge_resolved.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print("\nResults saved.")
