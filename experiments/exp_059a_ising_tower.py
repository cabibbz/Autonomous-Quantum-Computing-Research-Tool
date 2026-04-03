"""Sprint 059a: Conformal tower analysis for q=2 Ising (TFIM) with momentum resolution.

At criticality on periodic chain of length N, CFT predicts conformal towers:
  E_n - E_0 = (2*pi*v_s/N) * x_n

For Ising CFT (c=1/2), primaries and their towers:
  Identity tower: x = 0, 2, 2, ...      spin s=0, ±2, ...
  Sigma tower:    x = 1/8, 9/8, 9/8, ... spin s=0, ±1, ...
  Epsilon tower:  x = 1, 2, 2, ...        spin s=0, ±1, ...

We resolve momentum k from the translation operator T, handling degeneracies
by diagonalizing T within degenerate eigenspaces of H.
"""
import numpy as np
from scipy.sparse import csr_matrix, eye as sp_eye, kron as sp_kron
from scipy.sparse.linalg import eigsh
import json, time

def tensor_op(ops_list):
    result = csr_matrix(ops_list[0])
    for op in ops_list[1:]:
        result = sp_kron(result, csr_matrix(op), format='csr')
    return result

def tfim_hamiltonian_periodic(n, g):
    dim = 2**n
    sz = np.array([[1,0],[0,-1]], dtype=float)
    sx = np.array([[0,1],[1,0]], dtype=float)
    eye2 = np.eye(2)
    H = csr_matrix((dim, dim))
    for i in range(n):
        j = (i + 1) % n
        ops = [eye2] * n; ops[i] = sz; ops[j] = sz
        H = H - tensor_op(ops)
    for i in range(n):
        ops = [eye2] * n; ops[i] = sx
        H = H - g * tensor_op(ops)
    return H

def build_translation_operator(n, q):
    """T: |s_1 s_2 ... s_n> -> |s_2 s_3 ... s_n s_1>."""
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

def get_momentum_robust(evals, evecs, T, n):
    """Extract momentum k for each eigenstate, handling degeneracies.
    Within degenerate subspaces, diagonalize T to get definite momentum."""
    n_states = len(evals)
    momenta = np.zeros(n_states, dtype=int)
    tol = 1e-8 * abs(evals[-1] - evals[0])  # relative tolerance
    processed = [False] * n_states

    i = 0
    while i < n_states:
        if processed[i]:
            i += 1
            continue
        # Find degenerate group
        group = [i]
        j = i + 1
        while j < n_states and abs(evals[j] - evals[i]) < tol:
            group.append(j)
            j += 1

        # Diagonalize T in this subspace
        sub = evecs[:, group]
        T_sub = (sub.conj().T) @ (T @ sub)
        t_evals = np.linalg.eigvals(T_sub)

        # Sort T eigenvalues and assign momenta
        for idx_g, g_idx in enumerate(group):
            phase = np.angle(t_evals[idx_g])
            k = phase / (2 * np.pi / n)
            momenta[g_idx] = int(round(k)) % n
            processed[g_idx] = True

        i = j

    return momenta

# ========== q=2 Ising ==========
print("=" * 70)
print("CONFORMAL TOWER ANALYSIS: q=2 ISING (TFIM) at g_c=1.0")
print("=" * 70)

results = {}

for n in [8, 10, 12]:
    t0 = time.time()
    H = tfim_hamiltonian_periodic(n, g=1.0)
    T = build_translation_operator(n, q=2)

    n_eig = min(30, 2**n - 2)
    evals, evecs = eigsh(H, k=n_eig, which='SA')
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]
    dt = time.time() - t0

    momenta = get_momentum_robust(evals, evecs, T, n)

    gaps = evals - evals[0]
    delta1 = gaps[1]
    ratios = gaps / delta1

    print(f"\nn={n} (dim={2**n}), time={dt:.1f}s")
    print(f"  E0 = {evals[0]:.6f}, Δ₁ = {delta1:.6f}, Δ₁·N = {delta1*n:.4f}")
    print(f"  {'i':>3} {'R':>8} {'k':>3} {'spin':>5} {'x=R·x₁':>8} {'CFT tower':>25}")
    print(f"  {'-'*58}")

    # CFT Ising: x₁ = 1/8. Descendants at x₁+1 = 9/8 → R = 9
    # ε at x = 1 → R = 8. ε descendants at x=2 → R = 16
    for i in range(min(25, len(ratios))):
        k = int(momenta[i])
        spin = k if k <= n//2 else k - n
        x_val = ratios[i] * 0.125  # x₁ = 1/8

        # Assign tower
        tower = ""
        R = ratios[i]
        if abs(R) < 0.01: tower = "I primary (0)"
        elif abs(R - 1) < 0.1: tower = "σ primary (1/8)"
        elif abs(R - 8) < 0.3: tower = "ε primary (1)"
        elif abs(R - 9) < 0.4: tower = "L₋₁σ (9/8)"
        elif abs(R - 16) < 0.6: tower = "level-2"
        elif abs(R - 17) < 0.6: tower = "L₋₂σ or L₋₁ε"

        print(f"  {i:>3} {R:>8.3f} {k:>3} {spin:>+5d} {x_val:>8.4f} {tower:>25}")

    results[f"n={n}"] = {
        "n": n, "dim": 2**n, "time_s": round(dt, 2),
        "E0": round(float(evals[0]), 8),
        "delta1": round(float(delta1), 8),
        "delta1_N": round(float(delta1 * n), 6),
        "levels": [
            {"i": int(i), "R": round(float(ratios[i]), 6),
             "k": int(momenta[i]),
             "spin": int(momenta[i] if momenta[i] <= n//2 else momenta[i] - n),
             "x": round(float(ratios[i] * 0.125), 6)}
            for i in range(min(25, len(ratios)))
        ]
    }

# Tower structure summary at n=12
print("\n\n" + "=" * 70)
print("TOWER STRUCTURE SUMMARY (n=12)")
print("=" * 70)
data = results["n=12"]["levels"]

print("\n--- σ tower (primary at R=1, x=1/8) ---")
print("CFT prediction: primary at k=0, descendants at k=±1 (R=9), k=0,±2 (R=16)")
for l in data:
    if 0.9 < l["R"] < 1.1 or 8.5 < l["R"] < 9.5 or 15.0 < l["R"] < 17.5:
        print(f"  R={l['R']:.3f}, spin={l['spin']:+d}, x={l['x']:.4f}")

print("\n--- ε tower (primary at R=8, x=1) ---")
print("CFT prediction: primary at k=0, descendants at k=±1 (R=16)")
for l in data:
    if 7.5 < l["R"] < 8.5 or 15.5 < l["R"] < 16.5:
        print(f"  R={l['R']:.3f}, spin={l['spin']:+d}, x={l['x']:.4f}")

# Check descendant gap
print("\n--- Descendant gap check ---")
sigma_R = [l["R"] for l in data if 0.9 < l["R"] < 1.1]
desc_R = [l["R"] for l in data if 8.5 < l["R"] < 9.5]
if sigma_R and desc_R:
    gap = desc_R[0] - sigma_R[0]
    print(f"R(L₋₁σ) - R(σ) = {gap:.3f} (CFT prediction: 8.0 = 1/x₁)")

eps_R = [l["R"] for l in data if 7.5 < l["R"] < 8.5]
if eps_R:
    print(f"R(ε) = {eps_R[0]:.3f} (CFT prediction: 8.0)")

with open("results/sprint_059a_ising_tower.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nResults saved.")
