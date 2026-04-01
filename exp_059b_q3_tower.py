"""Sprint 059b: Conformal tower for q=3 Potts at g_c=1/3 with momentum resolution.

q=3 Potts CFT (c=4/5) known primaries:
  σ:  x = 2/15 ≈ 0.1333  (Z₃ spin field, doubly degenerate: σ, σ†)
  ε:  x = 4/5 = 0.8       (energy field)
  μ:  x = 4/3 ≈ 1.333     (disorder field, doubly degenerate)

Descendants: L₋₁σ at x = 2/15 + 1 = 17/15 ≈ 1.133, spin ±1
Gap ratio: R(L₋₁σ)/R(σ) = (17/15)/(2/15) = 17/2 = 8.5

So R(desc) - R(primary) = 8.5 - 1.0 = 7.5 = 1/x₁ = 15/2
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
    """H = -sum delta(s_i,s_{i+1}) - g*sum(X + X^dag)"""
    dim = q**n
    eye_q = np.eye(q)

    # X operator: |s> -> |s+1 mod q>
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = X + X.T

    # delta projector for 2 sites
    delta_2 = np.zeros((q**2, q**2))
    for s in range(q):
        delta_2[s*q+s, s*q+s] = 1.0

    H = csr_matrix((dim, dim))

    # Interaction bonds
    for i in range(n):
        j = (i + 1) % n
        if j == i + 1:  # adjacent
            left_dim = q**i if i > 0 else 1
            right_dim = q**(n-j-1) if n-j-1 > 0 else 1
            op = csr_matrix(delta_2)
            if i > 0:
                op = sp_kron(csr_matrix(np.eye(left_dim)), op, format='csr')
            if n-j-1 > 0:
                op = sp_kron(op, csr_matrix(np.eye(right_dim)), format='csr')
            H = H - op
        else:  # wrap: i=n-1, j=0
            for s in range(q):
                proj = np.zeros((q, q)); proj[s, s] = 1.0
                ops = [eye_q] * n; ops[0] = proj; ops[n-1] = proj
                H = H - tensor_op(ops)

    # Transverse field
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

def get_momentum_robust(evals, evecs, T, n):
    n_states = len(evals)
    momenta = np.zeros(n_states, dtype=int)
    tol = 1e-8 * abs(evals[-1] - evals[0])
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

# ========== q=3 Potts ==========
print("=" * 70)
print("CONFORMAL TOWER ANALYSIS: q=3 POTTS at g_c=1/3")
print("=" * 70)

q = 3
g_c = 1.0 / 3.0
results = {}

for n in [6, 8, 10]:
    t0 = time.time()
    dim = q**n
    print(f"\nBuilding n={n} (dim={dim})...", end=" ", flush=True)

    H = potts_hamiltonian_periodic(n, q, g_c)
    T = build_translation_operator(n, q)

    n_eig = min(30, dim - 2)
    evals, evecs = eigsh(H, k=n_eig, which='SA')
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]
    dt = time.time() - t0
    print(f"{dt:.1f}s")

    momenta = get_momentum_robust(evals, evecs, T, n)

    gaps = evals - evals[0]
    delta1 = gaps[1]  # Should be doubly degenerate (σ, σ†)
    ratios = gaps / delta1

    # x₁ = 2/15 for q=3
    x1_exact = 2.0/15.0

    print(f"  E0={evals[0]:.6f}, Δ₁={delta1:.6f}, Δ₁·N={delta1*n:.4f}")
    print(f"  {'i':>3} {'R':>8} {'k':>3} {'spin':>5} {'x=R·x₁':>8} {'degen':>5} {'CFT tower':>30}")
    print(f"  {'-'*68}")

    # Count degeneracies
    i = 0
    level_info = []
    while i < min(25, len(ratios)):
        degen = 1
        j = i + 1
        while j < len(ratios) and abs(ratios[j] - ratios[i]) < 0.02:
            degen += 1
            j += 1

        k = int(momenta[i])
        spin = k if k <= n//2 else k - n
        x_val = ratios[i] * x1_exact
        R = ratios[i]

        tower = ""
        if abs(R) < 0.01: tower = "Identity (0)"
        elif abs(R - 1.0) < 0.1: tower = f"σ,σ† primary (2/15)"
        elif abs(R - 6.0) < 0.4: tower = f"ε primary (4/5)"
        elif abs(R - 8.5) < 0.5: tower = f"L₋₁σ desc (17/15)"
        elif abs(R - 10.0) < 0.5: tower = f"μ,μ† (4/3)"
        elif abs(R - 14.5) < 0.6: tower = f"σ level-2"
        elif abs(R - 15.0) < 0.6: tower = f"L₋₁ε desc"

        # Get all spins in this group
        spins = []
        for ii in range(i, i + degen):
            if ii < len(momenta):
                kk = int(momenta[ii])
                ss = kk if kk <= n//2 else kk - n
                spins.append(ss)

        spin_str = ",".join(f"{s:+d}" for s in sorted(set(spins)))
        print(f"  {i:>3} {R:>8.3f} {spin_str:>7} {x_val:>8.4f} {degen:>5} {tower:>30}")

        level_info.append({
            "i": int(i), "R": round(float(R), 6),
            "x": round(float(x_val), 6), "degen": degen,
            "spins": sorted(set(int(s) for s in spins))
        })
        i = j

    results[f"n={n}"] = {
        "n": n, "dim": dim, "time_s": round(dt, 2),
        "E0": round(float(evals[0]), 8),
        "delta1": round(float(delta1), 8),
        "delta1_N": round(float(delta1 * n), 6),
        "levels": level_info
    }

# Summary
print("\n\n" + "=" * 70)
print("q=3 TOWER STRUCTURE SUMMARY (n=10)")
print("=" * 70)
print("CFT predictions (c=4/5, x₁=2/15):")
print("  σ,σ† at R=1.0, spin=0, degen=2")
print("  ε at R=6.0, spin=0, degen=1")
print("  L₋₁σ at R=8.5, spin=±1, degen=4 (2 primaries × 2 chiralities)")
print("  μ,μ† at R=10.0, spin=0, degen=2")
print("  Descendant gap: R(desc)-R(σ) = 7.5 = 1/x₁ = 15/2")

data = results["n=10"]["levels"]
sigma_R = [l for l in data if 0.9 < l["R"] < 1.1]
desc_R = [l for l in data if 8.0 < l["R"] < 9.0]
eps_R = [l for l in data if 5.5 < l["R"] < 6.5]
mu_R = [l for l in data if 9.5 < l["R"] < 10.5]

print(f"\nMeasured:")
if sigma_R: print(f"  σ: R={sigma_R[0]['R']:.3f}, degen={sigma_R[0]['degen']}, spins={sigma_R[0]['spins']}")
if eps_R: print(f"  ε: R={eps_R[0]['R']:.3f}, degen={eps_R[0]['degen']}, spins={eps_R[0]['spins']}")
if desc_R: print(f"  L₋₁σ: R={desc_R[0]['R']:.3f}, degen={desc_R[0]['degen']}, spins={desc_R[0]['spins']}")
if mu_R: print(f"  μ: R={mu_R[0]['R']:.3f}, degen={mu_R[0]['degen']}, spins={mu_R[0]['spins']}")
if sigma_R and desc_R:
    gap = desc_R[0]["R"] - sigma_R[0]["R"]
    print(f"\n  Descendant gap: {gap:.3f} (predicted: 7.5 = 1/x₁)")

with open("results/sprint_059b_q3_tower.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nResults saved.")
