"""Sprint 059c: Conformal towers for q=4 and q=5 Potts — novel regime.

q=4 (c≈1.0, marginal Ashkin-Teller):
  Spin fields: 3 = 1 pair (σ,σ³) + 1 self-conjugate (σ²)
  x₁ ≈ 0.117, so 1/x₁ ≈ 8.5
  Predicted descendant gap: R(L₋₁σ) - R(σ) = 1/x₁ ≈ 8.5

q=5 (c≈1.10, novel CFT):
  Spin fields: 4 = 2 pairs (σ,σ⁴) and (σ²,σ³)
  x₁ ≈ 0.101, so 1/x₁ ≈ 9.9
  Predicted descendant gap ≈ 9.9

KEY QUESTION: Do descendants appear with correct momentum and degeneracy,
confirming genuine CFT structure in the novel q>4 regime?
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

def analyze_tower(n_sites, q, g_c, x1_approx, results_dict, max_eig=30):
    """Run tower analysis for given q at g_c."""
    t0 = time.time()
    dim = q**n_sites
    print(f"\n  n={n_sites} (dim={dim})...", end=" ", flush=True)

    H = potts_hamiltonian_periodic(n_sites, q, g_c)
    T = build_translation_operator(n_sites, q)

    n_eig = min(max_eig, dim - 2)
    evals, evecs = eigsh(H, k=n_eig, which='SA')
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]
    dt = time.time() - t0
    print(f"{dt:.1f}s")

    momenta = get_momentum_robust(evals, evecs, T, n_sites)

    gaps = evals - evals[0]
    delta1 = gaps[1]
    ratios = gaps / delta1

    print(f"  E0={evals[0]:.6f}, Δ₁={delta1:.6f}, Δ₁·N={delta1*n_sites:.4f}")
    print(f"  {'i':>3} {'R':>8} {'spins':>10} {'x=R·x₁':>8} {'deg':>4} {'note':>25}")
    print(f"  {'-'*65}")

    # Group by degeneracy and extract
    level_info = []
    i = 0
    while i < min(n_eig, len(ratios)):
        degen = 1
        j = i + 1
        while j < len(ratios) and abs(ratios[j] - ratios[i]) < 0.03:
            degen += 1
            j += 1

        R = ratios[i]
        x_val = R * x1_approx
        spins = set()
        for ii in range(i, min(i+degen, len(momenta))):
            kk = int(momenta[ii])
            ss = kk if kk <= n_sites//2 else kk - n_sites
            spins.add(ss)

        # Identify level
        note = ""
        if abs(R) < 0.01: note = "GS"
        elif abs(R - 1.0) < 0.1: note = f"σ primary (x≈{x1_approx:.3f})"
        # Check for descendant: should be at R ≈ 1 + 1/x₁
        desc_R_pred = 1.0 + 1.0/x1_approx
        if abs(R - desc_R_pred) < 1.0 and R > 2:
            note = f"L₋₁σ? (pred R={desc_R_pred:.1f})"

        spin_str = ",".join(f"{s:+d}" for s in sorted(spins))
        print(f"  {i:>3} {R:>8.3f} {spin_str:>10} {x_val:>8.4f} {degen:>4} {note:>25}")

        level_info.append({
            "i": int(i), "R": round(float(R), 6),
            "x": round(float(x_val), 6), "degen": degen,
            "spins": sorted(int(s) for s in spins)
        })
        i = j

    results_dict[f"n={n_sites}"] = {
        "n": n_sites, "dim": dim, "time_s": round(dt, 2),
        "E0": round(float(evals[0]), 8),
        "delta1": round(float(delta1), 8),
        "delta1_N": round(float(delta1 * n_sites), 6),
        "levels": level_info
    }
    return level_info

all_results = {}

# ========== q=4 ==========
print("=" * 70)
print("CONFORMAL TOWER: q=4 POTTS at g_c=0.392")
print("Spin fields: 3 = (σ,σ³) pair + σ² self-conjugate")
print(f"x₁ ≈ 0.117, predicted desc gap R = 1 + 1/0.117 ≈ 9.5")
print("=" * 70)

all_results["q=4"] = {}
for n in [4, 6, 8]:  # q=4: 4^8 = 65536 manageable
    analyze_tower(n, q=4, g_c=0.392, x1_approx=0.117, results_dict=all_results["q=4"])

# Also do n=10 (4^10 = 1M) — time check
print("\n  Attempting n=10 (dim=1,048,576)...")
t0 = time.time()
try:
    analyze_tower(10, q=4, g_c=0.392, x1_approx=0.117, results_dict=all_results["q=4"], max_eig=25)
except Exception as e:
    print(f"  n=10 failed: {e}")
    print(f"  Time before failure: {time.time()-t0:.1f}s")

# ========== q=5 ==========
print("\n" + "=" * 70)
print("CONFORMAL TOWER: q=5 POTTS at g_c=0.441")
print("Spin fields: 4 = (σ,σ⁴) + (σ²,σ³) pairs")
print(f"x₁ ≈ 0.101, predicted desc gap R = 1 + 1/0.101 ≈ 10.9")
print("=" * 70)

all_results["q=5"] = {}
for n in [4, 6, 8]:  # q=5: 5^8 = 390625
    analyze_tower(n, q=5, g_c=0.441, x1_approx=0.101, results_dict=all_results["q=5"])

# ========== Summary ==========
print("\n\n" + "=" * 70)
print("DESCENDANT GAP ANALYSIS")
print("=" * 70)

for q_label, q_data in all_results.items():
    print(f"\n{q_label}:")
    for n_label, ndata in q_data.items():
        levels = ndata["levels"]
        sigma_R = [l for l in levels if 0.9 < l["R"] < 1.1]
        # Find level with spin ±1 above R=5
        desc_candidates = [l for l in levels if l["R"] > 5 and (-1 in l["spins"] or 1 in l["spins"])]
        if sigma_R and desc_candidates:
            desc = desc_candidates[0]
            gap = desc["R"] - sigma_R[0]["R"]
            print(f"  {n_label}: R(σ)={sigma_R[0]['R']:.3f}, R(desc)={desc['R']:.3f}, "
                  f"gap={gap:.3f}, spins={desc['spins']}, degen={desc['degen']}")

with open("results/sprint_059c_q4q5_tower.json", "w") as f:
    json.dump(all_results, f, indent=2)
print("\nResults saved.")
