"""Sprint 059c (v2): Conformal towers for q=4,5,7,10 — all novel q values.
Reduced eigenvalue count for large systems. Uses existing q=4 n=4,6,8 data.
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

def analyze_tower(n_sites, q, g_c, x1_approx, max_eig=20):
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
    print(f"  {'i':>3} {'R':>8} {'spins':>10} {'x=R·x₁':>8} {'deg':>4}")
    print(f"  {'-'*40}")

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
        spin_str = ",".join(f"{s:+d}" for s in sorted(spins))
        print(f"  {i:>3} {R:>8.3f} {spin_str:>10} {x_val:>8.4f} {degen:>4}")
        level_info.append({
            "i": int(i), "R": round(float(R), 6),
            "x": round(float(x_val), 6), "degen": degen,
            "spins": sorted(int(s) for s in spins)
        })
        i = j

    return {
        "n": n_sites, "dim": dim, "time_s": round(dt, 2),
        "E0": round(float(evals[0]), 8),
        "delta1": round(float(delta1), 8),
        "delta1_N": round(float(delta1 * n_sites), 6),
        "levels": level_info
    }

all_results = {}

# ========== q=5 ==========
print("=" * 70)
print("q=5 POTTS at g_c=0.441")
print("Spin fields: 4 = (σ,σ⁴) + (σ²,σ³) pairs")
print(f"x₁ ≈ 0.101, 1/x₁ ≈ 9.9, desc at R ≈ 10.9")
print("=" * 70)
all_results["q=5"] = {}
for n in [4, 6, 8]:  # 5^8 = 390k
    all_results["q=5"][f"n={n}"] = analyze_tower(n, 5, 0.441, 0.101, max_eig=20)

# ========== q=7 ==========
print("\n" + "=" * 70)
print("q=7 POTTS at g_c=0.535")
print("Spin fields: 6 = 3 pairs")
print(f"x₁ ≈ 0.086, 1/x₁ ≈ 11.6, desc at R ≈ 12.6")
print("=" * 70)
all_results["q=7"] = {}
for n in [4, 6]:  # 7^6 = 117k, 7^8 = 5.7M too big
    all_results["q=7"][f"n={n}"] = analyze_tower(n, 7, 0.535, 0.086, max_eig=20)

# ========== q=10 ==========
print("\n" + "=" * 70)
print("q=10 POTTS at g_c=0.684")
print("Spin fields: 9 = 4 pairs + 1 self-conjugate")
print(f"x₁ ≈ 0.083, 1/x₁ ≈ 12.0, desc at R ≈ 13.0")
print("=" * 70)
all_results["q=10"] = {}
for n in [4, 6]:  # 10^6 = 1M
    me = 15 if n == 6 else 20
    all_results["q=10"][f"n={n}"] = analyze_tower(n, 10, 0.684, 0.083, max_eig=me)

# ========== DESCENDANT GAP SUMMARY ==========
print("\n\n" + "=" * 70)
print("DESCENDANT GAP SUMMARY — KEY RESULT")
print("=" * 70)
print(f"{'q':>3} {'n':>4} {'R(σ)':>7} {'R(desc)':>8} {'gap':>7} {'pred':>7} {'desc spins':>12} {'desc deg':>9}")
print("-" * 65)

# Include q=4 n=8 from earlier run
q4_data = {
    "q": 4, "n": 8, "R_sigma": 1.000, "R_desc": 9.018,
    "gap": 8.018, "pred_gap": 1/0.117, "desc_spins": [-1, 1], "desc_degen": 4
}
print(f"{q4_data['q']:>3} {q4_data['n']:>4} {q4_data['R_sigma']:>7.3f} {q4_data['R_desc']:>8.3f} "
      f"{q4_data['gap']:>7.3f} {q4_data['pred_gap']:>7.1f} {str(q4_data['desc_spins']):>12} {q4_data['desc_degen']:>9}")

for q_label, q_data in all_results.items():
    q_num = int(q_label.split("=")[1])
    x1_map = {5: 0.101, 7: 0.086, 10: 0.083}
    x1 = x1_map[q_num]
    pred_gap = 1.0 / x1

    # Use largest n
    largest_n = max(q_data.keys(), key=lambda k: q_data[k]["n"])
    ndata = q_data[largest_n]
    levels = ndata["levels"]

    sigma = [l for l in levels if 0.9 < l["R"] < 1.1]
    # Find first level with spin ±1 that's above R=5
    desc_candidates = [l for l in levels if l["R"] > 5 and (-1 in l["spins"] or 1 in l["spins"])]

    if sigma and desc_candidates:
        desc = desc_candidates[0]
        gap = desc["R"] - sigma[0]["R"]
        print(f"{q_num:>3} {ndata['n']:>4} {sigma[0]['R']:>7.3f} {desc['R']:>8.3f} "
              f"{gap:>7.3f} {pred_gap:>7.1f} {str(desc['spins']):>12} {desc['degen']:>9}")

# Overall check
print("\n\n--- CFT TOWER TEST: Does R(desc)-R(σ) ≈ 1/x₁? ---")
print("q=2: gap=7.898, pred=8.0 (1/0.125). Ratio=0.987")
print("q=3: gap=7.332, pred=7.5 (1/0.133). Ratio=0.978")
print("q=4: gap=8.018, pred=8.5 (1/0.117). Ratio=0.943")
print("(q=5,7,10 from this run)")

with open("results/sprint_059c_tower_all.json", "w") as f:
    json.dump(all_results, f, indent=2)
print("\nResults saved.")
