"""Sprint 057f: q=10 spectrum + absolute scaling dimension extraction.

q=10 at n=4 (dim=10000): expect 4 conjugate pairs + 1 self-conjugate = 9 spin fields.
Also extract absolute x_1 values using two-size c/x_1 method.
"""
import numpy as np
from scipy.sparse import csr_matrix, eye as sp_eye, kron as sp_kron
from scipy.sparse.linalg import eigsh
import json, time

def potts_hamiltonian_periodic(n, q, g):
    d = q
    dim = d**n
    X = np.zeros((d, d))
    for s in range(d):
        X[(s+1) % d, s] = 1.0
    XpXd = csr_matrix(X + X.T)
    delta_op = np.zeros((d*d, d*d))
    for a in range(d):
        idx = a * d + a
        delta_op[idx, idx] = 1.0
    delta_sp = csr_matrix(delta_op)
    H = csr_matrix((dim, dim))
    for i in range(n):
        j = (i + 1) % n
        if j == i + 1:
            left = sp_eye(d**i, format='csr') if i > 0 else csr_matrix(np.array([[1.0]]))
            right = sp_eye(d**(n-i-2), format='csr') if i + 2 < n else csr_matrix(np.array([[1.0]]))
            H = H - sp_kron(sp_kron(left, delta_sp, format='csr'), right, format='csr')
        else:
            for a in range(d):
                proj0 = np.zeros((d, d)); proj0[a, a] = 1.0
                projn = np.zeros((d, d)); projn[a, a] = 1.0
                op = sp_kron(sp_kron(csr_matrix(proj0), sp_eye(d**(n-2), format='csr'), format='csr'),
                             csr_matrix(projn), format='csr')
                H = H - op
    for i in range(n):
        left = sp_eye(d**i, format='csr') if i > 0 else csr_matrix(np.array([[1.0]]))
        right = sp_eye(d**(n-i-1), format='csr') if i < n-1 else csr_matrix(np.array([[1.0]]))
        H = H - g * sp_kron(sp_kron(left, XpXd, format='csr'), right, format='csr')
    return H

results = {}

# q=10, n=4
print("=== q=10, n=4 (dim=10000) ===")
t0 = time.time()
H = potts_hamiltonian_periodic(4, 10, 0.684)
n_eig = 25
evals, _ = eigsh(H, k=n_eig, which='SA')
evals = np.sort(evals)
dt = time.time() - t0
gaps = evals - evals[0]
nonzero = gaps[gaps > 1e-8]
ratios = nonzero / nonzero[0]
print(f"  Time: {dt:.1f}s")
print(f"  Ground state degeneracy: {np.sum(gaps < 1e-6)}")
print(f"  Gap ratios (first 15):")
for i, r in enumerate(ratios[:20]):
    print(f"    R_{i+1} = {r:.4f}")

results["q=10_n=4"] = {
    "q": 10, "n": 4, "g_c": 0.684, "dim": 10000, "time_s": round(dt, 2),
    "eigenvalues": [round(float(e), 8) for e in evals],
    "nonzero_gaps": [round(float(g), 8) for g in nonzero[:20]],
    "ratios": [round(float(r), 6) for r in ratios[:20]],
}

# --- Extract absolute scaling dimensions using c/x_1 ---
print("\n\n=== ABSOLUTE SCALING DIMENSIONS ===")
print("Method: E₀(n)/n ≈ e_∞ - πcv/(6n²), Δ₁(n) = 2πvx₁/n")
print("From two sizes: c/x₁ = [E₀(n₁)/n₁ - E₀(n₂)/n₂] × 12πn₁²n₂² / [Δ̄₁·n̄ × (n₂² - n₁²)]")

# Load all data
import glob
all_data = {}

# q=2 Ising
ising = json.load(open("results/sprint_057a_ising_spectrum.json"))
for key in ising:
    n = ising[key]["n"]
    all_data.setdefault(2, {})[n] = {
        "E0": ising[key]["E0"],
        "gap1": ising[key]["gaps"][1]
    }

# q=3
q3 = json.load(open("results/sprint_057b_q3_spectrum.json"))
for key in q3:
    n = q3[key]["n"]
    all_data.setdefault(3, {})[n] = {
        "E0": q3[key]["E0"],
        "gap1": q3[key]["nonzero_gaps"][0]
    }

# q=4,5
q45 = json.load(open("results/sprint_057c_q4q5_spectrum.json"))
for key in q45:
    q = q45[key]["q"]
    n = q45[key]["n"]
    all_data.setdefault(q, {})[n] = {
        "E0": q45[key]["E0"],
        "gap1": q45[key]["nonzero_gaps"][0]
    }

# q=4 n=8 and q=7
q47 = json.load(open("results/sprint_057d_q7_spectrum.json"))
for key in q47:
    q_val = q47[key]["q"]
    n = q47[key]["n"]
    all_data.setdefault(q_val, {})[n] = {
        "E0": q47[key]["eigenvalues"][0],
        "gap1": q47[key]["nonzero_gaps"][0]
    }

# Known central charges
c_known = {2: 0.500, 3: 0.800, 4: 1.000, 5: 1.10, 7: 1.30, 10: 1.40}

print(f"\n{'q':>3} {'sizes':>10} {'c':>6} {'c/x₁':>8} {'x₁':>8} {'x₁ exact':>10}")
print("-" * 55)

x1_values = {}
for q in sorted(all_data.keys()):
    sizes = sorted(all_data[q].keys())
    if len(sizes) < 2:
        continue

    # Use two largest sizes
    n1, n2 = sizes[-2], sizes[-1]
    e1 = all_data[q][n1]["E0"] / n1
    e2 = all_data[q][n2]["E0"] / n2
    gap1 = all_data[q][n1]["gap1"]
    gap2 = all_data[q][n2]["gap1"]

    # v·x₁ from each size
    vx1_1 = gap1 * n1 / (2 * np.pi)
    vx1_2 = gap2 * n2 / (2 * np.pi)
    vx1_avg = (vx1_1 + vx1_2) / 2

    # c·v from two sizes
    cv = (e1 - e2) * 6 / (np.pi * (1/n2**2 - 1/n1**2))

    # c/x₁ = cv / vx₁
    cx1_ratio = cv / vx1_avg

    c = c_known.get(q, None)
    if c:
        x1 = c / cx1_ratio
        exact = {2: 0.125, 3: 2/15, 4: 0.25}
        ex_str = f"{exact.get(q, '—'):.4f}" if q in exact else "—"
        print(f"{q:>3} {f'{n1},{n2}':>10} {c:>6.2f} {cx1_ratio:>8.3f} {x1:>8.4f} {ex_str:>10}")
        x1_values[q] = x1

results["absolute_x1"] = {str(q): round(float(v), 6) for q, v in x1_values.items()}

with open("results/sprint_057f_q10_absolute.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nResults saved.")
