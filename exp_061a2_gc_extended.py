"""Sprint 061a2: Extended g_c scan for q=15,20 — wider range and n=4,6 pairs.

Previous scan showed no crossing for n=4,5 in [0.70, 0.90].
n=5 consistently above n=4 — could be odd-n effect with periodic BC.
Try n=4,6 pairs instead (even-even, matching Sprint 051-052 method).
q=15 n=6: dim=15^6=11.4M — TOO LARGE for exact diag.

Alternative: use n=4 Δ·N slope method to estimate g_c directly.
At g_c, Δ·N should show scale-invariance. For a single n, look for where
d(Δ·N)/dg changes character — inflection point approach.

Also try: use the fact that Δ·N at g_c should equal 2π·x₁·v_F/N → 2π·x₁.
If x₁ ≈ 0.08 (extrapolating from q=10), then Δ·N(g_c) ≈ 0.50.
But n=4 gives Δ·N ≈ 0.12 at g=0.80 — way below 0.50.
This means n=4 is FAR from the thermodynamic limit.

Strategy: Use g_c formula as best estimate, verify with n=4 spectrum ratios.
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
            left_dim = q**i if i > 0 else 1
            right_dim = q**(n-j-1) if n-j-1 > 0 else 1
            op = csr_matrix(delta_2)
            if i > 0:
                op = sp_kron(sp_eye(left_dim), op, format='csr')
            if n-j-1 > 0:
                op = sp_kron(op, sp_eye(right_dim), format='csr')
            H = H - op
        else:
            for s in range(q):
                proj = np.zeros((q, q)); proj[s, s] = 1.0
                ops = [eye_q] * n; ops[0] = proj; ops[n-1] = proj
                result = csr_matrix(ops[0])
                for op in ops[1:]:
                    result = sp_kron(result, csr_matrix(op), format='csr')
                H = H - result
    for i in range(n):
        ops = [eye_q] * n; ops[i] = XpXd
        result = csr_matrix(ops[0])
        for op in ops[1:]:
            result = sp_kron(result, csr_matrix(op), format='csr')
        H = H - g * result
    return H

def energy_gap_scan(q, n, g_values, n_eig=6):
    dim = q**n
    print(f"  q={q}, n={n}, dim={dim}")
    results = []
    for g in g_values:
        t0 = time.time()
        H = potts_hamiltonian_periodic(n, q, g)
        evals = eigsh(H, k=n_eig, which='SA', return_eigenvectors=False)
        evals.sort()
        gap = evals[1] - evals[0]
        delta_N = gap * n
        dt = time.time() - t0
        results.append({"g": g, "E0": float(evals[0]), "gap": float(gap),
                        "delta_N": float(delta_N), "time_s": round(dt, 2)})
        print(f"    g={g:.3f}: Δ·N={delta_N:.6f} ({dt:.1f}s)")
    return results

# ===== Strategy: compare n=4 Δ·N to known q values =====
# For q=10 (known g_c=0.684), Δ·N at n=4 was calibrated.
# Use the RATIO of Δ·N(n=4)/Δ·N(n=6) at g_c for known q to extrapolate.

# Actually, let's try a wider scan for q=15 n=4,5 to find the crossing
print("=" * 60)
print("q=15: Extended scan n=4,5 from g=0.30 to g=1.00")
print("=" * 60)

g_wide = np.arange(0.30, 1.02, 0.05)
print("\nn=4:")
r15_n4 = energy_gap_scan(15, 4, g_wide, n_eig=4)
print("\nn=5:")
r15_n5 = energy_gap_scan(15, 5, g_wide, n_eig=4)

# Check for crossing
print("\n  g      Δ·N(n=4)    Δ·N(n=5)    diff(n5-n4)")
for r4, r5 in zip(r15_n4, r15_n5):
    d = r5["delta_N"] - r4["delta_N"]
    print(f"  {r4['g']:.3f}  {r4['delta_N']:.6f}  {r5['delta_N']:.6f}  {d:+.6f}")

# Also try q=10 as calibration
print("\n" + "=" * 60)
print("q=10 calibration: n=4,5 scan near g_c=0.684")
print("=" * 60)
g_cal = np.arange(0.50, 0.85, 0.05)
print("\nn=4:")
r10_n4 = energy_gap_scan(10, 4, g_cal, n_eig=4)
print("\nn=5:")
r10_n5 = energy_gap_scan(10, 5, g_cal, n_eig=4)

print("\n  g      Δ·N(n=4)    Δ·N(n=5)    diff")
for r4, r5 in zip(r10_n4, r10_n5):
    d = r5["delta_N"] - r4["delta_N"]
    print(f"  {r4['g']:.3f}  {r4['delta_N']:.6f}  {r5['delta_N']:.6f}  {d:+.6f}")

# Find crossing for q=10
for i in range(len(r10_n4) - 1):
    d1 = r10_n5[i]["delta_N"] - r10_n4[i]["delta_N"]
    d2 = r10_n5[i+1]["delta_N"] - r10_n4[i+1]["delta_N"]
    if d1 * d2 <= 0:
        g1, g2 = r10_n4[i]["g"], r10_n4[i+1]["g"]
        gc = g1 - d1 * (g2 - g1) / (d2 - d1)
        print(f"\n  q=10 crossing (n=4,5): g_c_raw = {gc:.4f} (known: 0.684)")

all_results = {
    "q15_n4": r15_n4, "q15_n5": r15_n5,
    "q10_n4_cal": r10_n4, "q10_n5_cal": r10_n5,
}

with open("results/sprint_061a2_gc_extended.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print("\nResults saved.")
