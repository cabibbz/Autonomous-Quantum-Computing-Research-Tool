"""Sprint 061a: Energy gap g_c for q=15,20 Potts via Δ·N crossing.

q=15 at n=4,5,6: dim = 50625, 759375, 11390625 (n=6 too large → n=4,5)
q=20 at n=4,5: dim = 160000, 3200000 (n=5 borderline)

Formula prediction: g_c(15) ≈ 0.798, g_c(20) ≈ 0.922.
FSS correction for n=4,5 pairs: ~5% upward (calibrated from q=2,3).
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import eigsh
import json, time

def potts_hamiltonian_periodic(n, q, g):
    """H = -J*sum delta(s_i, s_j) - g*sum (X + X†)"""
    dim = q**n
    eye_q = np.eye(q)

    # Clock shift operator X and X+X†
    X = np.zeros((q, q))
    for s in range(q):
        X[(s+1) % q, s] = 1.0
    XpXd = X + X.T

    # Potts interaction: delta(s_i, s_j) = sum_s |s><s| ⊗ |s><s|
    delta_2 = np.zeros((q**2, q**2))
    for s in range(q):
        delta_2[s*q+s, s*q+s] = 1.0

    H = csr_matrix((dim, dim))

    # Nearest-neighbor interaction (periodic)
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
            # Wrap-around bond (site n-1 to site 0)
            for s in range(q):
                proj = np.zeros((q, q)); proj[s, s] = 1.0
                # Build tensor product manually
                ops = [eye_q] * n
                ops[0] = proj
                ops[n-1] = proj
                result = csr_matrix(ops[0])
                for op in ops[1:]:
                    result = sp_kron(result, csr_matrix(op), format='csr')
                H = H - result

    # Transverse field
    for i in range(n):
        ops = [eye_q] * n
        ops[i] = XpXd
        result = csr_matrix(ops[0])
        for op in ops[1:]:
            result = sp_kron(result, csr_matrix(op), format='csr')
        H = H - g * result

    return H

def energy_gap_scan(q, n, g_values, n_eig=4):
    """Compute E0, E1 at each g, return Δ·N."""
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
        results.append({"g": g, "E0": float(evals[0]), "E1": float(evals[1]),
                        "gap": float(gap), "delta_N": float(delta_N), "time_s": round(dt, 2)})
        print(f"    g={g:.3f}: Δ·N={delta_N:.4f} ({dt:.1f}s)")
    return results

def find_crossing(results_n1, results_n2):
    """Find where Δ·N curves cross via linear interpolation."""
    g_vals = [r["g"] for r in results_n1]
    dN1 = [r["delta_N"] for r in results_n1]
    dN2 = [r["delta_N"] for r in results_n2]

    diff = [d2 - d1 for d1, d2 in zip(dN1, dN2)]

    for i in range(len(diff) - 1):
        if diff[i] * diff[i+1] <= 0:
            # Linear interpolation
            g_cross = g_vals[i] - diff[i] * (g_vals[i+1] - g_vals[i]) / (diff[i+1] - diff[i])
            return g_cross
    return None

all_results = {}

# ===== q=15 =====
print("=" * 60)
print("q=15 POTTS — formula predicts g_c ≈ 0.798")
print("=" * 60)

# Coarse scan first
g_coarse = np.arange(0.60, 1.00, 0.05)
print("\nCoarse scan n=4:")
r15_n4_coarse = energy_gap_scan(15, 4, g_coarse, n_eig=4)

# Fine scan around predicted g_c
g_fine = np.arange(0.70, 0.90, 0.02)
print("\nFine scan n=4:")
r15_n4 = energy_gap_scan(15, 4, g_fine, n_eig=4)
print("\nFine scan n=5:")
t0 = time.time()
r15_n5 = energy_gap_scan(15, 5, g_fine, n_eig=4)
t_n5 = time.time() - t0
print(f"  Total n=5 time: {t_n5:.1f}s")

g_c_15_raw = find_crossing(r15_n4, r15_n5)
print(f"\n  Crossing (n=4,5): g_c_raw = {g_c_15_raw:.4f}" if g_c_15_raw else "\n  No crossing found!")

# Apply FSS correction (calibrated ~5% for small n pairs)
if g_c_15_raw:
    g_c_15_corr = g_c_15_raw * 1.048  # 4.8% correction for n=4,5
    print(f"  Corrected g_c(15) = {g_c_15_corr:.4f}")
    print(f"  Formula prediction: 0.798")
    print(f"  Deviation: {abs(g_c_15_corr - 0.798)/0.798*100:.1f}%")

all_results["q15"] = {
    "n4_fine": r15_n4,
    "n5_fine": r15_n5,
    "g_c_raw": g_c_15_raw,
    "g_c_corrected": g_c_15_corr if g_c_15_raw else None,
    "formula_prediction": 0.798,
}

# ===== q=20 =====
print("\n" + "=" * 60)
print("q=20 POTTS — formula predicts g_c ≈ 0.922")
print("=" * 60)

# n=4 only (n=5 is 3.2M dim, might be too slow)
g_fine_20 = np.arange(0.80, 1.05, 0.02)
print("\nFine scan n=4:")
r20_n4 = energy_gap_scan(20, 4, g_fine_20, n_eig=4)

# Try n=5 with a single timing test
print("\nTiming test n=5, single point...")
t0 = time.time()
try:
    r20_n5_test = energy_gap_scan(20, 5, [0.92], n_eig=4)
    t_test = time.time() - t0
    print(f"  Single point: {t_test:.1f}s")

    if t_test < 20:  # If fast enough, do full scan
        print("\nFull scan n=5:")
        r20_n5 = energy_gap_scan(20, 5, g_fine_20, n_eig=4)
        g_c_20_raw = find_crossing(r20_n4, r20_n5)
        if g_c_20_raw:
            g_c_20_corr = g_c_20_raw * 1.048
            print(f"\n  Crossing (n=4,5): g_c_raw = {g_c_20_raw:.4f}")
            print(f"  Corrected g_c(20) = {g_c_20_corr:.4f}")
        else:
            print("\n  No crossing found for q=20!")
            g_c_20_raw = None
            g_c_20_corr = None
            r20_n5 = r20_n5_test
    else:
        print(f"  Too slow ({t_test:.0f}s), using n=4 only")
        r20_n5 = r20_n5_test
        g_c_20_raw = None
        g_c_20_corr = None
except Exception as e:
    print(f"  n=5 failed: {e}")
    r20_n5 = []
    g_c_20_raw = None
    g_c_20_corr = None

all_results["q20"] = {
    "n4_fine": r20_n4,
    "n5_fine": r20_n5 if isinstance(r20_n5, list) else r20_n5,
    "g_c_raw": g_c_20_raw,
    "g_c_corrected": g_c_20_corr,
    "formula_prediction": 0.922,
}

# Summary
print("\n" + "=" * 60)
print("SUMMARY — g_c(q) extended")
print("=" * 60)
print(f"  q=15: formula={0.798:.3f}, measured={g_c_15_corr:.3f}" if g_c_15_raw else "  q=15: no crossing")
if g_c_20_corr:
    print(f"  q=20: formula={0.922:.3f}, measured={g_c_20_corr:.3f}")
else:
    print(f"  q=20: formula={0.922:.3f}, crossing not available (n=4 only)")

with open("results/sprint_061a_gc_large_q.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print("\nResults saved to results/sprint_061a_gc_large_q.json")
