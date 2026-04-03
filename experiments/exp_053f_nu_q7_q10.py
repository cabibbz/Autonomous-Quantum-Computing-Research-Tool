#!/usr/bin/env python3
"""Sprint 053f: ν(q=7) and ν(q=10) — push exact diag limits.

q=7: n=4 (2401) done. Try n=6 (117649) — sparse should work.
q=10: n=4 (10000). Try n=6 (1000000) — probably too large.

Also: refine q=4 with n=10 (4^10=1048576, too large) — stick with n=4,6,8.
"""
import numpy as np, json, time, warnings
warnings.filterwarnings('ignore')
from scipy.sparse.linalg import eigsh
from scipy.sparse import kron, csr_matrix
from scipy.interpolate import interp1d

def potts_hamiltonian(q, n, g, J=1.0):
    dim = q**n
    X = np.zeros((q, q), dtype=complex)
    for a in range(q):
        X[(a+1) % q, a] = 1.0
    Xphc = X + X.conj().T
    projectors = [csr_matrix(np.diag([1.0 if b == a else 0.0 for b in range(q)])) for a in range(q)]
    H = csr_matrix((dim, dim), dtype=complex)
    I = csr_matrix(np.eye(q))
    for i in range(n - 1):
        for a in range(q):
            op = csr_matrix(np.eye(1))
            for j in range(n):
                if j == i:
                    op = kron(op, projectors[a])
                elif j == i + 1:
                    op = kron(op, projectors[a])
                else:
                    op = kron(op, I)
            H += -J * op
    Xphc_sp = csr_matrix(Xphc)
    for i in range(n):
        op = csr_matrix(np.eye(1))
        for j in range(n):
            if j == i:
                op = kron(op, Xphc_sp)
            else:
                op = kron(op, I)
        H += -g * op
    return H

def energy_gap(q, n, g, J=1.0):
    H = potts_hamiltonian(q, n, g, J)
    vals = eigsh(H, k=min(6, H.shape[0]-1), which='SA', return_eigenvectors=False)
    vals = np.sort(vals.real)
    return vals[1] - vals[0]

b_cal = 0.86

def nu_corrected(s1, s2, n1, n2, b=0.86):
    corr = (1 + b/n2) / (1 + b/n1)
    eff = (s2 / s1) / corr
    if eff <= 1:
        return None
    return 1.0 / (np.log(eff) / np.log(n2/n1))

# === q=7, n=6 timing test ===
print("=== q=7, n=6 (dim=117649) timing ===", flush=True)
t0 = time.time()
try:
    gap = energy_gap(7, 6, 0.535)
    dt = time.time() - t0
    print(f"  gap={gap:.6f}, t={dt:.1f}s", flush=True)
    q7_n6_ok = dt < 30  # Need ~20 points, so 30s each = 10 min max
except Exception as e:
    print(f"  FAILED: {e}", flush=True)
    q7_n6_ok = False

results = {}

# === q=7 ===
if q7_n6_ok:
    g_c_q7 = 0.535
    g_values = np.linspace(0.40, 0.68, 29)

    for n in [4, 6]:
        dim = 7**n
        print(f"\n=== q=7 n={n} (dim={dim}): {len(g_values)} points ===", flush=True)
        t0 = time.time()
        gaps = []
        for g in g_values:
            if time.time() - t0 > 200:
                print(f"  Time limit at g={g:.3f}", flush=True)
                break
            gap = energy_gap(7, n, g)
            gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * n)})
        dt = time.time() - t0
        print(f"  Done: {len(gaps)} points in {dt:.1f}s", flush=True)
        results[f'q7_n{n}'] = gaps

    # Crossing
    if 'q7_n4' in results and 'q7_n6' in results:
        d4 = results['q7_n4']
        d6 = results['q7_n6']
        min_len = min(len(d4), len(d6))
        for j in range(min_len - 1):
            diff_j = d4[j]['gap_x_n'] - d6[j]['gap_x_n']
            diff_j1 = d4[j+1]['gap_x_n'] - d6[j+1]['gap_x_n']
            if diff_j * diff_j1 < 0:
                gc = d4[j]['g'] - diff_j * (d4[j+1]['g'] - d4[j]['g']) / (diff_j1 - diff_j)
                print(f"\n  q=7 crossing (4,6): g = {gc:.5f}", flush=True)

        # Slopes
        g_arr4 = np.array([p['g'] for p in d4])
        y_arr4 = np.array([p['gap_x_n'] for p in d4])
        g_arr6 = np.array([p['g'] for p in d6])
        y_arr6 = np.array([p['gap_x_n'] for p in d6])
        f4 = interp1d(g_arr4, y_arr4, kind='cubic')
        f6 = interp1d(g_arr6, y_arr6, kind='cubic')

        for gc_test in [0.52, 0.53, 0.535, 0.54, 0.55]:
            dg = 0.005
            if gc_test - dg < max(g_arr4.min(), g_arr6.min()) or gc_test + dg > min(g_arr4.max(), g_arr6.max()):
                continue
            s4 = float((f4(gc_test + dg) - f4(gc_test - dg)) / (2 * dg))
            s6 = float((f6(gc_test + dg) - f6(gc_test - dg)) / (2 * dg))
            ratio = s6 / s4
            nu_raw = np.log(6/4) / np.log(ratio)
            nu_c = nu_corrected(s4, s6, 4, 6)
            print(f"  g_c={gc_test:.3f}: s4={s4:.3f}, s6={s6:.3f}, raw ν={nu_raw:.4f}, corr ν={nu_c:.4f}", flush=True)

# === q=10, n=4 ===
print(f"\n=== q=10, n=4 (dim={10**4}) ===", flush=True)
g_c_q10 = 0.684
g_values_q10 = np.linspace(0.50, 0.85, 36)
t0 = time.time()
gaps_q10_n4 = []
for g in g_values_q10:
    gap = energy_gap(10, 4, g)
    gaps_q10_n4.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * 4)})
print(f"  Done: {len(gaps_q10_n4)} points in {time.time()-t0:.1f}s", flush=True)
results['q10_n4'] = gaps_q10_n4

# Try q=10, n=5 (dim=100000, should be feasible)
print(f"\n=== q=10, n=5 (dim={10**5}) timing ===", flush=True)
t0 = time.time()
try:
    gap = energy_gap(10, 5, 0.684)
    dt = time.time() - t0
    print(f"  gap={gap:.6f}, t={dt:.1f}s", flush=True)
    if dt < 30:
        print(f"\n=== q=10 n=5: {len(g_values_q10)} points ===", flush=True)
        t0 = time.time()
        gaps_q10_n5 = []
        for g in g_values_q10:
            if time.time() - t0 > 200:
                print(f"  Time limit", flush=True)
                break
            gap = energy_gap(10, 5, g)
            gaps_q10_n5.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * 5)})
        print(f"  Done: {len(gaps_q10_n5)} points in {time.time()-t0:.1f}s", flush=True)
        results['q10_n5'] = gaps_q10_n5

        # Crossing and slopes
        d4 = results['q10_n4']
        d5 = results['q10_n5']
        min_len = min(len(d4), len(d5))
        for j in range(min_len - 1):
            diff_j = d4[j]['gap_x_n'] - d5[j]['gap_x_n']
            diff_j1 = d4[j+1]['gap_x_n'] - d5[j+1]['gap_x_n']
            if diff_j * diff_j1 < 0:
                gc = d4[j]['g'] - diff_j * (d4[j+1]['g'] - d4[j]['g']) / (diff_j1 - diff_j)
                print(f"\n  q=10 crossing (4,5): g = {gc:.5f}", flush=True)

        g_arr4 = np.array([p['g'] for p in d4[:min_len]])
        y_arr4 = np.array([p['gap_x_n'] for p in d4[:min_len]])
        g_arr5 = np.array([p['g'] for p in d5])
        y_arr5 = np.array([p['gap_x_n'] for p in d5])
        f4 = interp1d(g_arr4, y_arr4, kind='cubic')
        f5 = interp1d(g_arr5, y_arr5, kind='cubic')

        for gc_test in [0.65, 0.67, 0.684, 0.70]:
            dg = 0.005
            if gc_test - dg < max(g_arr4.min(), g_arr5.min()) or gc_test + dg > min(g_arr4.max(), g_arr5.max()):
                continue
            s4 = float((f4(gc_test + dg) - f4(gc_test - dg)) / (2 * dg))
            s5 = float((f5(gc_test + dg) - f5(gc_test - dg)) / (2 * dg))
            ratio = s5 / s4
            nu_raw = np.log(5/4) / np.log(ratio)
            nu_c = nu_corrected(s4, s5, 4, 5)
            print(f"  g_c={gc_test:.3f}: s4={s4:.3f}, s5={s5:.3f}, raw ν={nu_raw:.4f}, corr ν={nu_c:.4f}", flush=True)

except Exception as e:
    print(f"  FAILED: {e}", flush=True)

# Save
with open('results/sprint_053f_nu_q7_q10.json', 'w') as f:
    json.dump({'sprint': '053f', 'b_cal': b_cal, 'data': results}, f, indent=2)

print("\nDone!", flush=True)
