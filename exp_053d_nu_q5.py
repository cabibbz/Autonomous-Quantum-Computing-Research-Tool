#!/usr/bin/env python3
"""Sprint 053d: ν(q=5) at true g_c ≈ 0.44 using validated methods.

Previous Sprint 045 found ν≈2.0 at g≈0.45, close to true g_c≈0.44.
Sprint 052 corrected g_c(5) = 0.441 (energy gap method).

Use exact diag for n=4,6 (q=5: 625, 15625) — n=8 would be 5^8=390625, too large.
The slope ratio from just 2 sizes is noisier, so we'll also try DMRG for larger sizes.

Methods:
1. Exact diag n=4,6: energy gap slopes and crossing
2. Refine g_c via gap crossings at finer resolution
3. Slope ratio → ν estimate with 1/N correction from q=3 calibration
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

# === Timing ===
print("=== Timing tests ===", flush=True)
for n in [4, 6]:
    t0 = time.time()
    gap = energy_gap(5, n, 0.44)
    dt = time.time() - t0
    print(f"  q=5 n={n} (dim={5**n}): gap={gap:.6f}, t={dt:.1f}s", flush=True)

# === q=5: energy gap scan near g_c ≈ 0.44 ===
g_c_approx = 0.441  # From Sprint 052
g_values = np.linspace(0.30, 0.58, 57)

results = {}
for n in [4, 6]:
    dim = 5**n
    print(f"\n=== q=5 n={n} (dim={dim}): {len(g_values)} points ===", flush=True)
    t0 = time.time()
    gaps = []
    for g in g_values:
        gap = energy_gap(5, n, g)
        gaps.append({'g': float(g), 'gap': float(gap), 'gap_x_n': float(gap * n)})
    dt = time.time() - t0
    print(f"  Done in {dt:.1f}s", flush=True)
    results[f'n{n}'] = gaps

# Save
with open('results/sprint_053d_nu_q5.json', 'w') as f:
    json.dump({'sprint': '053d', 'q': 5, 'data': results}, f, indent=2)

# === Find crossing point ===
print("\n=== Δ·N crossing (n=4,6) ===", flush=True)
d4, d6 = results['n4'], results['n6']
g_arr = np.array([p['g'] for p in d4])
y4 = np.array([p['gap_x_n'] for p in d4])
y6 = np.array([p['gap_x_n'] for p in d6])

diff = y4 - y6
for j in range(len(diff) - 1):
    if diff[j] * diff[j+1] < 0:
        gc = g_arr[j] - diff[j] * (g_arr[j+1] - g_arr[j]) / (diff[j+1] - diff[j])
        print(f"  Crossing at g = {gc:.5f}", flush=True)

# === Slope at various g_c candidates ===
print("\n=== Slope ratios at different g_c ===", flush=True)
f4 = interp1d(g_arr, y4, kind='cubic')
f6 = interp1d(g_arr, y6, kind='cubic')

for gc_test in [0.42, 0.43, 0.44, 0.441, 0.45, 0.46]:
    dg = 0.005
    if gc_test - dg < g_arr.min() or gc_test + dg > g_arr.max():
        continue
    s4 = (f4(gc_test + dg) - f4(gc_test - dg)) / (2 * dg)
    s6 = (f6(gc_test + dg) - f6(gc_test - dg)) / (2 * dg)
    if s4 > 0 and s6 > 0:
        ratio = s6 / s4
        nu_raw = np.log(6/4) / np.log(ratio)
        print(f"  g_c={gc_test:.3f}: s4={s4:.3f}, s6={s6:.3f}, ratio={ratio:.4f}, ν_raw={nu_raw:.4f}", flush=True)

# === Apply calibration from q=3 ===
# At q=3, n=4,6 raw ν = 0.986, true ν = 5/6 = 0.833
# Correction factor: ν_true = ν_raw * (5/6) / 0.986 ≈ ν_raw * 0.845
# But better: from corrected fit, 1/N extrapolation reduces by ~13%
# For (4,6) pair: raw 0.986 → true 0.833, ratio = 0.845
print("\n=== Calibrated ν estimates ===", flush=True)
cal_factor = (5/6) / 0.986  # q=3 calibration: true/raw for (4,6) pair
print(f"  Calibration factor (from q=3): {cal_factor:.4f}", flush=True)

gc_best = 0.441
dg = 0.005
s4 = (f4(gc_best + dg) - f4(gc_best - dg)) / (2 * dg)
s6 = (f6(gc_best + dg) - f6(gc_best - dg)) / (2 * dg)
ratio = s6 / s4
nu_raw = np.log(6/4) / np.log(ratio)
nu_cal = nu_raw * cal_factor
print(f"  At g_c={gc_best}: raw ν = {nu_raw:.4f}, calibrated ν = {nu_cal:.4f}", flush=True)

# === Also check with corrected power-law approach ===
# From q=3: slope = a * N^{1/ν} * (1 + b/N) with b ≈ 0.86
# This is a 3-parameter fit with only 2 data points → underdetermined
# But we can fix b=0.86 from q=3 calibration
print("\n=== Corrected power-law (b=0.86 from q=3) ===", flush=True)
b_cal = 0.86
for gc_test in [0.43, 0.44, 0.441, 0.45]:
    dg = 0.005
    if gc_test - dg < g_arr.min() or gc_test + dg > g_arr.max():
        continue
    s4 = float((f4(gc_test + dg) - f4(gc_test - dg)) / (2 * dg))
    s6 = float((f6(gc_test + dg) - f6(gc_test - dg)) / (2 * dg))

    # slope = a * N^{1/ν} * (1 + b/N)
    # s4 = a * 4^{1/ν} * (1 + b/4)
    # s6 = a * 6^{1/ν} * (1 + b/6)
    # ratio = (6/4)^{1/ν} * (1+b/6)/(1+b/4)
    corr_ratio = (1 + b_cal/6) / (1 + b_cal/4)
    # s6/s4 = (6/4)^{1/ν} * corr_ratio
    effective_ratio = (s6 / s4) / corr_ratio
    inv_nu = np.log(effective_ratio) / np.log(6/4)
    nu_corr = 1.0 / inv_nu
    print(f"  g_c={gc_test:.3f}: effective ratio={effective_ratio:.4f}, ν={nu_corr:.4f}", flush=True)

# Save final
with open('results/sprint_053d_nu_q5.json', 'w') as f:
    json.dump({
        'sprint': '053d', 'q': 5, 'g_c_sprint052': 0.441,
        'cal_factor': float(cal_factor), 'b_correction': b_cal,
        'data': results
    }, f, indent=2)

print("\nDone!", flush=True)
