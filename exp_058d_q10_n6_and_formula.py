"""Sprint 058d: q=10 n=6 periodic spectrum + c/x1 formula search.

q=10 at n=6 (dim=10^6) to get second size for c/x1 extraction.
Then test candidate formulas for c/x1(q).
"""
import numpy as np
from scipy.sparse import csr_matrix, eye as sp_eye, kron as sp_kron
from scipy.sparse.linalg import eigsh
from scipy.optimize import curve_fit
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

# ========== PART 1: q=10, n=6 spectrum ==========
print("=" * 70)
print("PART 1: q=10 Potts, n=6 periodic (dim=1,000,000)")
print("=" * 70)

q, n, g_c = 10, 6, 0.684
t0 = time.time()
print(f"Building Hamiltonian (dim={q**n})...")
H = potts_hamiltonian_periodic(n, q, g_c)
t_build = time.time() - t0
print(f"  Build time: {t_build:.1f}s, nnz={H.nnz}")

print("Running eigsh for 6 lowest eigenvalues...")
t1 = time.time()
evals, _ = eigsh(H, k=6, which='SA')
evals = np.sort(evals)
t_eig = time.time() - t1
print(f"  Eigsh time: {t_eig:.1f}s")

gaps = evals - evals[0]
nonzero = gaps[gaps > 1e-8]

print(f"  E0 = {evals[0]:.10f}")
print(f"  E0/n = {evals[0]/n:.10f}")
print(f"  First gap = {nonzero[0]:.10f}")
print(f"  Delta1*n = {nonzero[0]*n:.6f}")
ratios = nonzero / nonzero[0]
print(f"  Gap ratios: {[f'{r:.4f}' for r in ratios]}")

results["q10_n6"] = {
    "q": 10, "n": 6, "g_c": g_c,
    "E0": float(evals[0]),
    "E0_per_site": float(evals[0]/n),
    "eigenvalues": [float(e) for e in evals],
    "nonzero_gaps": [float(g) for g in nonzero],
    "build_time_s": round(t_build, 2),
    "eigsh_time_s": round(t_eig, 2),
}
total_p1 = time.time() - t0
print(f"  Part 1 total: {total_p1:.1f}s")


# ========== PART 2: Complete c/x1 table ==========
print("\n" + "=" * 70)
print("PART 2: Complete c/x₁(q) Table")
print("=" * 70)

# All data: {q: {n: (E0, gap1)}}
all_data = {}

# q=2
ising = json.load(open("results/sprint_057a_ising_spectrum.json"))
for key in ising:
    nn = ising[key]["n"]
    all_data.setdefault(2, {})[nn] = (ising[key]["E0"], ising[key]["gaps"][1])

# q=3
q3 = json.load(open("results/sprint_057b_q3_spectrum.json"))
for key in q3:
    nn = q3[key]["n"]
    all_data.setdefault(3, {})[nn] = (q3[key]["E0"], q3[key]["nonzero_gaps"][0])

# q=4,5 Sprint 057c
q45 = json.load(open("results/sprint_057c_q4q5_spectrum.json"))
for key in q45:
    qq = q45[key]["q"]; nn = q45[key]["n"]
    all_data.setdefault(qq, {})[nn] = (q45[key]["E0"], q45[key]["nonzero_gaps"][0])

# q=4 n=8, q=7 Sprint 057d
q47 = json.load(open("results/sprint_057d_q7_spectrum.json"))
for key in q47:
    qq = q47[key]["q"]; nn = q47[key]["n"]
    all_data.setdefault(qq, {})[nn] = (q47[key]["eigenvalues"][0], q47[key]["nonzero_gaps"][0])

# q=10 n=4 Sprint 057f
q10f = json.load(open("results/sprint_057f_q10_absolute.json"))
d10 = q10f["q=10_n=4"]
all_data.setdefault(10, {})[4] = (d10["eigenvalues"][0], d10["nonzero_gaps"][0])

# New Sprint 058 data
q4_new = json.load(open("results/sprint_058a_q4_n10.json"))
all_data.setdefault(4, {})[10] = (q4_new["E0"], q4_new["nonzero_gaps"][0])

q5_new = json.load(open("results/sprint_058b_q5_n8.json"))
all_data.setdefault(5, {})[8] = (q5_new["E0"], q5_new["nonzero_gaps"][0])

# q=10 n=6 from this experiment
all_data.setdefault(10, {})[6] = (float(evals[0]), float(nonzero[0]))

# c/x1 extraction — best pair (largest sizes) for each q
c_exact = {2: 0.500, 3: 0.800, 4: 1.000}
c_est = {5: 1.10, 7: 1.30, 10: 1.40}
c_all = {**c_exact, **c_est}

print(f"\n{'q':>3} {'pair':>10} {'c/x1':>8} {'2q':>5} {'overshoot':>10} {'x1':>8} {'c':>5}")
print("-" * 60)

cx1_best = {}
x1_best = {}

for q in sorted(all_data.keys()):
    sizes = sorted(all_data[q].keys())
    if len(sizes) < 2:
        continue

    # Use two largest sizes
    n1, n2 = sizes[-2], sizes[-1]
    e0_1, g1_1 = all_data[q][n1]
    e0_2, g1_2 = all_data[q][n2]

    vc = (e0_1/n1 - e0_2/n2) * 6 / (np.pi * (1/n2**2 - 1/n1**2))
    vx1_1 = g1_1 * n1 / (2 * np.pi)
    vx1_2 = g1_2 * n2 / (2 * np.pi)
    vx1_avg = (vx1_1 + vx1_2) / 2
    cx1 = vc / vx1_avg

    c = c_all[q]
    x1 = c / cx1
    overshoot = (cx1 - 2*q) / (2*q) * 100

    print(f"{q:>3} {f'{n1},{n2}':>10} {cx1:>8.4f} {2*q:>5} {overshoot:>9.1f}% {x1:>8.5f} {c:>5.2f}")
    cx1_best[q] = cx1
    x1_best[q] = x1

# ========== PART 3: Formula search for c/x1(q) ==========
print("\n" + "=" * 70)
print("PART 3: Formula Search for c/x₁(q)")
print("=" * 70)

qs = np.array(sorted(cx1_best.keys()), dtype=float)
cx1s = np.array([cx1_best[int(q)] for q in qs])

# Candidate 1: c/x1 = 2q
pred_2q = 2 * qs
rms_2q = np.sqrt(np.mean((cx1s - pred_2q)**2))
print(f"\n1. c/x₁ = 2q")
for q, m, p in zip(qs, cx1s, pred_2q):
    print(f"   q={int(q):>2}: measured={m:.3f}, predicted={p:.1f}, error={abs(m-p)/m*100:.1f}%")
print(f"   RMS = {rms_2q:.3f}")

# Candidate 2: c/x1 = aq + b
from numpy.polynomial import polynomial as P
coeffs_lin = np.polyfit(qs, cx1s, 1)
pred_lin = np.polyval(coeffs_lin, qs)
rms_lin = np.sqrt(np.mean((cx1s - pred_lin)**2))
print(f"\n2. c/x₁ = {coeffs_lin[0]:.4f}·q + {coeffs_lin[1]:.4f}")
for q, m, p in zip(qs, cx1s, pred_lin):
    print(f"   q={int(q):>2}: measured={m:.3f}, predicted={p:.3f}, error={abs(m-p)/m*100:.1f}%")
print(f"   RMS = {rms_lin:.3f}")

# Candidate 3: c/x1 = a*q^alpha
def power_law(q, a, alpha):
    return a * q**alpha
try:
    popt, pcov = curve_fit(power_law, qs, cx1s, p0=[2, 1])
    pred_pow = power_law(qs, *popt)
    rms_pow = np.sqrt(np.mean((cx1s - pred_pow)**2))
    print(f"\n3. c/x₁ = {popt[0]:.4f}·q^{popt[1]:.4f}")
    for q, m, p in zip(qs, cx1s, pred_pow):
        print(f"   q={int(q):>2}: measured={m:.3f}, predicted={p:.3f}, error={abs(m-p)/m*100:.1f}%")
    print(f"   RMS = {rms_pow:.3f}")
except:
    print("\n3. Power law fit failed")
    rms_pow = 999

# Candidate 4: c/x1 = 2q + a*ln(q)
def model_2q_log(q, a):
    return 2*q + a * np.log(q)
try:
    popt4, _ = curve_fit(model_2q_log, qs, cx1s, p0=[0.5])
    pred_4 = model_2q_log(qs, *popt4)
    rms_4 = np.sqrt(np.mean((cx1s - pred_4)**2))
    print(f"\n4. c/x₁ = 2q + {popt4[0]:.4f}·ln(q)")
    for q, m, p in zip(qs, cx1s, pred_4):
        print(f"   q={int(q):>2}: measured={m:.3f}, predicted={p:.3f}, error={abs(m-p)/m*100:.1f}%")
    print(f"   RMS = {rms_4:.3f}")
except:
    rms_4 = 999

# Candidate 5: c/x1 = q*(q+1)/... let's try c/x1 = q(q-1)/(q-2+epsilon)
# For q=2: 2*1/eps, q=3: 6/(1+eps), etc. Too many parameters.

# Candidate 6: c/x1 = 2q + a*(q-3)^+ (correction kicks in at q>=4)
def model_thresh(q, a, b):
    correction = np.where(q > 3.5, a * (q - 3)**b, 0)
    return 2*q + correction
try:
    # Only fit with q>=4 data
    popt6, _ = curve_fit(model_thresh, qs, cx1s, p0=[0.3, 0.7], maxfev=10000)
    pred_6 = model_thresh(qs, *popt6)
    rms_6 = np.sqrt(np.mean((cx1s - pred_6)**2))
    print(f"\n5. c/x₁ = 2q + {popt6[0]:.4f}·(q-3)^{popt6[1]:.4f} for q>3")
    for q, m, p in zip(qs, cx1s, pred_6):
        print(f"   q={int(q):>2}: measured={m:.3f}, predicted={p:.3f}, error={abs(m-p)/m*100:.1f}%")
    print(f"   RMS = {rms_6:.3f}")
except Exception as e:
    print(f"\n5. Threshold model failed: {e}")
    rms_6 = 999

# Candidate 7: FSS-corrected c/x1 — subtract systematic overshoot
# At q=2,3 we know exact c/x1. The overshoot is ~0.3% (q=2), ~0.3% (q=3).
# But at q=4, overshoot seems much larger. Is the intrinsic c/x1 exactly 2q?
# If correction ~ f(q)/N^2 and we use pair (N1,N2), the residual overshoot
# scales as 1/N2^2. Let's check.
print("\n\n--- Finite-size correction analysis ---")
print("All pairwise c/x₁ values to check convergence:")
for q in sorted(all_data.keys()):
    sizes = sorted(all_data[q].keys())
    pairs = []
    for i in range(len(sizes)-1):
        n1, n2 = sizes[i], sizes[i+1]
        e1, g1 = all_data[q][n1]
        e2, g2 = all_data[q][n2]
        vc = (e1/n1 - e2/n2) * 6 / (np.pi * (1/n2**2 - 1/n1**2))
        vx_1 = g1 * n1 / (2*np.pi)
        vx_2 = g2 * n2 / (2*np.pi)
        cx1 = vc / ((vx_1 + vx_2)/2)
        pairs.append((n1, n2, cx1))
    if pairs:
        pair_str = ", ".join([f"({n1},{n2}):{cx1:.4f}" for n1,n2,cx1 in pairs])
        print(f"  q={q}: {pair_str}  [2q={2*q}]")

# Extrapolate c/x1 to N->inf for q=4 (4 sizes)
print("\n--- Extrapolation of c/x₁ for q=4 (best data) ---")
q = 4
sizes = sorted(all_data[q].keys())
pair_cx1 = []
pair_inv = []
for i in range(len(sizes)-1):
    n1, n2 = sizes[i], sizes[i+1]
    e1, g1 = all_data[q][n1]
    e2, g2 = all_data[q][n2]
    vc = (e1/n1 - e2/n2) * 6 / (np.pi * (1/n2**2 - 1/n1**2))
    vx_1 = g1 * n1 / (2*np.pi)
    vx_2 = g2 * n2 / (2*np.pi)
    cx1 = vc / ((vx_1 + vx_2)/2)
    pair_cx1.append(cx1)
    pair_inv.append(1/n2**2)  # largest size in pair

pair_cx1 = np.array(pair_cx1)
pair_inv = np.array(pair_inv)

# Linear extrapolation in 1/N^2
if len(pair_cx1) >= 2:
    c_extrap = np.polyfit(pair_inv, pair_cx1, 1)
    cx1_inf = c_extrap[1]
    print(f"  Pairs: {list(zip([f'({sizes[i]},{sizes[i+1]})' for i in range(len(sizes)-1)], [f'{c:.4f}' for c in pair_cx1]))}")
    print(f"  Linear fit in 1/N²: c/x₁(∞) = {cx1_inf:.4f}")
    print(f"  Compared to 2q={2*q}: deviation = {(cx1_inf - 2*q)/2/q*100:.1f}%")

    # Also try 1/N extrapolation (logarithmic corrections → 1/ln(N) but approximate as 1/N)
    pair_invN = 1.0/np.array(sizes[1:])
    c_extrap2 = np.polyfit(pair_invN, pair_cx1, 1)
    cx1_inf2 = c_extrap2[1]
    print(f"  Linear fit in 1/N: c/x₁(∞) = {cx1_inf2:.4f}")
    print(f"  Compared to 2q={2*q}: deviation = {(cx1_inf2 - 2*q)/2/q*100:.1f}%")


# ========== PART 4: x₁(q) if c/x₁ = 2q ==========
print("\n" + "=" * 70)
print("PART 4: x₁(q) Under c/x₁ = 2q Hypothesis")
print("=" * 70)
print("If c/x₁ = 2q exactly, then x₁ = c/(2q):")
print(f"{'q':>3} {'c':>6} {'x₁=c/2q':>10} {'x₁ measured':>12} {'x₁ exact':>10}")
print("-" * 50)
for q_val in [2, 3, 4, 5, 7, 10]:
    c = c_all[q_val]
    x1_hyp = c / (2 * q_val)
    x1_meas = x1_best.get(q_val, None)
    exact_x1 = {2: 0.125, 3: 2/15}
    ex_str = f"{exact_x1[q_val]:.6f}" if q_val in exact_x1 else "—"
    m_str = f"{x1_meas:.6f}" if x1_meas else "—"
    print(f"{q_val:>3} {c:>6.3f} {x1_hyp:>10.6f} {m_str:>12} {ex_str:>10}")

print("\nWith c(q) ≈ 0.40·ln(q-1) + 0.55:")
print("  x₁(q) = [0.40·ln(q-1) + 0.55] / (2q)")
for q_val in [2, 3, 4, 5, 7, 10, 20, 50]:
    c_pred = 0.40 * np.log(q_val - 1) + 0.55
    x1_pred = c_pred / (2 * q_val)
    print(f"  q={q_val:>2}: c≈{c_pred:.3f}, x₁≈{x1_pred:.5f}")


# Save
results["cx1_best"] = {str(q): {"cx1": cx1_best[q], "x1": x1_best[q]} for q in cx1_best}
results["formulas"] = {
    "2q_rms": rms_2q,
    "linear_rms": rms_lin, "linear_coeffs": coeffs_lin.tolist(),
}
with open("results/sprint_058d_x1_formula.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved.")
