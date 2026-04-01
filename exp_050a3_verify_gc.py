#!/usr/bin/env python3
"""Sprint 050a3: Verify self-duality prediction g_c = J/q for our Potts Hamiltonian.

Our PottsChain: H = -J Σ δ(s_i,s_j) - g Σ (X + X†)
For q=2,3: X+X† = Σ_{k=1}^{q-1} X^k, so self-duality applies → g_c = J/q.
For q≥4: X+X† ≠ Σ_{k=1}^{q-1} X^k, self-duality is broken.

Test 1: q=2 entropy peak should be at g = J/2 = 0.5 (NOT 1.0 — that's TFIChain convention)
Test 2: q=3 entropy peak should be at g = J/3 ≈ 0.333
Test 3: q=4 — find where entropy peaks (no self-dual prediction)
"""
import numpy as np, json, time
from scipy.linalg import eigh


def exact_diag_potts(q, n, g, J=1.0):
    dim = q**n
    H = np.zeros((dim, dim))
    for i in range(n - 1):
        for s in range(dim):
            digits = [(s // q**k) % q for k in range(n)]
            if digits[i] == digits[i + 1]:
                H[s, s] += -J
    for i in range(n):
        for s in range(dim):
            digits = [(s // q**k) % q for k in range(n)]
            new_d = digits.copy()
            new_d[i] = (digits[i] + 1) % q
            s2 = sum(new_d[k] * q**k for k in range(n))
            H[s, s2] += -g
            new_d[i] = (digits[i] - 1) % q
            s2 = sum(new_d[k] * q**k for k in range(n))
            H[s, s2] += -g
    evals, evecs = eigh(H)
    gs_vec = evecs[:, 0]
    gs_mat = gs_vec.reshape(q**(n // 2), q**(n // 2))
    rho = gs_mat @ gs_mat.T
    rho_evals = np.linalg.eigvalsh(rho)
    rho_evals = rho_evals[rho_evals > 1e-15]
    S = -np.sum(rho_evals * np.log(rho_evals))
    return evals[0], S


# === Test 1: q=2, n=8 — where does entropy peak? ===
print("=== q=2 Potts (exact diag n=8, dim=256) ===", flush=True)
print("Self-dual prediction: g_c = J/q = 0.500", flush=True)
g_vals_q2 = np.arange(0.05, 1.51, 0.05)
s_q2 = []
for g in g_vals_q2:
    E, S = exact_diag_potts(2, 8, float(g))
    print(f"  g={g:.2f}: S={S:.6f}", flush=True)
    s_q2.append(S)
peak_idx = np.argmax(s_q2)
print(f"  *** Entropy PEAK at g={g_vals_q2[peak_idx]:.2f}, S={s_q2[peak_idx]:.6f} ***", flush=True)

# === Test 2: q=3, n=8 — entropy peak ===
# Already have this data from exp_050a2, but let me check dS/dg more carefully
print("\n=== q=3 Potts (exact diag n=8, dim=6561) ===", flush=True)
print("Self-dual prediction: g_c = J/q = 0.333", flush=True)
g_vals_q3 = np.arange(0.10, 0.60, 0.02)
s_q3 = []
for g in g_vals_q3:
    E, S = exact_diag_potts(3, 8, float(g))
    print(f"  g={g:.2f}: S={S:.6f}", flush=True)
    s_q3.append(S)

# Find where derivative is most negative (transition)
dSdg = np.gradient(s_q3, g_vals_q3)
min_idx = np.argmin(dSdg)
print(f"  *** Steepest drop at g={g_vals_q3[min_idx]:.2f}, dS/dg={dSdg[min_idx]:.3f} ***", flush=True)

# Find inflection (second derivative peak)
d2Sdg2 = np.gradient(dSdg, g_vals_q3)
infl_idx = np.argmax(np.abs(d2Sdg2))
print(f"  *** Inflection at g={g_vals_q3[infl_idx]:.2f} ***", flush=True)

# === Also do q=2, q=3 at n=6 for size comparison ===
print("\n=== q=2 n=6 for comparison ===", flush=True)
s_q2_n6 = []
for g in g_vals_q2:
    E, S = exact_diag_potts(2, 6, float(g))
    s_q2_n6.append(S)
peak_idx_6 = np.argmax(s_q2_n6)
print(f"  Entropy peak at g={g_vals_q2[peak_idx_6]:.2f}", flush=True)

print(f"\n=== Comparison: q=2 pseudo-critical point drift ===", flush=True)
dSdg_q2_n6 = np.gradient(s_q2_n6, g_vals_q2)
dSdg_q2_n8 = np.gradient(s_q2, g_vals_q2)
min6 = np.argmin(dSdg_q2_n6)
min8 = np.argmin(dSdg_q2_n8)
print(f"  n=6 steepest drop: g={g_vals_q2[min6]:.2f}", flush=True)
print(f"  n=8 steepest drop: g={g_vals_q2[min8]:.2f}", flush=True)
print(f"  Self-dual g_c = 0.500", flush=True)

# === Test: q=4, n=6 ===
print("\n=== q=4 Potts (exact diag n=6, dim=4096) ===", flush=True)
print("No self-dual prediction (X+X† ≠ Σ X^k for q≥4)", flush=True)
g_vals_q4 = np.arange(0.05, 1.01, 0.05)
s_q4 = []
for g in g_vals_q4:
    E, S = exact_diag_potts(4, 6, float(g))
    print(f"  g={g:.2f}: S={S:.6f}", flush=True)
    s_q4.append(S)

dSdg_q4 = np.gradient(s_q4, g_vals_q4)
min_q4 = np.argmin(dSdg_q4)
print(f"  *** Steepest drop at g={g_vals_q4[min_q4]:.2f} ***", flush=True)

# === Test: q=5, n=6 ===
print("\n=== q=5 Potts (exact diag n=6, dim=15625) ===", flush=True)
g_vals_q5 = np.arange(0.05, 0.51, 0.05)
s_q5 = []
for g in g_vals_q5:
    t0 = time.time()
    E, S = exact_diag_potts(5, 6, float(g))
    dt = time.time() - t0
    print(f"  g={g:.2f}: S={S:.6f}, t={dt:.1f}s", flush=True)
    s_q5.append(S)
    if dt > 30:
        print("  Too slow, stopping", flush=True)
        break

if len(s_q5) > 3:
    dSdg_q5 = np.gradient(s_q5, g_vals_q5[:len(s_q5)])
    min_q5 = np.argmin(dSdg_q5)
    print(f"  *** Steepest drop at g={g_vals_q5[min_q5]:.2f} ***", flush=True)

# Save all results
results = {
    'q2_n8': [{'g': float(g), 'S': float(s)} for g, s in zip(g_vals_q2, s_q2)],
    'q3_n8': [{'g': float(g), 'S': float(s)} for g, s in zip(g_vals_q3, s_q3)],
    'q4_n6': [{'g': float(g), 'S': float(s)} for g, s in zip(g_vals_q4, s_q4)],
    'q5_n6': [{'g': float(g), 'S': float(s)} for g, s in zip(g_vals_q5[:len(s_q5)], s_q5)],
}
with open('results/sprint_050a3_gc_verify.json', 'w') as f:
    json.dump({'sprint': '050a3', 'method': 'exact_diag', 'data': results}, f, indent=2)

print("\nDone!", flush=True)
