#!/usr/bin/env python3
"""Sprint 049d: Wide g scan for q=3 Potts to find the actual critical point.

S=0.047 at g=1.0 is confirmed by exact diag. The entropy is suspiciously low and
doesn't grow with n. Scan g from 0 to 2 to find where entropy peaks (= critical point).
Also do exact diag at n=8,10 for validation and faster exploration.
"""
import numpy as np, json, time, warnings
warnings.filterwarnings('ignore')

# === Exact diag for q=3 Potts ===
def exact_diag_potts(q, n, g, J=1.0):
    """Full exact diag for q-state Potts, n sites, OBC."""
    dim = q**n
    H = np.zeros((dim, dim))

    # Coupling: -J Σ_a P_a⊗P_a
    for i in range(n - 1):
        for s in range(dim):
            digits = [(s // q**k) % q for k in range(n)]
            if digits[i] == digits[i + 1]:
                H[s, s] += -J

    # Field: -g (X + X†)
    for i in range(n):
        for s in range(dim):
            digits = [(s // q**k) % q for k in range(n)]
            # X: shift state at site i by +1
            new_d = digits.copy()
            new_d[i] = (digits[i] + 1) % q
            s2 = sum(new_d[k] * q**k for k in range(n))
            H[s, s2] += -g
            # X†: shift by -1
            new_d[i] = (digits[i] - 1) % q
            s2 = sum(new_d[k] * q**k for k in range(n))
            H[s, s2] += -g

    from scipy.linalg import eigh
    evals, evecs = eigh(H)
    gs = evecs[:, 0]

    # Half-chain entropy
    gs_mat = gs.reshape(q**(n // 2), q**(n // 2))
    rho = gs_mat @ gs_mat.T
    rho_evals = np.linalg.eigvalsh(rho)
    rho_evals = rho_evals[rho_evals > 1e-15]
    S = -np.sum(rho_evals * np.log(rho_evals))

    return evals[0], S

# === Wide g scan with exact diag (fast, correct) ===
print("=== q=3 Potts: Exact diag, wide g scan ===", flush=True)
g_values = np.concatenate([
    np.arange(0.0, 0.10, 0.02),
    np.arange(0.10, 0.50, 0.05),
    np.arange(0.50, 2.01, 0.10),
])

results = {}
for n in [6, 8, 10]:
    print(f"\n  n={n} (dim={3**n}):", flush=True)
    results[n] = []
    t_start = time.time()
    for g in g_values:
        t0 = time.time()
        E, S = exact_diag_potts(3, n, float(g))
        dt = time.time() - t0
        print(f"    g={g:.2f}: S={S:.6f}, E={E:.6f}, t={dt:.1f}s", flush=True)
        results[n].append({'g': float(g), 'E': float(E), 'S_half': float(S), 'time': dt})
        if dt > 30 and n == 10:
            print(f"    Too slow at n={n}, stopping", flush=True)
            break
    print(f"  Total time: {time.time() - t_start:.0f}s", flush=True)

# Save
with open('results/sprint_049d_potts_wide.json', 'w') as f:
    json.dump({'sprint': '049d', 'q': 3, 'method': 'exact_diag',
               'data': {str(k): v for k, v in results.items()}}, f, indent=2)

# === Find entropy peak ===
print("\n=== Entropy peak analysis ===", flush=True)
for n, data in results.items():
    gs = [d['g'] for d in data]
    Ss = [d['S_half'] for d in data]
    peak_idx = np.argmax(Ss)
    print(f"  n={n}: S_max = {Ss[peak_idx]:.6f} at g = {gs[peak_idx]:.2f}", flush=True)

    # dS/dg
    gs_a = np.array(gs)
    Ss_a = np.array(Ss)
    dSdg = np.gradient(Ss_a, gs_a)
    # Where does dS/dg change sign?
    for i in range(len(dSdg) - 1):
        if dSdg[i] > 0 and dSdg[i + 1] < 0:
            g_peak = gs_a[i] + (gs_a[i+1]-gs_a[i]) * dSdg[i] / (dSdg[i] - dSdg[i+1])
            print(f"         dS/dg sign change at g ≈ {g_peak:.3f}", flush=True)

# === Also scan q=2 TFIM with exact diag for comparison ===
print("\n=== q=2 TFIM: Exact diag comparison ===", flush=True)
for n in [8, 10]:
    data = []
    for g in [0.5, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.5]:
        E, S = exact_diag_potts(2, n, g)
        print(f"  n={n} g={g:.2f}: S={S:.4f}", flush=True)
        data.append({'g': g, 'S': S})

# === q=4 Potts exact diag ===
print("\n=== q=4 Potts: Exact diag ===", flush=True)
for n in [4, 6]:
    print(f"  n={n} (dim={4**n}):", flush=True)
    for g in np.arange(0.0, 2.01, 0.20):
        t0 = time.time()
        E, S = exact_diag_potts(4, n, float(g))
        dt = time.time() - t0
        print(f"    g={g:.2f}: S={S:.6f}, t={dt:.1f}s", flush=True)
        if dt > 30:
            print(f"    Too slow, stopping", flush=True)
            break

print("\nDone!", flush=True)
