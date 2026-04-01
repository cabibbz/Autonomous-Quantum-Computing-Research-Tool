#!/usr/bin/env python3
"""Sprint 050a2: Find q=3 Potts g_c using exact diag (n=6,8) + Sprint 049e DMRG (n=12).

DMRG at q=3 without symmetry conservation is too slow (353s/point at n=16 chi=40).
Use exact diag for n=6,8 (dim 729, 6561) and existing n=12 DMRG data.
Dense g scan near g=0.20-0.35 where the transition lives.
"""
import numpy as np, json, time, warnings
warnings.filterwarnings('ignore')
from scipy.linalg import eigh


def exact_diag_potts(q, n, g, J=1.0):
    """Full exact diag for q-state Potts."""
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


# Dense g values in the critical region
g_values = [0.05, 0.10, 0.14, 0.16, 0.18, 0.20, 0.21, 0.22, 0.23, 0.24,
            0.25, 0.26, 0.27, 0.28, 0.30, 0.32, 0.34, 0.36, 0.40, 0.50, 0.70, 1.00]

results = {}

# Exact diag for n=6 and n=8
for n in [6, 8]:
    print(f"\n=== q=3 Potts n={n} (exact diag, dim={3**n}) ===", flush=True)
    results[n] = []
    for g in g_values:
        t0 = time.time()
        E, S = exact_diag_potts(3, n, float(g))
        dt = time.time() - t0
        print(f"  g={g:.2f}: S={S:.6f}, E/n={E/n:.4f}, t={dt:.1f}s", flush=True)
        results[n].append({'g': float(g), 'E': float(E), 'S_half': float(S), 'time': dt})

# Save
with open('results/sprint_050a_q3_gc.json', 'w') as f:
    json.dump({'sprint': '050a', 'q': 3, 'method': 'exact_diag + DMRG',
               'data': {str(k): v for k, v in results.items()}}, f, indent=2)

# Load Sprint 049e DMRG data for n=12
try:
    with open('results/sprint_049e_potts_critical.json') as f:
        old = json.load(f)
    n12_data = {round(d['g'], 2): d['S_half'] for d in old['data'].get('12', [])}
    n8_dmrg = {round(d['g'], 2): d['S_half'] for d in old['data'].get('8', [])}
    print(f"\nLoaded n=12 DMRG data: {len(n12_data)} points", flush=True)
except:
    n12_data = {}
    n8_dmrg = {}

# Validate exact diag vs DMRG at n=8
print("\n=== Validation: exact diag vs DMRG at n=8 ===", flush=True)
n8_exact = {round(d['g'], 2): d['S_half'] for d in results.get(8, [])}
for g in sorted(set(n8_exact.keys()) & set(n8_dmrg.keys())):
    diff = abs(n8_exact[g] - n8_dmrg[g])
    print(f"  g={g:.2f}: exact={n8_exact[g]:.6f}, DMRG={n8_dmrg[g]:.6f}, diff={diff:.2e}", flush=True)

# === Central charge analysis ===
print("\n=== Central charge c(g) from S(n) = (c/6)ln(n) + const ===", flush=True)
all_data = {}
# n=6 exact
for d in results.get(6, []):
    g_r = round(d['g'], 2)
    if g_r not in all_data:
        all_data[g_r] = {}
    all_data[g_r][6] = d['S_half']
# n=8 exact
for d in results.get(8, []):
    g_r = round(d['g'], 2)
    if g_r not in all_data:
        all_data[g_r] = {}
    all_data[g_r][8] = d['S_half']
# n=12 DMRG
for g, s in n12_data.items():
    if g in all_data:
        all_data[g][12] = s

# Also add n=16 timing test point
# S=1.109430 at g=0.26, n=16, chi=60
all_data.setdefault(0.26, {})[16] = 1.109430
# And n=16 chi=40 point
all_data.setdefault(0.18, {})[16] = 1.101584

print(f"\n{'g':>6s} | {'c(all)':>7s} | {'c(6,8)':>7s} | {'c(8,12)':>7s} | S(6)    S(8)    S(12)", flush=True)
print("-" * 80, flush=True)

c_results = []
for g in sorted(all_data.keys()):
    ns = sorted(all_data[g].keys())
    Ss = {n: all_data[g][n] for n in ns}

    # Full fit if >=3 sizes
    c_all = None
    if len(ns) >= 3:
        ns_a = np.array(ns, dtype=float)
        Ss_a = np.array([Ss[n] for n in ns])
        coeffs = np.polyfit(np.log(ns_a), Ss_a, 1)
        c_all = 6 * coeffs[0]

    # Pairwise
    c_68 = None
    if 6 in Ss and 8 in Ss:
        c_68 = 6 * (Ss[8] - Ss[6]) / (np.log(8) - np.log(6))
    c_812 = None
    if 8 in Ss and 12 in Ss:
        c_812 = 6 * (Ss[12] - Ss[8]) / (np.log(12) - np.log(8))

    s_str = "  ".join([f"S({n})={Ss[n]:.4f}" for n in ns])
    c_all_str = f"{c_all:.3f}" if c_all is not None else "  --  "
    c_68_str = f"{c_68:.3f}" if c_68 is not None else "  --  "
    c_812_str = f"{c_812:.3f}" if c_812 is not None else "  --  "
    print(f"  {g:.2f} | {c_all_str:>7s} | {c_68_str:>7s} | {c_812_str:>7s} | {s_str}", flush=True)

    if c_all is not None:
        c_results.append((g, c_all, c_68, c_812))

# Find g_c where c is closest to 4/5
if c_results:
    print(f"\n=== Finding g_c: where c → 4/5 = 0.800 ===", flush=True)
    # Use c(8,12) which is the larger-n pairwise (more reliable)
    c812_vals = [(g, c812) for g, _, _, c812 in c_results if c812 is not None]
    if c812_vals:
        best_g, best_c = min(c812_vals, key=lambda x: abs(x[1] - 0.800))
        print(f"  c(8,12) closest to 0.800: g={best_g:.2f}, c={best_c:.3f}", flush=True)

    # Interpolate: find g where c(8,12) = 0.800
    for i in range(len(c812_vals) - 1):
        g1, c1 = c812_vals[i]
        g2, c2 = c812_vals[i + 1]
        if (c1 - 0.800) * (c2 - 0.800) < 0:
            g_interp = g1 + (g2 - g1) * (c1 - 0.800) / (c1 - c2)
            print(f"  c(8,12) = 0.800 interpolated at g_c ≈ {g_interp:.4f}", flush=True)

    # Also find dS/dg peak (pseudo-critical point)
    print(f"\n=== dS/dg peak (pseudo-critical point) ===", flush=True)
    for n in [6, 8]:
        gs = [d['g'] for d in results[n]]
        Ss = [d['S_half'] for d in results[n]]
        dSdg = np.gradient(Ss, gs)
        # Most negative dS/dg (entropy drops fastest)
        min_idx = np.argmin(dSdg)
        print(f"  n={n}: max |dS/dg| at g={gs[min_idx]:.2f}, dS/dg={dSdg[min_idx]:.3f}", flush=True)
        # Steepest decline is on the disordered side of g_c

print("\nDone!", flush=True)
