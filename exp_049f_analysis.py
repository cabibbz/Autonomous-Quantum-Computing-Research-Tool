#!/usr/bin/env python3
"""Sprint 049f: Full analysis of all Sprint 049 data.

1. q=2 TFIM: entropy FSS → ν, central charge c
2. q=2 TFIM: correlation length → ν
3. q=3 Potts: locate true g_c from entropy scaling
4. Compare and synthesize
"""
import numpy as np, json

# ================================================================
# 1. q=2 TFIM Entropy FSS
# ================================================================
print("=" * 60)
print("1. q=2 TFIM: Entropy FSS")
print("=" * 60)

with open('results/sprint_049b_entropy_fss.json') as f:
    ent_data = json.load(f)

q2 = ent_data['q2']
sizes_q2 = sorted([int(k) for k in q2.keys()])

# Table of S at key g values
print("\nS(n, g) table:")
key_gs = [0.90, 0.95, 0.97, 0.99, 1.00, 1.01, 1.03, 1.05, 1.10]
header = "  g    " + "  ".join(f"n={n:3d}" for n in sizes_q2)
print(header)
for g in key_gs:
    row = f"  {g:.2f}  "
    for n in sizes_q2:
        data = q2[str(n)]
        match = [d for d in data if abs(d['g'] - g) < 0.001]
        if match:
            row += f" {match[0]['S_half']:.4f}"
        else:
            row += "    -- "
    print(row)

# Central charge from S at g_c=1.0
print("\nCentral charge c from S(g_c=1.0) vs ln(n):")
ns, Ss = [], []
for n in sizes_q2:
    data = q2[str(n)]
    match = [d for d in data if abs(d['g'] - 1.0) < 0.001]
    if match:
        ns.append(n)
        Ss.append(match[0]['S_half'])
        print(f"  n={n:3d}: S = {match[0]['S_half']:.6f}")

ns_a = np.array(ns, dtype=float)
Ss_a = np.array(Ss)
coeffs = np.polyfit(np.log(ns_a), Ss_a, 1)
c_fit = 6 * coeffs[0]
print(f"\n  c = {c_fit:.4f} (exact Ising c = 0.5000)")
print(f"  Fit: S = {coeffs[0]:.4f} * ln(n) + {coeffs[1]:.4f}")

# Also fit c from pairs
print("\n  Pairwise c estimates:")
for i in range(len(ns) - 1):
    c_pair = 6 * (Ss_a[i+1] - Ss_a[i]) / (np.log(ns_a[i+1]) - np.log(ns_a[i]))
    print(f"    n={ns[i]},{ns[i+1]}: c = {c_pair:.4f}")

# dS/dg peak location (pseudo-critical point)
print("\nPseudo-critical point from dS/dg peak:")
for n in sizes_q2:
    data = q2[str(n)]
    gs = np.array([d['g'] for d in data])
    Ss = np.array([d['S_half'] for d in data])
    order = np.argsort(gs)
    gs, Ss = gs[order], Ss[order]
    dSdg = np.gradient(Ss, gs)
    peak = np.argmin(dSdg)  # Most negative slope
    print(f"  n={n:3d}: g_pseudo = {gs[peak]:.3f}, dS/dg = {dSdg[peak]:.3f}")

# FSS collapse for ν
print("\nEntropy FSS collapse quality:")
g_all, S_all, n_all = [], [], []
for n in sizes_q2:
    for d in q2[str(n)]:
        g_all.append(d['g']); S_all.append(d['S_half']); n_all.append(n)
g_all = np.array(g_all); S_all = np.array(S_all); n_all = np.array(n_all)

def collapse_quality(nu, g_c, g_arr, y_arr, n_arr):
    sizes = sorted(set(n_arr))
    curves = {}
    for s in sizes:
        mask = (n_arr == s)
        x = (g_arr[mask] - g_c) * s**(1.0/nu)
        y = y_arr[mask]
        order = np.argsort(x)
        curves[s] = (x[order], y[order])
    total_err, count = 0.0, 0
    for i, s1 in enumerate(sizes):
        for s2 in sizes[i+1:]:
            x1, y1 = curves[s1]
            x2, y2 = curves[s2]
            x_min = max(x1.min(), x2.min())
            x_max = min(x1.max(), x2.max())
            if x_min >= x_max: continue
            xg = np.linspace(x_min, x_max, 50)
            y1i = np.interp(xg, x1, y1)
            y2i = np.interp(xg, x2, y2)
            diff = (y1i - y2i)**2
            norm = (y1i**2 + y2i**2) / 2
            total_err += np.sum(diff / (norm + 1e-15))
            count += len(xg)
    return total_err / max(count, 1)

for nu, label in [(0.5, 'mean-field'), (0.8, '≈Potts q=3'), (1.0, 'Ising exact'),
                   (1.2, ''), (1.5, ''), (2.0, '')]:
    q_val = collapse_quality(nu, 1.0, g_all, S_all, n_all)
    print(f"  ν={nu:.1f} ({label:14s}): quality = {q_val:.6f}")

nu_scan = np.linspace(0.3, 3.0, 55)
qualities = [collapse_quality(nu, 1.0, g_all, S_all, n_all) for nu in nu_scan]
best = nu_scan[np.argmin(qualities)]
print(f"  Best ν = {best:.2f} (quality = {min(qualities):.6f})")

# ================================================================
# 2. q=2 TFIM Correlation Length
# ================================================================
print("\n" + "=" * 60)
print("2. q=2 TFIM: Correlation Length ξ")
print("=" * 60)

with open('results/sprint_049a_tfim_xi.json') as f:
    xi_data = json.load(f)

for n_str, data in xi_data['sizes'].items():
    n = int(n_str)
    print(f"\n  n={n}:")
    ok = [d for d in data if d['status'] == 'ok']
    for d in ok:
        xi_str = f"{d['xi']:.2f}" if d['xi'] < 1e5 else ">>n"
        print(f"    g={d['g']:.2f}: ξ={xi_str:>8s}, R²={d['R2']:.3f}, S={d['S_half']:.4f}")

# ν from ξ on disordered side
print("\nν extraction from ξ(g) on disordered side (g > 1):")
for n_str, data in xi_data['sizes'].items():
    n = int(n_str)
    ok = [d for d in data if d['status'] == 'ok' and d['g'] > 1.01
          and d['xi'] < n/2 and d['R2'] > 0.9 and d['xi'] > 0.5]
    if len(ok) < 3:
        print(f"  n={n}: insufficient points ({len(ok)})")
        continue
    gs = np.array([d['g'] for d in ok])
    xis = np.array([d['xi'] for d in ok])
    log_dg = np.log(gs - 1.0)
    log_xi = np.log(xis)
    coeffs = np.polyfit(log_dg, log_xi, 1)
    nu = -coeffs[0]
    print(f"  n={n}: ν = {nu:.3f} from {len(ok)} points (g ∈ [{gs[0]:.2f}, {gs[-1]:.2f}])")

    # Comparison with exact ξ = 1/ln(g)
    xi_exact = 1.0 / np.log(gs)
    ratios = xis / xi_exact
    print(f"         ξ/ξ_exact: {', '.join(f'{r:.3f}' for r in ratios)}")

# ================================================================
# 3. q=3 Potts: True Critical Point
# ================================================================
print("\n" + "=" * 60)
print("3. q=3 Potts: True Critical Point")
print("=" * 60)

with open('results/sprint_049e_potts_critical.json') as f:
    potts_data = json.load(f)

sizes_q3 = sorted([int(k) for k in potts_data['data'].keys()])
print(f"Available sizes: {sizes_q3}")

# Entropy table
print("\nS(n, g) table for q=3 Potts:")
g_focus = [0.15, 0.20, 0.24, 0.26, 0.28, 0.30, 0.32, 0.36, 0.40, 0.50, 0.80]
for g in g_focus:
    row = f"  g={g:.2f}: "
    for n in sizes_q3:
        data = potts_data['data'][str(n)]
        match = [d for d in data if abs(d['g'] - g) < 0.005]
        if match:
            row += f"S({n})={match[0]['S_half']:.4f}  "
    print(row)

# Central charge estimates from pairs
print("\nCentral charge c estimates from S(n1) vs S(n2):")
if len(sizes_q3) >= 2:
    n1, n2 = sizes_q3[0], sizes_q3[1]
    for g in g_focus:
        S_vals = {}
        for n in [n1, n2]:
            data = potts_data['data'][str(n)]
            match = [d for d in data if abs(d['g'] - g) < 0.005]
            if match:
                S_vals[n] = match[0]['S_half']
        if len(S_vals) == 2:
            c = 6 * (S_vals[n2] - S_vals[n1]) / (np.log(n2) - np.log(n1))
            print(f"  g={g:.2f}: c({n1},{n2}) = {c:.3f}", end="")
            if abs(c - 0.8) < 0.2:
                print(f"  ← near c=4/5 for q=3 Potts!")
            else:
                print()

# Identify g_c: where c is closest to 4/5
print(f"\n  q=3 Potts exact central charge: c = 4/5 = 0.800")

# ================================================================
# 4. Implications for Prior Results
# ================================================================
print("\n" + "=" * 60)
print("4. Implications for Prior Potts Results")
print("=" * 60)

print("""
Prior sprints (038-048) assumed g_c ≈ 1.0 for q=3 Potts based on MI-CV crossings.
Exact diag + DMRG now show:
  - At g=1.0: S = 0.047 (independent of n → gapped, disordered phase)
  - At g=0.15: S = ln(3) = 1.099 (GHZ/ordered phase)
  - Entropy drops sharply at g ≈ 0.25-0.35
  - S grows with n up to g ≈ 0.35 → critical region

The MI-CV crossings near g ≈ 0.9-1.0 were measuring a CROSSOVER in the
disordered phase, not the actual phase transition.

Impact on prior results:
  - g_c scaling law g_c(q) = 0.87*(q-3)^{-0.85} may be wrong
  - ν(q) extraction for q≥3 needs revision at true g_c values
  - The qualitative MI-CV signatures (crossings) are real but not at g_c

Possible cause: the Hamiltonian normalization differs from the self-dual form.
Our coupling: -J Σ_a P_a⊗P_a,  field: -g(X+X†) = -g Σ_{k=1}^{q-1} τ^k
The self-dual point depends on q: g_c = f(q, J), NOT simply g_c = J.
""")

print("Analysis complete!")
