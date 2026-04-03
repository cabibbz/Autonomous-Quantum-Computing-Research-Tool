#!/usr/bin/env python3
"""Sprint 054i: Final central charge analysis with correct q=2 data."""
import numpy as np, json

# === Correct data ===
data = {
    2: {'g_c': 1.0, 'sizes': [16, 24, 32, 48, 64],
        'S': [0.423409, 0.460207, 0.485718, 0.521066, 0.545827],
        'c_exact': 0.500, 'source': 'Sprint 049 (TFIChain)'},
    3: {'g_c': 1/3, 'sizes': [8, 12, 16, 24],
        'S': [0.516117, 0.579230, 0.622697, 0.682416],
        'c_exact': 0.800, 'source': 'Sprint 054a (PottsChain, chi=80)'},
    4: {'g_c': 0.392, 'sizes': [8, 12, 16, 24],
        'S': [0.646366, 0.729571, 0.788153, 0.871215],
        'c_exact': 1.000, 'source': 'Sprint 054f (PottsChain, chi=40)'},
    5: {'g_c': 0.441, 'sizes': [8, 12, 16],
        'S': [0.760529, 0.854768, 0.918754],
        'c_exact': None, 'source': 'Sprint 054f (PottsChain, chi=40)'},
}

print("=" * 70)
print("CENTRAL CHARGE c(q) AT TRUE CRITICAL POINTS")
print("=" * 70)

results = {}
for q, d in data.items():
    ns = np.array(d['sizes'], dtype=float)
    Ss = np.array(d['S'], dtype=float)
    ln_n = np.log(ns)

    A = np.vstack([ln_n, np.ones(len(ln_n))]).T
    slope, intercept = np.linalg.lstsq(A, Ss, rcond=None)[0]
    c_full = 6 * slope

    c_pairs = []
    for i in range(len(ns) - 1):
        c_pair = 6 * (Ss[i+1] - Ss[i]) / (ln_n[i+1] - ln_n[i])
        c_pairs.append(float(c_pair))

    print(f"\nq={q} (g_c={d['g_c']:.4f}):")
    print(f"  Full fit c = {c_full:.4f}", end="")
    if d['c_exact']:
        print(f"  (CFT: {d['c_exact']:.3f})")
    else:
        print("  (no CFT prediction)")
    for i, cp in enumerate(c_pairs):
        n1, n2 = d['sizes'][i], d['sizes'][i+1]
        s = f"  c({n1},{n2}) = {cp:.4f}"
        if d['c_exact']:
            s += f"  (overshoot: {cp/d['c_exact']*100-100:+.1f}%)"
        print(s)

    results[str(q)] = {
        'g_c': d['g_c'], 'c_full': float(c_full),
        'c_pairwise': c_pairs, 'c_exact': d['c_exact'],
        'sizes': d['sizes'], 'S': d['S'],
    }

# === Summary table ===
print("\n" + "=" * 70)
print("SUMMARY TABLE")
print("=" * 70)
print(f"{'q':>3} {'g_c':>7} {'c(raw)':>8} {'c(CFT)':>8} {'Overshoot':>10} {'Trend':>12}")
print("-" * 54)
for q in [2, 3, 4, 5]:
    r = results[str(q)]
    cp_last = r['c_pairwise'][-1]
    c_exact = r['c_exact']
    if c_exact:
        overshoot = f"{cp_last/c_exact*100-100:+.1f}%"
    else:
        overshoot = "N/A"
    # Check if converging
    if len(r['c_pairwise']) >= 2:
        diff = r['c_pairwise'][-1] - r['c_pairwise'][-2]
        trend = f"{'↓' if diff < 0 else '↑'} {abs(diff):.4f}"
    else:
        trend = "—"
    print(f"{q:3d} {r['g_c']:7.4f} {cp_last:8.4f} {str(c_exact or '?'):>8} {overshoot:>10} {trend:>12}")

# === Chi convergence ===
print("\n" + "=" * 70)
print("CHI CONVERGENCE (q=3, n=16)")
print("=" * 70)
print("  chi=20: S=0.622695")
print("  chi=40: S=0.622697")
print("  chi=60: S=0.622697 (from 054a)")
print("  chi=80: S=0.622697 (from 054a)")
print("  Conclusion: entropy fully converged at chi=20. Chi is NOT the overshoot source.")

# === Exact diag cross-check ===
print("\n" + "=" * 70)
print("EXACT DIAG CROSS-CHECK (q=3, n=8)")
print("=" * 70)
print("  Exact: E=-8.64004396, S=0.516117")
print("  DMRG:  E=-8.64004396, S=0.516117")
print("  ΔE = 2.1e-14, ΔS = 0.000000")
print("  DMRG agrees with exact diag to machine precision.")

# === c(q) growth pattern ===
print("\n" + "=" * 70)
print("c(q) GROWTH PATTERN")
print("=" * 70)
print("Raw largest-pair estimates:")
qs = [2, 3, 4, 5]
cs = [results[str(q)]['c_pairwise'][-1] for q in qs]
for q, c in zip(qs, cs):
    print(f"  q={q}: c ≈ {c:.4f}")

# If overshoot ratio is ~(1 + A/(q-1)) at these sizes:
print("\nOvershoot-corrected estimates (using q=2,3 calibration):")
# q=2 at (n=48,64): +3.4%. q=3 at (n=16,24): +10.5%
# If overshoot ~ B*q at (n=16,24):
# q=2 at (n=16,24): c(16,24)=0.5462, overshoot = +9.2%
q2_os_16_24 = 6*(0.460207-0.423409)/(np.log(24)-np.log(16))
print(f"  q=2 c(16,24) = {q2_os_16_24:.4f}, overshoot = {q2_os_16_24/0.5*100-100:+.1f}%")
print(f"  q=3 c(16,24) = 0.8837, overshoot = +10.5%")
print(f"  q=4 c(16,24) = 1.2291, overshoot = +22.9%")
print(f"  q=5 c(12,16) = 1.3345 (smaller sizes, more overshoot expected)")

# Save
with open('results/sprint_054i_final_analysis.json', 'w') as f:
    json.dump({'sprint': '054i', 'results': results,
               'chi_convergence': 'entropy converged at chi=20',
               'exact_diag_check': 'DMRG matches to machine precision',
               'conclusion': 'c(q) increases with q, overshoot increases with q'}, f, indent=2)

print("\nDone.", flush=True)
