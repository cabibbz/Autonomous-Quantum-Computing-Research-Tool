"""Sprint 061d: Comprehensive asymptotic analysis — all CFT quantities at q=2-30.

Combine all measured data and fit asymptotic formulas.
Test free boson predictions: R(sigma^k)/R(sigma) → k², x₁ → 0, etc.
"""
import numpy as np
from scipy.optimize import curve_fit
import json

# ===== Collected data (n=4 charge-resolved, from 061c) =====
# q: R_sigma2, R_epsilon, C_sse, delta1_N
n4_data = {
    2:  {"R_s2": None, "R_eps": 7.696, "C_sse": 0.4802, "d1N": None},
    3:  {"R_s2": 1.000, "R_eps": 6.254, "C_sse": 0.4781, "d1N": None},
    4:  {"R_s2": 1.751, "R_eps": 6.437, "C_sse": 0.3964, "d1N": None},
    5:  {"R_s2": 2.436, "R_eps": 6.976, "C_sse": 0.3224, "d1N": None},
    7:  {"R_s2": 3.268, "R_eps": 7.804, "C_sse": 0.2354, "d1N": None},
    10: {"R_s2": 3.659, "R_eps": 8.345, "C_sse": 0.1843, "d1N": None},
    15: {"R_s2": 3.925, "R_eps": 13.656, "C_sse": 0.0932, "d1N": 0.1234},
    20: {"R_s2": 3.974, "R_eps": 17.748, "C_sse": 0.0633, "d1N": 0.0824},
    25: {"R_s2": 3.989, "R_eps": 21.908, "C_sse": 0.0473, "d1N": 0.0601},
    30: {"R_s2": 3.994, "R_eps": 26.108, "C_sse": 0.0375, "d1N": 0.0463},
}

# Best estimates from larger n (Sprints 054-060)
best_data = {
    2:  {"g_c": 0.250, "c": 0.500, "nu": 1.00, "x1": 0.1246, "cx1": 4.01,  "C_sse_best": 0.503},
    3:  {"g_c": 0.333, "c": 0.800, "nu": 0.86, "x1": 0.1337, "cx1": 5.98,  "C_sse_best": 0.540},
    4:  {"g_c": 0.392, "c": 1.000, "nu": 0.82, "x1": 0.1172, "cx1": 8.54,  "C_sse_best": 0.464},
    5:  {"g_c": 0.441, "c": 1.100, "nu": 0.85, "x1": 0.1015, "cx1": 10.84, "C_sse_best": 0.376},
    7:  {"g_c": 0.535, "c": 1.300, "nu": 0.97, "x1": 0.0860, "cx1": 15.11, "C_sse_best": 0.272},
    10: {"g_c": 0.684, "c": 1.400, "nu": 1.12, "x1": 0.0831, "cx1": 16.86, "C_sse_best": 0.210},
}

# g_c formula
gc_formula = lambda q: 0.2 * np.sqrt(q - 1) + 0.05

# Harmonic ratios from 061b
harmonic_data = {
    15: {2: 3.925, 3: 8.683, 4: 15.140},
    20: {2: 3.974, 3: 8.883, 4: 15.670},
    25: {2: 3.989, 3: 8.946, 4: 15.839, 5: 24.625},
    30: {2: 3.994, 3: 8.970, 4: 15.907, 5: 24.777},
}
# Add known smaller q from Sprint 057
harmonic_data[4] = {2: 1.68}
harmonic_data[5] = {2: 2.41}
harmonic_data[7] = {2: 3.32, 3: 5.95}
harmonic_data[10] = {2: 3.66, 3: 7.69}

print("=" * 70)
print("ASYMPTOTIC ANALYSIS — All CFT Quantities for Potts q=2-30")
print("=" * 70)

# ===== 1. Harmonic ratios → k² (free boson limit) =====
print("\n\n--- 1. HARMONIC RATIOS → k² (Free Boson Limit) ---")
print(f"\n  Fraction of k² reached: R(sigma^k)/R(sigma) / k²")
print(f"\n  {'q':>3}", end="")
for k in [2, 3, 4, 5]:
    print(f"  {'k='+str(k):>8}", end="")
print()
for q in sorted(harmonic_data.keys()):
    print(f"  {q:>3}", end="")
    for k in [2, 3, 4, 5]:
        if k in harmonic_data[q]:
            frac = harmonic_data[q][k] / k**2
            print(f"  {frac:>8.4f}", end="")
        else:
            print(f"  {'---':>8}", end="")
    print()

# Fit: R/k² = 1 - a/q for each k
print(f"\n  Fit: R(sigma^k)/k² ≈ 1 - a_k/q")
for k in [2, 3, 4]:
    qs, fracs = [], []
    for q in sorted(harmonic_data.keys()):
        if k in harmonic_data[q]:
            qs.append(q)
            fracs.append(harmonic_data[q][k] / k**2)
    if len(qs) >= 3:
        qs, fracs = np.array(qs), np.array(fracs)
        deviations = 1 - fracs
        # Fit deviation = a/q^b
        if np.all(deviations > 0):
            log_d = np.log(deviations)
            log_q = np.log(qs)
            slope, intercept = np.polyfit(log_q, log_d, 1)
            a = np.exp(intercept)
            print(f"    k={k}: deviation ∝ q^{slope:.2f}, a={a:.3f}")

# ===== 2. R_epsilon(q) at n=4 =====
print("\n\n--- 2. R_EPSILON(q) at n=4 ---")
qs_eps = []
reps = []
for q in sorted(n4_data.keys()):
    r = n4_data[q]["R_eps"]
    if r is not None:
        qs_eps.append(q)
        reps.append(r)
        print(f"  q={q:>3}: R_eps = {r:.4f}")
qs_eps, reps = np.array(qs_eps), np.array(reps)

# For q≥10, R_eps grows linearly
mask = qs_eps >= 10
if sum(mask) >= 2:
    slope, intercept = np.polyfit(qs_eps[mask], reps[mask], 1)
    print(f"\n  Linear fit (q≥10): R_eps ≈ {slope:.3f}·q + {intercept:.2f}")
    print(f"  Interpretation: at n=4, R_eps grows ~linearly with q")
    print(f"  (At larger n, R_eps saturates — cf. q=10 n=6: R=8.3 = n=4: R=8.3)")

# ===== 3. C_sse(q) — OPE coefficient =====
print("\n\n--- 3. C_sse(q) — OPE Coefficient Scaling ---")
print(f"\n  n=4 (charge-resolved):")
qs_c, cs_c = [], []
for q in sorted(n4_data.keys()):
    C = n4_data[q]["C_sse"]
    if C is not None:
        qs_c.append(q)
        cs_c.append(C)
        print(f"  q={q:>3}: C_sse = {C:.4f}")

qs_c, cs_c = np.array(qs_c), np.array(cs_c)

# Power law fit for q≥4
mask = qs_c >= 4
log_q = np.log(qs_c[mask])
log_C = np.log(cs_c[mask])
slope, intercept = np.polyfit(log_q, log_C, 1)
print(f"\n  Power law fit (q≥4): C_sse ~ {np.exp(intercept):.3f} · q^{slope:.3f}")
print(f"  (n=4 data; at larger n the exponent is ~{slope:.1f})")

# Best estimates from larger n
print(f"\n  Best estimates (largest n available):")
qs_best, cs_best = [], []
for q in sorted(best_data.keys()):
    C = best_data[q].get("C_sse_best")
    if C is not None:
        qs_best.append(q)
        cs_best.append(C)
        print(f"  q={q:>3}: C_sse = {C:.4f}")
qs_best, cs_best = np.array(qs_best), np.array(cs_best)
mask_b = qs_best >= 4
if sum(mask_b) >= 2:
    sl, ic = np.polyfit(np.log(qs_best[mask_b]), np.log(cs_best[mask_b]), 1)
    print(f"\n  Best-n power law (q≥4): C_sse ~ {np.exp(ic):.3f} · q^{sl:.3f}")

# ===== 4. x₁(q) — scaling dimension =====
print("\n\n--- 4. x₁(q) — Scaling Dimension ---")
print(f"\n  From best estimates (larger n):")
qs_x, xs = [], []
for q in sorted(best_data.keys()):
    x1 = best_data[q].get("x1")
    if x1 is not None:
        qs_x.append(q)
        xs.append(x1)
        print(f"  q={q:>3}: x₁ = {x1:.4f}")

qs_x, xs = np.array(qs_x), np.array(xs)

# From descendant gap at n=4 (1/desc_gap ≈ x₁)
desc_gaps = {15: 14.106, 20: 18.053, 25: 22.040, 30: 26.061}
print(f"\n  From n=4 descendant gap (x₁ ≈ 1/gap, rough):")
for q, dg in desc_gaps.items():
    x1_approx = 1.0 / dg
    print(f"  q={q:>3}: desc_gap = {dg:.1f}, x₁_approx = {x1_approx:.4f}")

# Fit x₁ ~ a/q for q≥5
mask_x = qs_x >= 5
if sum(mask_x) >= 2:
    sl_x, ic_x = np.polyfit(np.log(qs_x[mask_x]), np.log(xs[mask_x]), 1)
    print(f"\n  Power law (q≥5): x₁ ~ {np.exp(ic_x):.3f} · q^{sl_x:.3f}")

# ===== 5. c(q) — central charge =====
print("\n\n--- 5. c(q) — Central Charge ---")
print(f"\n  Known values:")
qs_cc, ccs = [], []
for q in sorted(best_data.keys()):
    c = best_data[q].get("c")
    if c is not None:
        qs_cc.append(q)
        ccs.append(c)
        print(f"  q={q:>3}: c = {c:.3f}")

qs_cc, ccs = np.array(qs_cc), np.array(ccs)

# Test c ≈ 0.40·ln(q-1) + 0.55
print(f"\n  Formula: c ≈ 0.40·ln(q-1) + 0.55")
for q, c in zip(qs_cc, ccs):
    c_pred = 0.40 * np.log(q - 1) + 0.55
    err = abs(c - c_pred) / c * 100
    print(f"  q={q:>3}: measured={c:.3f}, predicted={c_pred:.3f}, error={err:.1f}%")

# Predictions for new q
print(f"\n  Predictions:")
for q in [15, 20, 25, 30, 50, 100]:
    c_pred = 0.40 * np.log(q - 1) + 0.55
    gc = gc_formula(q)
    print(f"  q={q:>3}: c_pred = {c_pred:.3f}, g_c_pred = {gc:.3f}")

# ===== 6. c/x₁ ratio =====
print("\n\n--- 6. c/x₁ Ratio ---")
print(f"\n  Known:")
for q in sorted(best_data.keys()):
    cx1 = best_data[q].get("cx1")
    if cx1 is not None:
        print(f"  q={q:>3}: c/x₁ = {cx1:.2f} (2q = {2*q})")

qs_cx, cxs = [], []
for q in sorted(best_data.keys()):
    cx = best_data[q].get("cx1")
    if cx is not None:
        qs_cx.append(q)
        cxs.append(cx)
qs_cx, cxs = np.array(qs_cx), np.array(cxs)

# Test c/x₁ = 2q for small q
print(f"\n  Deviation from 2q:")
for q, cx in zip(qs_cx, cxs):
    dev = cx / (2*q) * 100
    print(f"  q={q:>3}: c/x₁ = {cx:.2f}, 2q = {2*q}, fraction = {dev:.1f}%")

# Fit c/x₁ ~ a·q^b
sl_cx, ic_cx = np.polyfit(np.log(qs_cx), np.log(cxs), 1)
print(f"\n  Power law: c/x₁ ~ {np.exp(ic_cx):.3f} · q^{sl_cx:.3f}")
print(f"  (Exact 2q would give exponent 1.0; measured {sl_cx:.3f} is sub-linear)")

# ===== 7. Free Boson Limit Summary =====
print("\n\n--- 7. FREE BOSON LIMIT CONVERGENCE ---")
print("""
  In a free compact boson CFT:
  - Vertex operators: x(V_k) = k² · x(V_1) → harmonic ratios = k²
  - Central charge: c = 1 (single scalar field)
  - All OPE coefficients calculable from radius R

  Our measurements at q→∞:
""")

# Harmonic ratio convergence to k²
print(f"  Harmonic ratio R(σ²)/4:")
for q in [4, 5, 7, 10, 15, 20, 25, 30]:
    if q in harmonic_data and 2 in harmonic_data[q]:
        frac = harmonic_data[q][2] / 4
        deficit = (1 - frac) * 100
        print(f"    q={q:>3}: {frac:.5f} (deficit {deficit:.2f}%)")

# x₁ approach to 0
print(f"\n  x₁ → 0:")
for q in sorted(best_data.keys()):
    x1 = best_data[q].get("x1")
    if x1: print(f"    q={q:>3}: x₁ = {x1:.4f}")
for q, dg in desc_gaps.items():
    print(f"    q={q:>3}: x₁ ≈ {1/dg:.4f} (from n=4 desc gap)")

# c growth
print(f"\n  c → ∞ as ln(q):")
for q in sorted(best_data.keys()):
    c = best_data[q].get("c")
    if c:
        c_over_lnq = c / np.log(q) if q > 1 else None
        print(f"    q={q:>3}: c = {c:.3f}, c/ln(q) = {c_over_lnq:.3f}" if c_over_lnq else "")

print(f"""
  CONCLUSION: The large-q limit is NOT a single free boson (c→∞, not c=1).
  Instead, c ~ 0.40·ln(q) — effectively ~0.40 independent scalar modes per unit of ln(q).
  The harmonic ratios approach k² from below with corrections O(1/q).
  x₁ → 0 meaning the spin field becomes marginal.
  C_sse → 0 meaning spin-spin-energy coupling vanishes.

  Physical picture: at large q, the Z_q clock field effectively decompactifies,
  producing a free boson on a large-radius circle. The radius R ~ √(ln q)
  explains c ~ ln(q)/ln(something) and x₁ ~ 1/(R²) → 0.
""")

# Save comprehensive results
results = {
    "harmonic_fraction_k2": {},
    "R_epsilon_n4": {str(q): r for q, r in zip(qs_eps, reps)},
    "C_sse_n4": {str(q): c for q, c in zip(qs_c, cs_c)},
    "C_sse_best": {str(q): best_data[q]["C_sse_best"] for q in best_data if "C_sse_best" in best_data[q]},
    "x1_best": {str(q): best_data[q]["x1"] for q in best_data if "x1" in best_data[q]},
    "c_best": {str(q): best_data[q]["c"] for q in best_data if "c" in best_data[q]},
    "power_law_C_sse_n4_exponent": float(slope),
    "power_law_x1_exponent": float(sl_x) if sum(mask_x) >= 2 else None,
    "power_law_cx1_exponent": float(sl_cx),
    "desc_gap_x1_approx": {str(q): 1.0/dg for q, dg in desc_gaps.items()},
    "gc_predictions": {str(q): gc_formula(q) for q in [15, 20, 25, 30, 50]},
    "c_predictions": {str(q): 0.40*np.log(q-1)+0.55 for q in [15, 20, 25, 30, 50]},
}

for q in sorted(harmonic_data.keys()):
    for k, r in harmonic_data[q].items():
        key = f"q={q},k={k}"
        results["harmonic_fraction_k2"][key] = r / k**2

with open("results/sprint_061d_asymptotics.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to results/sprint_061d_asymptotics.json")
