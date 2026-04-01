"""Sprint 058c: Comprehensive x1(q) extraction from multi-size periodic chain data.

Method: CFT on periodic chain gives
  E0(N)/N = e_inf - pi*v*c/(6*N^2) + O(1/N^3)
  Delta_1(N) = 2*pi*v*x1/N * (1 + a/N^2 + ...)

From two sizes n1, n2:
  v*c = [E0(n1)/n1 - E0(n2)/n2] * 6 / [pi * (1/n2^2 - 1/n1^2)]
  v*x1 = Delta_1(N) * N / (2*pi)
  => c/x1 = (v*c) / (v*x1)

With known c, get x1 = c / (c/x1).
"""
import numpy as np
import json

# Collect all data: {q: {n: (E0, gap1)}}
all_data = {}

# q=2 (Ising, g_c=1.0, TFIChain convention)
ising = json.load(open("results/sprint_057a_ising_spectrum.json"))
for key in ising:
    n = ising[key]["n"]
    all_data.setdefault(2, {})[n] = (ising[key]["E0"], ising[key]["gaps"][1])

# q=3
q3 = json.load(open("results/sprint_057b_q3_spectrum.json"))
for key in q3:
    n = q3[key]["n"]
    all_data.setdefault(3, {})[n] = (q3[key]["E0"], q3[key]["nonzero_gaps"][0])

# q=4,5 from Sprint 057c
q45 = json.load(open("results/sprint_057c_q4q5_spectrum.json"))
for key in q45:
    q = q45[key]["q"]
    n = q45[key]["n"]
    all_data.setdefault(q, {})[n] = (q45[key]["E0"], q45[key]["nonzero_gaps"][0])

# q=4 n=8 and q=7 from Sprint 057d
q47 = json.load(open("results/sprint_057d_q7_spectrum.json"))
for key in q47:
    q_val = q47[key]["q"]
    n = q47[key]["n"]
    all_data.setdefault(q_val, {})[n] = (q47[key]["eigenvalues"][0], q47[key]["nonzero_gaps"][0])

# q=10 from Sprint 057f
q10 = json.load(open("results/sprint_057f_q10_absolute.json"))
data_q10 = q10["q=10_n=4"]
all_data.setdefault(10, {})[4] = (data_q10["eigenvalues"][0], data_q10["nonzero_gaps"][0])

# New data from Sprint 058
q4_new = json.load(open("results/sprint_058a_q4_n10.json"))
all_data.setdefault(4, {})[10] = (q4_new["E0"], q4_new["nonzero_gaps"][0])

q5_new = json.load(open("results/sprint_058b_q5_n8.json"))
all_data.setdefault(5, {})[8] = (q5_new["E0"], q5_new["nonzero_gaps"][0])

# Known central charges
c_exact = {2: 0.500, 3: 0.800, 4: 1.000}
c_est = {5: 1.10, 7: 1.30, 10: 1.40}
c_all = {**c_exact, **c_est}

print("=" * 80)
print("COMPREHENSIVE x₁(q) EXTRACTION")
print("=" * 80)

print("\n--- Available data ---")
for q in sorted(all_data.keys()):
    sizes = sorted(all_data[q].keys())
    print(f"  q={q}: sizes={sizes}")
    for n in sizes:
        e0, g1 = all_data[q][n]
        print(f"    n={n}: E0/n={e0/n:.8f}, Delta1={g1:.8f}, Delta1*n={g1*n:.6f}")

# --- Method 1: Pairwise c/x1 extraction ---
print("\n--- Method 1: Pairwise c/x₁ from consecutive size pairs ---")
print(f"{'q':>3} {'n1,n2':>10} {'v*c':>10} {'v*x1_1':>10} {'v*x1_2':>10} {'c/x1':>8} {'x1 (c known)':>14} {'x1 exact':>10}")
print("-" * 90)

x1_results = {}

for q in sorted(all_data.keys()):
    sizes = sorted(all_data[q].keys())
    pairs = []
    for i in range(len(sizes) - 1):
        n1, n2 = sizes[i], sizes[i+1]
        e0_1, g1_1 = all_data[q][n1]
        e0_2, g1_2 = all_data[q][n2]

        # v*c from Casimir energy
        vc = (e0_1/n1 - e0_2/n2) * 6 / (np.pi * (1/n2**2 - 1/n1**2))

        # v*x1 from each size
        vx1_1 = g1_1 * n1 / (2 * np.pi)
        vx1_2 = g1_2 * n2 / (2 * np.pi)
        vx1_avg = (vx1_1 + vx1_2) / 2

        # c/x1
        cx1 = vc / vx1_avg

        c = c_all.get(q)
        x1 = c / cx1 if c else None

        exact_x1 = {2: 0.125, 3: 2/15}
        ex_str = f"{exact_x1[q]:.6f}" if q in exact_x1 else "—"
        x1_str = f"{x1:.6f}" if x1 else "—"

        print(f"{q:>3} {f'{n1},{n2}':>10} {vc:>10.6f} {vx1_1:>10.6f} {vx1_2:>10.6f} {cx1:>8.4f} {x1_str:>14} {ex_str:>10}")
        pairs.append({"n1": n1, "n2": n2, "vc": vc, "vx1_avg": vx1_avg, "cx1": cx1, "x1": x1})

    x1_results[q] = pairs

# --- Method 2: Multi-size extrapolation ---
print("\n\n--- Method 2: Multi-size fit ---")
print("Fit Delta1*N = 2*pi*v*x1 * (1 + a/N^2) and E0/N = e_inf - pi*v*c/(6*N^2)")

from numpy.polynomial import polynomial as P

for q in sorted(all_data.keys()):
    sizes = sorted(all_data[q].keys())
    if len(sizes) < 3:
        continue

    ns = np.array(sizes, dtype=float)
    e0s = np.array([all_data[q][n][0] for n in sizes])
    g1s = np.array([all_data[q][n][1] for n in sizes])

    # Delta1 * N vs 1/N^2 — linear fit to extrapolate to N->inf
    delta_n = g1s * ns
    inv_n2 = 1.0 / ns**2

    # Fit Delta1*N = A + B/N^2
    coeffs = np.polyfit(inv_n2, delta_n, 1)
    A_delta = coeffs[1]  # intercept = 2*pi*v*x1 at N->inf
    print(f"\nq={q}: sizes={sizes}")
    print(f"  Delta1*N data: {[f'{d:.6f}' for d in delta_n]}")
    print(f"  Linear fit: Delta1*N = {A_delta:.6f} + {coeffs[0]:.4f}/N^2")
    print(f"  => v*x1 (extrapolated) = {A_delta/(2*np.pi):.6f}")

    # E0/N vs 1/N^2 — linear fit to get e_inf and v*c
    e0_n = e0s / ns
    coeffs_e = np.polyfit(inv_n2, e0_n, 1)
    e_inf = coeffs_e[1]
    vc_fit = -coeffs_e[0] * 6 / np.pi
    vx1_fit = A_delta / (2 * np.pi)
    cx1_fit = vc_fit / vx1_fit

    c = c_all.get(q)
    x1_fit = c / cx1_fit if c else None

    print(f"  E0/N fit: e_inf={e_inf:.8f}, v*c={vc_fit:.6f}")
    print(f"  c/x1 = {cx1_fit:.4f}")
    if x1_fit:
        exact_x1 = {2: 0.125, 3: 2/15}
        err = ""
        if q in exact_x1:
            err = f" (error: {abs(x1_fit - exact_x1[q])/exact_x1[q]*100:.1f}%)"
        print(f"  x1 = {x1_fit:.6f} (using c={c}){err}")

# --- Summary table ---
print("\n\n" + "=" * 80)
print("SUMMARY: x₁(q) Best Estimates")
print("=" * 80)
print(f"{'q':>3} {'c':>6} {'x1 (best)':>12} {'x1 exact':>12} {'method':>25} {'c/x1':>8}")
print("-" * 75)

exact_x1_all = {2: 0.125, 3: 2/15}
summary = {}

for q in sorted(all_data.keys()):
    sizes = sorted(all_data[q].keys())
    pairs = x1_results[q]

    # Use largest size pair for best estimate
    best = pairs[-1] if pairs else None
    c = c_all.get(q)

    if best and best["x1"]:
        x1 = best["x1"]
        cx1 = best["cx1"]
        method = f"pair ({best['n1']},{best['n2']})"
    else:
        x1 = None
        cx1 = best["cx1"] if best else None
        method = "—"

    ex_str = f"{exact_x1_all[q]:.6f}" if q in exact_x1_all else "—"
    x1_str = f"{x1:.6f}" if x1 else "—"
    cx1_str = f"{cx1:.4f}" if cx1 else "—"

    print(f"{q:>3} {c:>6.2f} {x1_str:>12} {ex_str:>12} {method:>25} {cx1_str:>8}")
    summary[q] = {"c": c, "x1": x1, "cx1": cx1, "method": method, "sizes": sizes}

# Save results
results = {
    "all_data": {str(q): {str(n): {"E0": all_data[q][n][0], "gap1": all_data[q][n][1]}
                          for n in all_data[q]}
                 for q in all_data},
    "pairwise": {str(q): [{"n1": p["n1"], "n2": p["n2"],
                            "vc": p["vc"], "vx1_avg": p["vx1_avg"],
                            "cx1": p["cx1"], "x1": p["x1"]}
                           for p in x1_results[q]]
                 for q in x1_results},
    "summary": {str(q): {"c": summary[q]["c"], "x1": summary[q]["x1"],
                          "cx1": summary[q]["cx1"]}
                for q in summary}
}

with open("results/sprint_058c_x1_extraction.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nResults saved to results/sprint_058c_x1_extraction.json")
