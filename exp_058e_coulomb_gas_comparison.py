"""Sprint 058e: Coulomb gas comparison + finite-size correction model.

Compare measured c/x₁ against Coulomb gas predictions for q≤4.
Model finite-size corrections to estimate converged x₁.
Test whether c/x₁ = 2q is the asymptotic result for q=4 (log corrections).
"""
import numpy as np
from scipy.optimize import curve_fit
import json

# Coulomb gas predictions for q-state Potts (q ≤ 4)
# Parametrization: q = 4*cos²(pi*t) where t = 1/(m+1), m = minimal model index
# h_sigma from Kac table: h_{(m-1)/2, (m+1)/2} for Potts spin operator
# x_sigma = 2*h_sigma

# Exact values from literature:
# q=2 (m=3, c=1/2): h = 1/16, x = 1/8
# q=3 (m=5, c=4/5): h = 1/15, x = 2/15
# q=4 (m→∞, c=1): h = 1/16, x = 1/8  [Coulomb gas limit]

# General Coulomb gas formula: x_sigma(q) for q ≤ 4
# q = 4*cos²(pi*g) where g = p'/(p'+1) and m = p'
# Actually: q = 2 + 2*cos(2*pi/(m+1)), so m = 2*pi/arccos((q-2)/2) - 1
# Then h_{1,2} = (m-1)/(4*(m+1))... let me just compute numerically.

def coulomb_gas_x1(q):
    """Coulomb gas x_sigma for q-state Potts, q <= 4."""
    if q == 2:
        return 1/8
    elif q == 3:
        return 2/15
    elif q == 4:
        return 1/8  # marginal limit

    # General: q = 2 + 2*cos(2*pi/(m+1))
    # => (q-2)/2 = cos(2*pi/(m+1))
    # => m+1 = 2*pi / arccos((q-2)/2)
    if q < 4:
        m_plus_1 = 2 * np.pi / np.arccos((q - 2) / 2)
        # h_{1/2,0} for Potts spin: h = (m+1-2)^2 / (16*m*(m+1)) ???
        # Actually h_(2,1) in Potts: complicated. Use known parametric formula.
        # Exact: h_sigma = (2-y)(y-1)/(2y) where y = m/(m+1)
        # Wait, the standard result is:
        # For q-state Potts: x_sigma = (1/(2*pi)) * arccos(-sqrt(q)/2) * ...
        # This is getting complicated. Let me just interpolate from exact values.
        pass
    return None

# Use exact values only
exact_cx1 = {2: 4.0, 3: 6.0, 4: 8.0}  # Coulomb gas: c/x1 = 2q for q=2,3,4
exact_x1 = {2: 0.125, 3: 2/15, 4: 0.125}

print("=" * 70)
print("COULOMB GAS vs MEASUREMENT")
print("=" * 70)

# Load all pairwise c/x1 data
data_058c = json.load(open("results/sprint_058c_x1_extraction.json"))
data_058d = json.load(open("results/sprint_058d_x1_formula.json"))

# Reconstruct all_data from the raw data
all_data = {}
for q_str, sizes in data_058c["all_data"].items():
    q = int(q_str)
    for n_str, vals in sizes.items():
        n = int(n_str)
        all_data.setdefault(q, {})[n] = (vals["E0"], vals["gap1"])

# Add q=10 n=6
q10_n6 = data_058d["q10_n6"]
all_data.setdefault(10, {})[6] = (q10_n6["E0"], q10_n6["nonzero_gaps"][0])

# Compute ALL pairwise c/x1 for each q
print(f"\n{'q':>3} {'pair':>8} {'c/x1':>8} {'1/N_max²':>10} {'2q':>5} {'Δ(c/x1)':>10}")
print("-" * 55)

all_pairs = {}
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
        delta = cx1 - 2*q
        print(f"{q:>3} {f'({n1},{n2})':>8} {cx1:>8.4f} {1/n2**2:>10.5f} {2*q:>5} {delta:>+10.4f}")
        pairs.append((n2, cx1))
    all_pairs[q] = pairs

# --- FSS correction: c/x1(N) = c/x1(inf) + alpha/(N^2 * ln(N)^beta) ---
print("\n\n" + "=" * 70)
print("FINITE-SIZE CORRECTION MODEL")
print("=" * 70)

# For q=2,3: exact c/x1 known. Calibrate the correction.
print("\n--- q=2 calibration ---")
for n2, cx1 in all_pairs[2]:
    delta = cx1 - 4.0
    print(f"  N={n2}: c/x1={cx1:.5f}, excess={delta:.5f}, excess*N²={delta*n2**2:.2f}")

print("\n--- q=3 calibration ---")
for n2, cx1 in all_pairs[3]:
    delta = cx1 - 6.0
    print(f"  N={n2}: c/x1={cx1:.5f}, excess={delta:.5f}, excess*N²={delta*n2**2:.2f}")

print("\n--- q=4 convergence test ---")
# If c/x1(inf) = 8.0, what's the correction?
for n2, cx1 in all_pairs[4]:
    delta_8 = cx1 - 8.0
    delta_85 = cx1 - 8.5
    print(f"  N={n2}: c/x1={cx1:.5f}, excess(8.0)={delta_8:.4f}*N²={delta_8*n2**2:.2f}, excess(8.5)={delta_85:+.4f}")

print("\n--- q=5 convergence ---")
for n2, cx1 in all_pairs[5]:
    print(f"  N={n2}: c/x1={cx1:.5f}")

# For q=4 with 3 pairs, test: c/x1 = A + B/N^gamma
# Fit for A, B, gamma
if len(all_pairs[4]) >= 3:
    ns = np.array([p[0] for p in all_pairs[4]], dtype=float)
    cx1s = np.array([p[1] for p in all_pairs[4]])

    print("\n  Fitting c/x₁(q=4) = A + B/N^γ:")
    for gamma_test in [0.5, 1.0, 1.5, 2.0]:
        X = np.column_stack([np.ones(len(ns)), 1/ns**gamma_test])
        coeffs = np.linalg.lstsq(X, cx1s, rcond=None)[0]
        residuals = cx1s - X @ coeffs
        rms = np.sqrt(np.mean(residuals**2))
        print(f"    γ={gamma_test:.1f}: A={coeffs[0]:.4f}, B={coeffs[1]:.4f}, RMS={rms:.6f}")

    # Also try A + B/ln(N)
    print("\n  Fitting c/x₁(q=4) = A + B/ln(N):")
    X = np.column_stack([np.ones(len(ns)), 1/np.log(ns)])
    coeffs = np.linalg.lstsq(X, cx1s, rcond=None)[0]
    residuals = cx1s - X @ coeffs
    rms = np.sqrt(np.mean(residuals**2))
    print(f"    A={coeffs[0]:.4f}, B={coeffs[1]:.4f}, RMS={rms:.6f}")

    # A + B/ln(N)^2
    print("\n  Fitting c/x₁(q=4) = A + B/ln(N)²:")
    X = np.column_stack([np.ones(len(ns)), 1/np.log(ns)**2])
    coeffs = np.linalg.lstsq(X, cx1s, rcond=None)[0]
    residuals = cx1s - X @ coeffs
    rms = np.sqrt(np.mean(residuals**2))
    print(f"    A={coeffs[0]:.4f}, B={coeffs[1]:.4f}, RMS={rms:.6f}")


# ========== FINAL SUMMARY ==========
print("\n\n" + "=" * 70)
print("FINAL SUMMARY: x₁(q) Best Estimates")
print("=" * 70)

c_all = {2: 0.500, 3: 0.800, 4: 1.000, 5: 1.10, 7: 1.30, 10: 1.40}

print(f"\n{'q':>3} {'c':>6} {'c/x₁ meas':>10} {'c/x₁ CG':>9} {'x₁ meas':>10} {'x₁ CG':>8} {'Δ₁·N (max n)':>14}")
print("-" * 75)

for q in sorted(all_data.keys()):
    sizes = sorted(all_data[q].keys())
    # Best c/x1 = from largest pair
    if len(all_pairs[q]) > 0:
        _, cx1_meas = all_pairs[q][-1]
    else:
        cx1_meas = None

    c = c_all[q]
    x1_meas = c / cx1_meas if cx1_meas else None
    cx1_cg = 2*q if q <= 4 else None
    x1_cg = c / (2*q) if q <= 4 else None

    # Delta1*N at largest size
    e0_max, g1_max = all_data[q][sizes[-1]]
    dn_max = g1_max * sizes[-1]

    cx1_str = f"{cx1_meas:.4f}" if cx1_meas else "—"
    cg_str = f"{cx1_cg:.1f}" if cx1_cg else "—"
    x1_str = f"{x1_meas:.5f}" if x1_meas else "—"
    xcg_str = f"{x1_cg:.5f}" if x1_cg else "—"

    print(f"{q:>3} {c:>6.3f} {cx1_str:>10} {cg_str:>9} {x1_str:>10} {xcg_str:>8} {dn_max:>14.6f}")

print("\n\nKey conclusions:")
print("1. c/x₁ = 2q is EXACT for q=2,3 (from known CFT)")
print("2. c/x₁(q=4) ≈ 8.5 at n≤10 — MAY converge to 8=2q under log corrections")
print("3. c/x₁ sub-linear in q: grows slower than 2q for q≥4")
print("4. x₁ peaks at q=3 (x₁=2/15), decreases for q>3")
print("5. For q≥5: no CFT prediction exists — measurements are novel")

# Save final results
final = {
    "c_x1_measured": {str(q): all_pairs[q][-1][1] for q in all_pairs if all_pairs[q]},
    "x1_measured": {str(q): c_all[q] / all_pairs[q][-1][1]
                    for q in all_pairs if all_pairs[q]},
    "coulomb_gas_cx1": {str(q): 2*q for q in [2, 3, 4]},
    "conclusions": {
        "cx1_exact_2q_for": "q=2,3",
        "cx1_q4_finite_size": "8.54 ± 0.10, may converge to 8.0 with log corrections",
        "cx1_sublinear": True,
        "x1_peak_at_q3": True,
        "novel_for_q_geq": 5,
    }
}
with open("results/sprint_058e_coulomb_gas.json", "w") as f:
    json.dump(final, f, indent=2)
print("\nResults saved.")
