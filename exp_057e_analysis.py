"""Sprint 057e: Analysis of scaling dimension ratios across q.

Compile all spectral data and look for universal patterns.
"""
import numpy as np
import json

# Best available data (largest n for each q)
data = {
    2: {"n": 12, "R_harmonics": [], "R_energy": 7.966,
        "degeneracy": "1 (σ self-conj)",
        "notes": "σ is Z₂ self-conjugate"},
    3: {"n": 10, "R_harmonics": [], "R_energy": 6.192,
        "degeneracy": "2 (σ, σ†)",
        "notes": ""},
    4: {"n": 8, "R_harmonics": [1.677], "R_energy": 6.578,
        "degeneracy": "2+1 (σ,σ†,σ²)",
        "notes": "σ² self-conjugate"},
    5: {"n": 6, "R_harmonics": [2.409], "R_energy": 7.143,
        "degeneracy": "2+2 (σ,σ†,σ²,σ²†)",
        "notes": ""},
    7: {"n": 6, "R_harmonics": [3.323, 5.948], "R_energy": 8.072,
        "degeneracy": "2+2+2",
        "notes": "3 conjugate pairs"},
}

print("=" * 70)
print("SCALING DIMENSION RATIO TABLE")
print("=" * 70)
print(f"{'q':>3} {'n':>3} {'R(σ²/σ)':>10} {'R(σ³/σ)':>10} {'R(ε/σ)':>10} {'Spin degen.':>15}")
print("-" * 70)
for q in [2, 3, 4, 5, 7]:
    d = data[q]
    r2 = f"{d['R_harmonics'][0]:.3f}" if d['R_harmonics'] else "—"
    r3 = f"{d['R_harmonics'][1]:.3f}" if len(d['R_harmonics']) > 1 else "—"
    re = f"{d['R_energy']:.3f}"
    print(f"{q:>3} {d['n']:>3} {r2:>10} {r3:>10} {re:>10} {d['degeneracy']:>15}")

# Test formula: x(σ^k)/x(σ) = ?
# Hypothesis 1: k² (free boson)
# Hypothesis 2: sin²(πk/q)/sin²(π/q) (quantum group)
# Hypothesis 3: (1-cos(2πk/q))/(1-cos(2π/q)) (dispersion)
# Hypothesis 4: empirical fit

print("\n\n" + "=" * 70)
print("HARMONIC RATIO FORMULA TEST")
print("=" * 70)

for q in [4, 5, 7]:
    d = data[q]
    print(f"\nq={q} (n={d['n']}):")
    floor_pairs = (q-1)//2
    for k_idx, k in enumerate(range(2, floor_pairs + 1 + (1 if q % 2 == 0 else 0))):
        if k_idx >= len(d['R_harmonics']):
            break
        measured = d['R_harmonics'][k_idx]

        # Predictions
        free_boson = k**2
        sin2 = np.sin(np.pi*k/q)**2 / np.sin(np.pi/q)**2
        disp = (1 - np.cos(2*np.pi*k/q)) / (1 - np.cos(2*np.pi/q))
        # Chebyshev U_{k-1}(cos(π/q)) related
        cheby = (np.sin(np.pi*k/q) / np.sin(np.pi/q))**2

        print(f"  σ^{k}/σ: measured = {measured:.3f}")
        print(f"    k²           = {free_boson:.3f}  (err {abs(free_boson-measured)/measured*100:.1f}%)")
        print(f"    sin²(πk/q)/sin²(π/q) = {sin2:.3f}  (err {abs(sin2-measured)/measured*100:.1f}%)")
        print(f"    dispersion   = {disp:.3f}  (err {abs(disp-measured)/measured*100:.1f}%)")

# Energy field ratio trend
print("\n\n" + "=" * 70)
print("ENERGY FIELD RATIO R(ε/σ)")
print("=" * 70)
for q in [2, 3, 4, 5, 7]:
    d = data[q]
    # Known exact values for q=2,3
    exact = {2: 8.0, 3: 6.0}
    ex_str = f" (exact: {exact[q]})" if q in exact else ""
    print(f"  q={q}: R = {d['R_energy']:.3f}{ex_str}")

# Number of spin fields before energy
n_spin = {2: 1, 3: 2, 4: 3, 5: 4, 7: 6}
print("\n  Number of spin field levels before energy:")
for q in [2, 3, 4, 5, 7]:
    print(f"    q={q}: {n_spin[q]} spin fields (= q-1)")

# Central charge connection
# c(q) ~ ln(q), and the spin sector grows as q-1
# Effective degrees of freedom grow with q
print("\n\n" + "=" * 70)
print("KEY FINDING: Operator Content Structure")
print("=" * 70)
print("""
The Z_q Potts CFT spectrum has:
  - (q-1) spin fields organized as ⌊(q-1)/2⌋ conjugate pairs
    (plus 1 self-conjugate if q even)
  - Harmonic ratios x(σ^k)/x(σ) grow subquadratically with k
  - Energy field appears AFTER all spin harmonics
  - R(ε/σ) has minimum at q=3 (=6), increases for q>3

The growing number of primary fields below the energy scale
is consistent with c(q) ~ ln(q) — more primaries means
more effective massless modes, hence higher central charge.
""")

# Save summary
summary = {
    "data": {str(q): {
        "n": data[q]["n"],
        "R_harmonics": data[q]["R_harmonics"],
        "R_energy": data[q]["R_energy"],
        "degeneracy": data[q]["degeneracy"],
    } for q in data},
    "n_spin_fields": n_spin,
    "key_finding": "Z_q Potts has (q-1) spin fields as primaries below energy field",
}
with open("results/sprint_057e_analysis.json", "w") as f:
    json.dump(summary, f, indent=2)
print("Analysis saved.")
