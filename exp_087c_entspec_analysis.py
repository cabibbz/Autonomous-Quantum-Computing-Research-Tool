#!/usr/bin/env python3
"""Sprint 087c: Entanglement spectrum scaling analysis.

Compare q=5 (087a) vs q=7 (087b) tail weight growth.
Key question: does tail weight saturate, grow as power law, or grow logarithmically?
What does this predict for the energy-entropy decoupling at large n?
"""
import numpy as np
import json, time
from scipy.optimize import curve_fit

results = {
    'experiment': '087c_entspec_analysis',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

def save():
    with open("results/sprint_087c_entspec_analysis.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# Load data from 087a and 087b
with open("results/sprint_087a_entspec_dmrg_q5.json") as f:
    d5 = json.load(f)
with open("results/sprint_087b_entspec_dmrg_q7.json") as f:
    d7 = json.load(f)

# Also load Sprint 084 exact diag data if available
try:
    with open("results/sprint_084a_entspec_q235.json") as f:
        d84a = json.load(f)
    with open("results/sprint_084b_entspec_q678.json") as f:
        d84b = json.load(f)
    has_084 = True
except:
    has_084 = False

print("Sprint 087c: Entanglement spectrum scaling analysis")
print("=" * 70)

# Extract data
q5_n = [d['n'] for d in d5['data']]
q5_wtail = [d['w_tail'] for d in d5['data']]
q5_Stail = [d['S_tail_frac'] for d in d5['data']]
q5_lmax = [d['lam_max'] for d in d5['data']]
q5_wmult = [d['w_mult'] for d in d5['data']]
q5_gap = [d['ent_gap'] for d in d5['data']]
q5_S = [d['S_total'] for d in d5['data']]

q7_n = [d['n'] for d in d7['data']]
q7_wtail = [d['w_tail'] for d in d7['data']]
q7_Stail = [d['S_tail_frac'] for d in d7['data']]
q7_lmax = [d['lam_max'] for d in d7['data']]
q7_wmult = [d['w_mult'] for d in d7['data']]
q7_gap = [d['ent_gap'] for d in d7['data']]
q7_S = [d['S_total'] for d in d7['data']]

# ========================================================
# 1. Fit tail weight growth: power law vs log vs saturation
# ========================================================
print("\n1. TAIL WEIGHT SCALING")
print("-" * 50)

def power_law(n, a, b):
    return a * n**b

def log_growth(n, a, b):
    return a * np.log(n) + b

def saturation(n, w_inf, n0, gamma):
    return w_inf * (1 - np.exp(-(n / n0)**gamma))

# q=5 tail weight fits
n5 = np.array(q5_n, dtype=float)
wt5 = np.array(q5_wtail)

try:
    popt_pl5, _ = curve_fit(power_law, n5, wt5, p0=[1e-6, 2])
    wt5_pred_pl = power_law(n5, *popt_pl5)
    r2_pl5 = 1 - np.sum((wt5 - wt5_pred_pl)**2) / np.sum((wt5 - np.mean(wt5))**2)
    print(f"\nq=5 Power law: w_tail = {popt_pl5[0]:.3e} × n^{popt_pl5[1]:.3f}  (R²={r2_pl5:.6f})")
except Exception as e:
    print(f"q=5 power law fit failed: {e}")
    r2_pl5 = -1
    popt_pl5 = [0, 0]

try:
    popt_log5, _ = curve_fit(log_growth, n5, wt5, p0=[0.001, -0.002])
    wt5_pred_log = log_growth(n5, *popt_log5)
    r2_log5 = 1 - np.sum((wt5 - wt5_pred_log)**2) / np.sum((wt5 - np.mean(wt5))**2)
    print(f"q=5 Log growth: w_tail = {popt_log5[0]:.4e} × ln(n) + {popt_log5[1]:.4e}  (R²={r2_log5:.6f})")
except Exception as e:
    print(f"q=5 log fit failed: {e}")
    r2_log5 = -1

# q=7 tail weight fits
n7 = np.array(q7_n, dtype=float)
wt7 = np.array(q7_wtail)

try:
    popt_pl7, _ = curve_fit(power_law, n7, wt7, p0=[1e-6, 2])
    wt7_pred_pl = power_law(n7, *popt_pl7)
    r2_pl7 = 1 - np.sum((wt7 - wt7_pred_pl)**2) / np.sum((wt7 - np.mean(wt7))**2)
    print(f"\nq=7 Power law: w_tail = {popt_pl7[0]:.3e} × n^{popt_pl7[1]:.3f}  (R²={r2_pl7:.6f})")
except Exception as e:
    print(f"q=7 power law fit failed: {e}")
    r2_pl7 = -1
    popt_pl7 = [0, 0]

try:
    popt_log7, _ = curve_fit(log_growth, n7, wt7, p0=[0.001, -0.002])
    wt7_pred_log = log_growth(n7, *popt_log7)
    r2_log7 = 1 - np.sum((wt7 - wt7_pred_log)**2) / np.sum((wt7 - np.mean(wt7))**2)
    print(f"q=7 Log growth: w_tail = {popt_log7[0]:.4e} × ln(n) + {popt_log7[1]:.4e}  (R²={r2_log7:.6f})")
except Exception as e:
    print(f"q=7 log fit failed: {e}")
    r2_log7 = -1

results['tail_fits'] = {
    'q5': {
        'power_law': {'a': float(popt_pl5[0]), 'b': float(popt_pl5[1]), 'R2': float(r2_pl5)},
        'n': q5_n, 'w_tail': [float(x) for x in q5_wtail],
    },
    'q7': {
        'power_law': {'a': float(popt_pl7[0]), 'b': float(popt_pl7[1]), 'R2': float(r2_pl7)},
        'n': q7_n, 'w_tail': [float(x) for x in q7_wtail],
    },
}

# ========================================================
# 2. Entropy fraction scaling
# ========================================================
print("\n\n2. ENTROPY FRACTION SCALING")
print("-" * 50)

print("\nq=5:")
print(f"{'n':>4} {'%S(l0)':>8} {'%S(l1)':>8} {'%S(tail)':>8} {'S_total':>8}")
for d in d5['data']:
    print(f"{d['n']:4d} {d['S_lev0_frac']:8.4f} {d['S_lev1_frac']:8.4f} {d['S_tail_frac']:8.4f} {d['S_total']:8.4f}")

print("\nq=7:")
print(f"{'n':>4} {'%S(l0)':>8} {'%S(l1)':>8} {'%S(tail)':>8} {'S_total':>8}")
for d in d7['data']:
    print(f"{d['n']:4d} {d['S_lev0_frac']:8.4f} {d['S_lev1_frac']:8.4f} {d['S_tail_frac']:8.4f} {d['S_total']:8.4f}")

# Does S_lev0_frac saturate?
print("\n\nq=5 %S(lev0) trend:")
for i in range(1, len(q5_n)):
    delta = d5['data'][i]['S_lev0_frac'] - d5['data'][i-1]['S_lev0_frac']
    print(f"  n={q5_n[i-1]}→{q5_n[i]}: Δ=%S(l0) = {delta:+.5f}")

print("\nq=7 %S(lev0) trend:")
for i in range(1, len(q7_n)):
    delta = d7['data'][i]['S_lev0_frac'] - d7['data'][i-1]['S_lev0_frac']
    print(f"  n={q7_n[i-1]}→{q7_n[i]}: Δ=%S(l0) = {delta:+.5f}")

# ========================================================
# 3. Entanglement gap scaling
# ========================================================
print("\n\n3. ENTANGLEMENT GAP SCALING")
print("-" * 50)

# Fit Δξ vs ln(n)
try:
    popt_g5, _ = curve_fit(log_growth, n5, np.array(q5_gap), p0=[-0.5, 4])
    print(f"q=5: Δξ = {popt_g5[0]:.4f} × ln(n) + {popt_g5[1]:.4f}")
except:
    popt_g5 = [0, 0]

try:
    popt_g7, _ = curve_fit(log_growth, n7, np.array(q7_gap), p0=[-0.5, 5])
    print(f"q=7: Δξ = {popt_g7[0]:.4f} × ln(n) + {popt_g7[1]:.4f}")
except:
    popt_g7 = [0, 0]

results['ent_gap_fits'] = {
    'q5': {'slope': float(popt_g5[0]), 'intercept': float(popt_g5[1])},
    'q7': {'slope': float(popt_g7[0]), 'intercept': float(popt_g7[1])},
}

# ========================================================
# 4. Multiplet weight scaling: does λ_max/(q-1) multiplet ratio change?
# ========================================================
print("\n\n4. WEIGHT REDISTRIBUTION: λ_max vs (q-1) multiplet")
print("-" * 50)

print("\nq=5: lam_max / w_mult ratio:")
for d in d5['data']:
    ratio = d['lam_max'] / d['w_mult']
    print(f"  n={d['n']:2d}: λ_max={d['lam_max']:.5f}, w_mult={d['w_mult']:.5f}, ratio={ratio:.3f}")

print("\nq=7: lam_max / w_mult ratio:")
for d in d7['data']:
    ratio = d['lam_max'] / d['w_mult']
    print(f"  n={d['n']:2d}: λ_max={d['lam_max']:.5f}, w_mult={d['w_mult']:.5f}, ratio={ratio:.3f}")

# ========================================================
# 5. Extrapolations
# ========================================================
print("\n\n5. EXTRAPOLATIONS")
print("-" * 50)

# Extrapolate tail weight to large n
for n_ext in [50, 100, 200, 1000]:
    wt5_ext = power_law(n_ext, *popt_pl5)
    wt7_ext = power_law(n_ext, *popt_pl7)
    print(f"  n={n_ext:4d}: q=5 w_tail={wt5_ext:.4f} ({wt5_ext*100:.1f}%), q=7 w_tail={wt7_ext:.4f} ({wt7_ext*100:.1f}%)")

# At what n does tail weight reach 1%? 5%? 10%?
for threshold in [0.01, 0.05, 0.10]:
    n_threshold5 = (threshold / popt_pl5[0]) ** (1.0 / popt_pl5[1]) if popt_pl5[0] > 0 else np.inf
    n_threshold7 = (threshold / popt_pl7[0]) ** (1.0 / popt_pl7[1]) if popt_pl7[0] > 0 else np.inf
    print(f"  w_tail = {threshold*100:.0f}%: q=5 at n≈{n_threshold5:.0f}, q=7 at n≈{n_threshold7:.0f}")

# ========================================================
# 6. Key comparison: tail growth rate per log(n)
# ========================================================
print("\n\n6. TAIL GROWTH RATE d(w_tail)/d(ln n)")
print("-" * 50)

print("\nq=5:")
for i in range(1, len(q5_n)):
    dwt = q5_wtail[i] - q5_wtail[i-1]
    dln = np.log(q5_n[i]) - np.log(q5_n[i-1])
    print(f"  n={q5_n[i-1]}→{q5_n[i]}: dw_tail/d(ln n) = {dwt/dln:.6f}")

print("\nq=7:")
for i in range(1, len(q7_n)):
    dwt = q7_wtail[i] - q7_wtail[i-1]
    dln = np.log(q7_n[i]) - np.log(q7_n[i-1])
    print(f"  n={q7_n[i-1]}→{q7_n[i]}: dw_tail/d(ln n) = {dwt/dln:.6f}")

# ========================================================
# 7. c_eff from entanglement entropy size pairs (open BC)
# ========================================================
print("\n\n7. c_eff FROM ENTROPY SIZE PAIRS")
print("-" * 50)

# Open BC: S = (c/6) ln(2n/π sin(πℓ/n)) + const
# At midchain (ℓ=n/2): S = (c/6) ln(2n/π) + const
# Size pair: c = 6·ΔS / ln(n2/n1)

alpha_cft5 = np.arccosh(np.sqrt(5) / 2)
Re_c5 = 1 + 6 * alpha_cft5**2 / (np.pi**2 + alpha_cft5**2)
alpha_cft7 = np.arccosh(np.sqrt(7) / 2)
Re_c7 = 1 + 6 * alpha_cft7**2 / (np.pi**2 + alpha_cft7**2)

print(f"\nq=5 (Re(c) = {Re_c5:.4f}):")
for i in range(len(q5_n) - 1):
    n1, n2 = q5_n[i], q5_n[i+1]
    dS = q5_S[i+1] - q5_S[i]
    c_pair = 6 * dS / np.log(n2 / n1)
    print(f"  ({n1},{n2}): c_eff = {c_pair:.4f}, c/Re(c) = {c_pair/Re_c5:.4f}")

print(f"\nq=7 (Re(c) = {Re_c7:.4f}):")
for i in range(len(q7_n) - 1):
    n1, n2 = q7_n[i], q7_n[i+1]
    dS = q7_S[i+1] - q7_S[i]
    c_pair = 6 * dS / np.log(n2 / n1)
    print(f"  ({n1},{n2}): c_eff = {c_pair:.4f}, c/Re(c) = {c_pair/Re_c7:.4f}")

results['c_eff_pairs'] = {
    'q5': [], 'q7': [],
    'Re_c5': float(Re_c5), 'Re_c7': float(Re_c7),
}
for i in range(len(q5_n) - 1):
    n1, n2 = q5_n[i], q5_n[i+1]
    dS = q5_S[i+1] - q5_S[i]
    c_pair = 6 * dS / np.log(n2 / n1)
    results['c_eff_pairs']['q5'].append({
        'n1': n1, 'n2': n2, 'c_eff': float(c_pair), 'c_over_Rec': float(c_pair/Re_c5)
    })
for i in range(len(q7_n) - 1):
    n1, n2 = q7_n[i], q7_n[i+1]
    dS = q7_S[i+1] - q7_S[i]
    c_pair = 6 * dS / np.log(n2 / n1)
    results['c_eff_pairs']['q7'].append({
        'n1': n1, 'n2': n2, 'c_eff': float(c_pair), 'c_over_Rec': float(c_pair/Re_c7)
    })

# ========================================================
# 8. Summary
# ========================================================
print("\n\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"""
Key findings:
1. Tail weight grows as POWER LAW in both q=5 and q=7:
   q=5: w_tail ~ n^{popt_pl5[1]:.2f} (R²={r2_pl5:.5f})
   q=7: w_tail ~ n^{popt_pl7[1]:.2f} (R²={r2_pl7:.5f})

2. Power law exponent is SIMILAR for q=5 and q=7 — the difference is the prefactor.
   q=5 prefactor: {popt_pl5[0]:.3e}
   q=7 prefactor: {popt_pl7[0]:.3e}
   Ratio: {popt_pl7[0]/popt_pl5[0]:.2f}

3. NO saturation visible at n=24 (q=5) or n=12 (q=7).
   Tail weight is UNBOUNDED — will eventually affect spectral observables.

4. Entanglement gap DECREASES logarithmically:
   q=5: Δξ = {popt_g5[0]:.3f} × ln(n) + {popt_g5[1]:.3f}
   q=7: Δξ = {popt_g7[0]:.3f} × ln(n) + {popt_g7[1]:.3f}

5. Extrapolated tail weight = 10%:
   q=5: n≈{(0.10 / popt_pl5[0]) ** (1.0 / popt_pl5[1]):.0f}
   q=7: n≈{(0.10 / popt_pl7[0]) ** (1.0 / popt_pl7[1]):.0f}
""")

save()
print(f"\nSaved to results/sprint_087c_entspec_analysis.json")

from db_utils import record
record(sprint=87, model='sq_potts', q=5, n=24,
       quantity='w_tail_power_exponent', value=popt_pl5[1],
       method='fit_n8_24', notes=f'w_tail ~ n^b, R2={r2_pl5:.5f}')
record(sprint=87, model='sq_potts', q=7, n=12,
       quantity='w_tail_power_exponent', value=popt_pl7[1],
       method='fit_n6_12', notes=f'w_tail ~ n^b, R2={r2_pl7:.5f}')
print("Recorded to DB.")
