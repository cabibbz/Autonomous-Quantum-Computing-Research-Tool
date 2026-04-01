#!/usr/bin/env python3
"""Sprint 056a: c(q) formula fitting and predictions.
Known exact: c(2)=0.5, c(3)=0.8, c(4)=1.0.
Measured: c(5) ≈ 1.10 ± 0.10.
Goal: fit formulas, predict c(7) and c(10) for testing.
"""
import numpy as np, json
from scipy.optimize import curve_fit

# Data: exact CFT values for q=2,3,4, measured for q=5
q_exact = np.array([2, 3, 4])
c_exact = np.array([0.500, 0.800, 1.000])
q_meas = np.array([5])
c_meas = np.array([1.10])  # central estimate from Sprint 055
c_meas_err = np.array([0.10])

q_all = np.array([2, 3, 4, 5])
c_all = np.array([0.500, 0.800, 1.000, 1.10])

# --- Coulomb gas analytic continuation ---
# Standard: sqrt(q) = -2*cos(pi*g), c = 1 - 6*(1-g)^2/g
# q=2: g=3/4, q=3: g=5/6, q=4: g=1
# For q>4: g becomes complex. Let's compute anyway.
def coulomb_gas_c(q):
    """Coulomb gas central charge. Real for q<=4, complex for q>4."""
    sq = np.sqrt(q)
    g = np.arccos(-sq / 2) / np.pi  # complex for q>4
    return 1 - 6 * (1 - g)**2 / g

print("=== Coulomb Gas (exact for q<=4) ===")
for q in [2, 3, 4]:
    c = coulomb_gas_c(q)
    print(f"  q={q}: c = {c.real:.6f} (exact: {c_exact[q-2]:.3f})")

print("\n=== Coulomb Gas analytic continuation q>4 (complex) ===")
for q in [5, 7, 10]:
    c = coulomb_gas_c(q)
    print(f"  q={q}: c = {c.real:.4f} + {c.imag:.4f}i  |c| = {abs(c):.4f}")

# --- Fit 1: Quadratic through exact points ---
# 3 exact points uniquely determine a quadratic
coeffs_quad = np.polyfit([2, 3, 4], [0.5, 0.8, 1.0], 2)
print(f"\n=== Quadratic fit (exact q=2,3,4) ===")
print(f"  c(q) = {coeffs_quad[0]:.4f}*q^2 + {coeffs_quad[1]:.4f}*q + {coeffs_quad[2]:.4f}")
for q in [2, 3, 4, 5, 7, 10]:
    c = np.polyval(coeffs_quad, q)
    print(f"  q={q}: c = {c:.4f}")

# --- Fit 2: Logarithmic c(q) = a*ln(q) + b ---
def log_func(q, a, b):
    return a * np.log(q) + b

popt_log, _ = curve_fit(log_func, q_all, c_all)
print(f"\n=== Logarithmic fit c(q) = {popt_log[0]:.4f}*ln(q) + {popt_log[1]:.4f} ===")
for q in [2, 3, 4, 5, 7, 10]:
    print(f"  q={q}: c = {log_func(q, *popt_log):.4f}")
residuals_log = c_all - log_func(q_all, *popt_log)
print(f"  Residuals: {residuals_log}")

# --- Fit 3: Power law c(q) = a*(q-1)^b + d ---
def power_func(q, a, b):
    return a * (q - 1)**b

# Two-param through q=2,3,4 + 5
popt_pow, _ = curve_fit(power_func, q_all, c_all, p0=[0.5, 0.5])
print(f"\n=== Power law c(q) = {popt_pow[0]:.4f}*(q-1)^{popt_pow[1]:.4f} ===")
for q in [2, 3, 4, 5, 7, 10]:
    print(f"  q={q}: c = {power_func(q, *popt_pow):.4f}")
residuals_pow = c_all - power_func(q_all, *popt_pow)
print(f"  Residuals: {residuals_pow}")

# --- Fit 4: c(q) = 1 - 6/[m(m+1)] with m = f(q), analytic continuation ---
# For q<=4: m+1 = pi/arccos(sqrt(q)/2)
# For q>4: try m+1 = pi/arccosh(sqrt(q)/2) (real analytic continuation)
def c_from_m(m):
    return 1 - 6 / (m * (m + 1))

def m_analytic(q):
    """Analytic continuation: m+1 = pi/arccos(sqrt(q)/2) for q<=4,
    m+1 = pi/arccosh(sqrt(q)/2) for q>4 (taking real part)."""
    sq2 = np.sqrt(q) / 2
    if sq2 <= 1:
        return np.pi / np.arccos(sq2) - 1
    else:
        return np.pi / np.arccosh(sq2) - 1

print(f"\n=== Minimal model analytic continuation (m -> real for q>4) ===")
for q in [2, 3, 4, 4.5, 5, 7, 10, 20]:
    m = m_analytic(q)
    c = c_from_m(m)
    print(f"  q={q:5.1f}: m={m:.4f}, c = {c:.4f}")

# --- Fit 5: Coulomb gas with REAL continuation ---
# c = 1 - 6*(1-g)^2/g where g = arccos(-sqrt(q)/2)/pi for q<=4
# For q>4: use g_real = 1 + alpha/pi where alpha = arccosh(sqrt(q)/2)
# Then c = 1 + 6*alpha^2 / (pi^2 + pi*alpha) ... but this has issues.
# Instead: direct real continuation via g(q) = Re[arccos(-sqrt(q)/2)/pi]
# At q=4, g=1. For q>4, g continues past 1.
# Alternative: g = 1 + (1/pi)*arccosh(sqrt(q)/2)
def coulomb_real_cont(q):
    if q <= 4:
        g = np.arccos(-np.sqrt(q) / 2) / np.pi
    else:
        g = 1 + np.arccosh(np.sqrt(q) / 2) / np.pi
    return 1 - 6 * (1 - g)**2 / g

print(f"\n=== Coulomb gas REAL continuation (g>1 for q>4) ===")
for q in [2, 3, 4, 5, 7, 10, 20]:
    c = coulomb_real_cont(q)
    print(f"  q={q:3d}: c = {c:.4f}")

# --- Fit 6: Simple formula c(q) = (q-1)/q * something? ---
# c(2)=0.5=1/2, c(3)=0.8=4/5, c(4)=1=1.
# Try c = 1 - a/(q-1)^b
def inv_power(q, a, b):
    return 1 - a / (q - 1)**b

popt_inv, _ = curve_fit(inv_power, q_all, c_all, p0=[0.5, 1.0])
print(f"\n=== Inverse power c(q) = 1 - {popt_inv[0]:.4f}/(q-1)^{popt_inv[1]:.4f} ===")
for q in [2, 3, 4, 5, 7, 10, 20]:
    print(f"  q={q:3d}: c = {inv_power(q, *popt_inv):.4f}")

# --- Summary: predictions for q=7 and q=10 ---
print("\n" + "="*60)
print("PREDICTIONS FOR q=7 AND q=10")
print("="*60)
predictions = {}
formulas = {
    'quadratic': lambda q: np.polyval(coeffs_quad, q),
    'logarithmic': lambda q: log_func(q, *popt_log),
    'power_law': lambda q: power_func(q, *popt_pow),
    'minimal_model_cont': lambda q: c_from_m(m_analytic(q)),
    'coulomb_real_cont': coulomb_real_cont,
    'inv_power': lambda q: inv_power(q, *popt_inv),
}

for q_pred in [5, 7, 10]:
    print(f"\nq={q_pred}:")
    preds = {}
    for name, func in formulas.items():
        try:
            c = float(func(q_pred))
            print(f"  {name:25s}: c = {c:.4f}")
            preds[name] = c
        except:
            print(f"  {name:25s}: FAILED")
    predictions[q_pred] = preds

# Chi-squared for each formula against all 4 data points
print("\n=== Goodness of fit (chi^2 on q=2-5) ===")
for name, func in formulas.items():
    try:
        c_pred = np.array([func(q) for q in q_all])
        chi2 = np.sum((c_pred - c_all)**2 / np.array([0.001, 0.001, 0.001, 0.10])**2)
        print(f"  {name:25s}: chi^2 = {chi2:.2f}")
    except:
        pass

# Save results
results = {
    'experiment': '056a',
    'description': 'c(q) formula fitting and predictions',
    'data_used': {
        'q': [2, 3, 4, 5],
        'c_exact': [0.500, 0.800, 1.000, None],
        'c_measured': [None, None, None, 1.10],
        'c_error': [None, None, None, 0.10],
    },
    'formulas': {
        'quadratic': {'coeffs': coeffs_quad.tolist(), 'form': 'a*q^2 + b*q + c'},
        'logarithmic': {'params': popt_log.tolist(), 'form': 'a*ln(q) + b'},
        'power_law': {'params': popt_pow.tolist(), 'form': 'a*(q-1)^b'},
        'inv_power': {'params': popt_inv.tolist(), 'form': '1 - a/(q-1)^b'},
    },
    'predictions': {str(k): v for k, v in predictions.items()},
}
with open('results/exp_056a.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved results/exp_056a.json")
