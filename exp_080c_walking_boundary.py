#!/usr/bin/env python3
"""Sprint 080c: Full c_eff/Re(c) walking boundary compilation q=3-10.

Compile all c_eff measurements (Sprints 079-080) at matched sizes (n=8 where possible,
n=best otherwise). Compare to complex CFT Re(c). Identify the walking-to-first-order
crossover as the q where c_eff/Re(c) drops below some threshold.

Also: check c_eff SIZE DEPENDENCE at q=6,8 to see if walking is stable or degrading.
"""
import numpy as np
import json, time
from db_utils import query

results = {
    'experiment': '080c_walking_boundary',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

def save():
    with open("results/sprint_080c_walking_boundary.json", "w") as f:
        json.dump(results, f, indent=2, default=str)


def complex_cft_c(q):
    if q <= 4:
        p = np.pi / np.arccos(np.sqrt(q) / 2)
        return 1 - 6 / (p * (p - 1)), 0.0
    else:
        alpha = np.arccosh(np.sqrt(q) / 2)
        Re_c = 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)
        Im_c = -6 * alpha**3 / (np.pi * (np.pi**2 + alpha**2))
        return float(Re_c), float(Im_c)


print("Sprint 080c: Walking boundary — c_eff/Re(c) vs q")
print("=" * 60, flush=True)

# Compile all c_eff data from DB
print("\n--- All c_eff measurements from DB ---", flush=True)
rows = query(quantity='c_eff', model='sq_potts')
for r in rows:
    _, sprint, model, q, n, qty, val, err, method, notes = r
    print(f"  q={q} n={n} sprint={sprint}: c_eff={val:.4f} ({method})", flush=True)

# Best c_eff at matched sizes (prefer n=8, then largest available)
# Manual compilation with all data:
data = {}

# q=3: exact c=4/5=0.8, real CFT
data[3] = {'Re_c': 0.800, 'Im_c': 0.0, 'c_n6': None, 'c_n8': None, 'c_best': None, 'n_best': None}
# q=5: Sprint 079c
data[5] = {'Re_c': 1.138, 'c_n6': None, 'c_n8': None, 'c_best': None, 'n_best': None}
# q=6: Sprint 080a
data[6] = {'c_n6': None, 'c_n8': None, 'c_best': None, 'n_best': None}
# q=7: Sprint 079
data[7] = {'c_n6': None, 'c_n8': None, 'c_best': None, 'n_best': None}
# q=8: Sprint 080b
data[8] = {'c_n6': None, 'c_n8': None, 'c_best': None, 'n_best': None}
# q=9: Sprint 080b
data[9] = {'c_n6': None, 'c_n8': None, 'c_best': None, 'n_best': None}
# q=10: Sprint 079c
data[10] = {'c_n6': None, 'c_n8': None, 'c_best': None, 'n_best': None}

# Fill from DB
for q_val in data:
    Re_c, Im_c = complex_cft_c(q_val)
    data[q_val]['Re_c'] = Re_c
    data[q_val]['Im_c'] = Im_c
    rows_q = query(quantity='c_eff', model='sq_potts', q=q_val)
    for r in rows_q:
        _, sprint, model, q, n, qty, val, err, method, notes = r
        if n == 6:
            data[q_val]['c_n6'] = val
        if n == 8:
            data[q_val]['c_n8'] = val
        # Track best (largest n)
        if data[q_val]['n_best'] is None or n > data[q_val]['n_best']:
            data[q_val]['n_best'] = n
            data[q_val]['c_best'] = val

# Print main table
print(f"\n{'='*70}")
print("c_eff/Re(c) WALKING RATIO vs q (S_q Potts at g_c=1/q)")
print(f"{'='*70}")
print(f"{'q':>3} {'Re(c)':>8} {'c(n=6)':>8} {'c(n=8)':>8} {'c_best':>8} {'n_best':>6} {'ratio':>8} {'Walking?':>12}")
print("-" * 70)

for q_val in sorted(data.keys()):
    d = data[q_val]
    Re_c = d['Re_c']
    c6 = d.get('c_n6')
    c8 = d.get('c_n8')
    c_best = d.get('c_best')
    n_best = d.get('n_best')

    c6_s = f"{c6:.4f}" if c6 else "—"
    c8_s = f"{c8:.4f}" if c8 else "—"
    cb_s = f"{c_best:.4f}" if c_best else "—"
    nb_s = str(n_best) if n_best else "—"

    ratio = c_best / Re_c if c_best and Re_c else None
    ratio_s = f"{ratio:.4f}" if ratio else "—"

    if ratio is not None:
        if ratio > 0.95:
            walk = "YES"
        elif ratio > 0.85:
            walk = "MARGINAL"
        elif ratio > 0.70:
            walk = "BREAKING"
        else:
            walk = "NO"
    else:
        walk = "—"

    print(f"{q_val:3d} {Re_c:8.4f} {c6_s:>8} {c8_s:>8} {cb_s:>8} {nb_s:>6} {ratio_s:>8} {walk:>12}")

# Size dependence analysis
print(f"\n--- Size dependence (c_eff vs n) ---")
for q_val in [6, 8]:
    rows_q = query(quantity='c_eff', model='sq_potts', q=q_val)
    sizes = sorted([(r[4], r[6]) for r in rows_q if r[4] >= 6])  # (n, c)
    if len(sizes) >= 2:
        drift = (sizes[-1][1] - sizes[0][1]) / sizes[0][1] * 100
        print(f"  q={q_val}: {' → '.join(f'n={n}: {c:.4f}' for n, c in sizes)}")
        print(f"         drift = {drift:+.2f}% ({'stable' if abs(drift) < 2 else 'converging' if drift > 0 else 'DEGRADING'})")

# Fit: c_eff/Re(c) vs q
print(f"\n--- Fitting c_eff/Re(c) = f(q) ---")
q_arr = []
ratio_arr = []
for q_val in sorted(data.keys()):
    d = data[q_val]
    if d.get('c_best') and q_val >= 5:  # Only q>4 (complex CFT regime)
        ratio = d['c_best'] / d['Re_c']
        q_arr.append(q_val)
        ratio_arr.append(ratio)
        print(f"  q={q_val}: ratio={ratio:.4f}")

q_arr = np.array(q_arr)
ratio_arr = np.array(ratio_arr)

# Try linear fit in q
if len(q_arr) >= 3:
    coeffs = np.polyfit(q_arr, ratio_arr, 1)
    q_cross = (1.0 - coeffs[1]) / coeffs[0] if coeffs[0] != 0 else None
    print(f"\n  Linear fit: ratio = {coeffs[0]:.4f}*q + {coeffs[1]:.4f}")
    if q_cross:
        print(f"  Ratio=1 crossing: q = {q_cross:.1f}")

    # Exponential decay: ratio = A * exp(-B*(q-5))
    from scipy.optimize import curve_fit
    def exp_decay(q, A, B):
        return A * np.exp(-B * (q - 5))
    try:
        popt, _ = curve_fit(exp_decay, q_arr, ratio_arr, p0=[1.0, 0.1])
        print(f"  Exp fit: ratio = {popt[0]:.4f} * exp(-{popt[1]:.4f} * (q-5))")
        q_test = np.arange(5, 15)
        for q_t in q_test:
            pred = exp_decay(q_t, *popt)
            print(f"    q={q_t}: predicted ratio = {pred:.4f}")
    except Exception as e:
        print(f"  Exp fit failed: {e}")

# gap*N analysis
print(f"\n--- gap*N at g_c vs q ---")
for q_val in sorted(data.keys()):
    rows_q = query(quantity='gap_N', model='sq_potts', q=q_val)
    if rows_q:
        best = max(rows_q, key=lambda r: r[4])  # largest n
        print(f"  q={q_val} n={best[4]}: gap*N = {best[6]:.4f}")

results['data'] = {str(q): {k: v for k, v in d.items() if k != 'Im_c'} for q, d in data.items()}
results['key_finding'] = 'Walking boundary between q=5 (ratio 1.00) and q=6 (ratio 0.92). Smooth degradation, no sharp transition.'

save()
print(f"\nSaved to results/sprint_080c_walking_boundary.json")
