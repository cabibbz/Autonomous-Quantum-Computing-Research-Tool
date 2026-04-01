"""
Sprint 041d: Cross-q MI-CV analysis.
Compile q=2,3,4,5 data and compare crossing points, slopes, CV patterns.
"""

import json, os
import numpy as np

print("=== Sprint 041d: Cross-q MI-CV Analysis ===\n")

# q=5 data (this sprint)
q5_n8 = {
    0.50: 0.1409, 0.80: 0.3616, 0.90: 0.4384, 0.95: 0.4790,
    1.00: 0.5213, 1.05: 0.5651, 1.10: 0.6103, 1.30: 0.7941
}
q5_n12 = {0.50: 0.1206, 0.80: 0.3766, 1.00: 0.5794, 1.10: 0.7017}

# Load q=4 data
try:
    with open('results/sprint_040a2_potts_q4_n8.json') as f:
        d4_n8 = json.load(f)
    q4_n8 = {float(k): v['cv'] for k, v in d4_n8['data'].items() if v.get('cv') is not None}
except: q4_n8 = {}

try:
    with open('results/sprint_040b_potts_q4_n12.json') as f:
        d4_n12 = json.load(f)
    q4_n12 = {float(k): v['cv'] for k, v in d4_n12['data'].items() if v.get('cv') is not None}
except: q4_n12 = {}

# Also try fill data
for fname in ['results/sprint_040b2_potts_q4_n12_fill.json']:
    try:
        with open(fname) as f:
            d = json.load(f)
        for k, v in d['data'].items():
            if v.get('cv') is not None:
                q4_n12[float(k)] = v['cv']
    except: pass

print("=== q=5 n=8 vs n=12 crossing detection ===")
common_g = sorted(set(q5_n8.keys()) & set(q5_n12.keys()))
for g in common_g:
    diff = q5_n12[g] - q5_n8[g]
    sign = "↑" if diff > 0 else "↓"
    print(f"  g={g:.2f}: n=8 CV={q5_n8[g]:.4f}, n=12 CV={q5_n12[g]:.4f}, diff={diff:+.4f} {sign}")

# Find crossing point by interpolation
g_vals = sorted(common_g)
diffs = [q5_n12[g] - q5_n8[g] for g in g_vals]
print(f"\n  Diffs: {[f'{d:+.4f}' for d in diffs]}")
for i in range(len(diffs)-1):
    if diffs[i] * diffs[i+1] < 0:
        g_cross = g_vals[i] + (-diffs[i]) * (g_vals[i+1] - g_vals[i]) / (diffs[i+1] - diffs[i])
        print(f"  CROSSING at g ≈ {g_cross:.3f} (between g={g_vals[i]:.2f} and g={g_vals[i+1]:.2f})")

print(f"\n=== q=4 data ===")
print(f"  n=8: {sorted(q4_n8.items())}")
print(f"  n=12: {sorted(q4_n12.items())}")

print(f"\n=== Crossing point trend with q ===")
# Known crossings from prior sprints
crossings = {
    'q=2 (TFIM)': {'g_c_n8_12': 0.93, 'note': 'Sprint 037'},
    'q=3 (Potts)': {'g_c_n8_12': 0.923, 'note': 'Sprint 039'},
    'q=4 (clock)': {'g_c_n8_12': 0.893, 'note': 'Sprint 040'},
}
for name, info in crossings.items():
    print(f"  {name}: g_c ≈ {info['g_c_n8_12']:.3f} ({info['note']})")

# q=5 crossing interpolation
g_cross_q5 = None
for i in range(len(diffs)-1):
    if diffs[i] * diffs[i+1] < 0:
        g_cross_q5 = g_vals[i] + (-diffs[i]) * (g_vals[i+1] - g_vals[i]) / (diffs[i+1] - diffs[i])
print(f"  q=5 (clock): g_c ≈ {g_cross_q5:.3f} (this sprint)")

print(f"\n=== Slope at g=1.0 (transition region) ===")
# Slope = CV change per unit g around g=1.0
# For q=5 n=8: use points around g=1.0
q5_n8_slope = (q5_n8[1.05] - q5_n8[0.95]) / 0.1
print(f"  q=5 n=8: slope ≈ {q5_n8_slope:.2f}")
if 0.95 in q4_n8 and 1.05 in q4_n8:
    q4_n8_slope = (q4_n8[1.05] - q4_n8[0.95]) / 0.1
    print(f"  q=4 n=8: slope ≈ {q4_n8_slope:.2f}")

print(f"\n=== CV magnitude comparison at g=1.0 ===")
if 1.0 in q4_n8:
    print(f"  q=4 n=8: CV={q4_n8[1.0]:.4f}")
print(f"  q=5 n=8: CV={q5_n8[1.0]:.4f}")
if 1.0 in q4_n12:
    print(f"  q=4 n=12: CV={q4_n12[1.0]:.4f}")
print(f"  q=5 n=12: CV={q5_n12[1.0]:.4f}")

# Save analysis
analysis = {
    'experiment': '041d',
    'q5_n8': q5_n8,
    'q5_n12': q5_n12,
    'crossing_q5': g_cross_q5,
    'crossing_trend': {'q2': 0.93, 'q3': 0.923, 'q4': 0.893, 'q5': g_cross_q5},
    'cv_at_g1': {
        'q5_n8': q5_n8[1.0],
        'q5_n12': q5_n12[1.0],
    }
}
with open('results/sprint_041d_cross_q_analysis.json', 'w') as f:
    json.dump(analysis, f, indent=2, default=str)
print("\nResults saved.")
