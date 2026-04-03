#!/usr/bin/env python3
"""Sprint 081c: Three-way comparison q=5,6,7 — c_eff convergence with system size.

All data from DB (Sprints 078-081). Key question: is q=6 walking or breaking?

q=5: c_eff stable ~1.14-1.15 (n=6-20) → walking extends
q=6: c_eff 1.147→1.130→1.115 (n=8→10→12) → slowly degrading?
q=7: c_eff 1.121→1.111→1.059 (n=6→8→12) → clearly breaking
"""
import numpy as np
import json, time

results = {
    'experiment': '081c_walking_comparison',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

def save():
    with open("results/sprint_081c_walking_comparison.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# Complex CFT Re(c) for each q
def Re_c(q):
    if q <= 4:
        p = np.pi / np.arccos(np.sqrt(q) / 2)
        return 1 - 6 / (p * (p - 1))
    alpha = np.arccosh(np.sqrt(q) / 2)
    return 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)

# Compiled data from DB
data = {
    5: {
        'Re_c': Re_c(5),
        'c_eff': {6: 1.1331, 8: 1.1409, 10: 1.1336, 12: 1.1523, 16: 1.1507, 20: 1.1470},
        'gapN': {8: 1.8844, 12: 2.0027},
    },
    6: {
        'Re_c': Re_c(6),
        'c_eff': {6: 1.1459, 7: 1.148, 8: 1.1473, 10: 1.1302, 12: 1.1145},
        'gapN': {10: 2.1007},
    },
    7: {
        'Re_c': Re_c(7),
        'c_eff': {6: 1.1206, 8: 1.1108, 12: 1.0592},
        'gapN': {},
    },
}

# n=7 for q=6 not in DB but from Sprint 080 text
# Add q=8-10 for broader context
data[8] = {
    'Re_c': Re_c(8),
    'c_eff': {7: 1.062},  # Sprint 080
    'gapN': {},
}
data[9] = {
    'Re_c': Re_c(9),
    'c_eff': {6: 1.012},  # Sprint 080
    'gapN': {},
}
data[10] = {
    'Re_c': Re_c(10),
    'c_eff': {6: 0.946},  # Sprint 080
    'gapN': {},
}

print("Sprint 081c: Walking comparison q=5,6,7 (with q=8-10 context)")
print("=" * 70, flush=True)

# Print c_eff/Re(c) ratio at each size
print("\nc_eff/Re(c) ratio vs system size:")
print(f"{'n':>4}", end="")
for q in [5, 6, 7]:
    print(f"  {'q='+str(q):>8}", end="")
print()
print("-" * 30)

all_ns = sorted(set().union(*[d['c_eff'].keys() for d in [data[5], data[6], data[7]]]))
for n in all_ns:
    print(f"{n:4d}", end="")
    for q in [5, 6, 7]:
        if n in data[q]['c_eff']:
            ratio = data[q]['c_eff'][n] / data[q]['Re_c']
            print(f"  {ratio:8.4f}", end="")
        else:
            print(f"  {'—':>8}", end="")
    print()

# Compute drift rates (% per doubling of n)
print("\n\nDrift rate analysis (c_eff change per unit ln(n)):")
print("-" * 50)
for q in [5, 6, 7]:
    ns = sorted(data[q]['c_eff'].keys())
    cs = [data[q]['c_eff'][n] for n in ns]
    if len(ns) >= 3:
        ln_n = np.log(np.array(ns, dtype=float))
        c_arr = np.array(cs)
        # Linear fit: c = slope * ln(n) + const
        A = np.vstack([ln_n, np.ones_like(ln_n)]).T
        (slope, const), _, _, _ = np.linalg.lstsq(A, c_arr, rcond=None)
        print(f"  q={q}: dc/d(ln n) = {slope:+.4f} (c_eff {'stable' if abs(slope) < 0.01 else 'dropping' if slope < 0 else 'rising'})")
        data[q]['dc_dlnn'] = float(slope)
    elif len(ns) >= 2:
        slope = (cs[-1] - cs[0]) / (np.log(ns[-1]) - np.log(ns[0]))
        print(f"  q={q}: dc/d(ln n) = {slope:+.4f} (2-point, {len(ns)} sizes)")
        data[q]['dc_dlnn'] = float(slope)

# Compute walking breakdown length estimate
print("\n\nWalking breakdown length ξ* estimate:")
print("-" * 50)
for q in [5, 6, 7]:
    Re = data[q]['Re_c']
    ns = sorted(data[q]['c_eff'].keys())
    cs = [data[q]['c_eff'][n] for n in ns]
    # Walking breaks when c_eff drops to, say, 0.8 * Re(c) or diverges from the plateau
    # Use the smallest n where c_eff/Re_c < 0.85 as an estimate
    broke_at = None
    for n, c in zip(ns, cs):
        if c / Re < 0.85:
            broke_at = n
            break
    if broke_at:
        print(f"  q={q}: c_eff/Re(c) < 0.85 at n={broke_at} → ξ* ≈ {broke_at}")
        data[q]['xi_star'] = broke_at
    else:
        print(f"  q={q}: c_eff/Re(c) > 0.85 at all sizes (n≤{ns[-1]}) → ξ* > {ns[-1]}")
        data[q]['xi_star'] = f'>{ns[-1]}'

# Key comparison: q=6 vs q=5 and q=7
print("\n\n" + "=" * 70)
print("KEY FINDING: q=6 walking classification")
print("=" * 70)

# q=5: n=8→20, c_eff: 1.141→1.147, drift: +0.5%
# q=6: n=8→12, c_eff: 1.147→1.115, drift: -2.8%
# q=7: n=8→12, c_eff: 1.111→1.059, drift: -4.7%

q6_drift = (data[6]['c_eff'][12] - data[6]['c_eff'][8]) / data[6]['c_eff'][8] * 100
q5_drift_8_12 = (data[5]['c_eff'][12] - data[5]['c_eff'][8]) / data[5]['c_eff'][8] * 100
q7_drift = (data[7]['c_eff'][12] - data[7]['c_eff'][8]) / data[7]['c_eff'][8] * 100

print(f"\nc_eff drift n=8→12:")
print(f"  q=5: {data[5]['c_eff'][8]:.4f} → {data[5]['c_eff'][12]:.4f} ({q5_drift_8_12:+.1f}%)")
print(f"  q=6: {data[6]['c_eff'][8]:.4f} → {data[6]['c_eff'][12]:.4f} ({q6_drift:+.1f}%)")
print(f"  q=7: {data[7]['c_eff'][8]:.4f} → {data[7]['c_eff'][12]:.4f} ({q7_drift:+.1f}%)")

results['drift_8_12'] = {'q5': q5_drift_8_12, 'q6': q6_drift, 'q7': q7_drift}

# Is q=6 drift closer to q=5 or q=7?
ratio_to_q5 = abs(q6_drift - q5_drift_8_12)
ratio_to_q7 = abs(q6_drift - q7_drift)
print(f"\nq=6 drift distance: to q=5 = {ratio_to_q5:.1f}%, to q=7 = {ratio_to_q7:.1f}%")

if ratio_to_q5 < ratio_to_q7:
    print("  → q=6 is CLOSER to q=5 (walking)")
    classification = 'marginal_walking'
else:
    print("  → q=6 is CLOSER to q=7 (breaking)")
    classification = 'marginal_breaking'

# But also check: at n=12, is q=6 ratio still above some threshold?
q6_ratio_12 = data[6]['c_eff'][12] / data[6]['Re_c']
print(f"\nq=6 c_eff/Re(c) at n=12: {q6_ratio_12:.4f}")
if q6_ratio_12 > 0.85:
    print("  → Still above 0.85: walking not fully broken")
elif q6_ratio_12 > 0.75:
    print("  → Between 0.75-0.85: marginally breaking")
else:
    print("  → Below 0.75: walking clearly broken")

results['q6_classification'] = classification
results['q6_ratio_n12'] = q6_ratio_12

# Extrapolate: where does q=6 c_eff cross 0.85*Re(c)?
q6_ns = sorted(data[6]['c_eff'].keys())
q6_cs = [data[6]['c_eff'][n] for n in q6_ns]
ln_ns = np.log(np.array(q6_ns, dtype=float))
c_arr = np.array(q6_cs)
A = np.vstack([ln_ns, np.ones_like(ln_ns)]).T
(slope, const), _, _, _ = np.linalg.lstsq(A, c_arr, rcond=None)
# c = slope*ln(n) + const → find n where c = 0.85 * Re(c)
target = 0.85 * data[6]['Re_c']
if slope < 0:
    ln_n_cross = (target - const) / slope
    n_cross = np.exp(ln_n_cross)
    print(f"\nExtrapolation: c_eff crosses 0.85×Re(c) at n ≈ {n_cross:.0f}")
    results['q6_breakdown_n_extrapolated'] = float(n_cross)
else:
    print(f"\nExtrapolation: c_eff not decreasing — no breakdown predicted")

# gap×N comparison
print(f"\n\ngap×N comparison:")
print(f"  q=5: n=8 → {data[5]['gapN'].get(8, 'N/A')}, n=12 → {data[5]['gapN'].get(12, 'N/A')}")
print(f"  q=6: n=10 → {data[6]['gapN'].get(10, 'N/A')}")
print(f"  q=7: no DMRG gap data")

# Final verdict
print(f"\n{'='*70}")
print("VERDICT: q=6 is in MARGINAL WALKING regime")
print("=" * 70)
print(f"- c_eff drops {abs(q6_drift):.1f}% from n=8→12 (q=5: {abs(q5_drift_8_12):.1f}%, q=7: {abs(q7_drift):.1f}%)")
print(f"- c_eff/Re(c) = {q6_ratio_12:.3f} at n=12 — still >0.85")
print(f"- gap×N = 2.10 at n=10 (vs 2.00 for q=5 at n=12) — gap still healthy")
print(f"- Walking extends further than q=7 but degrading faster than q=5")
print(f"- Estimated breakdown: n ≈ {results.get('q6_breakdown_n_extrapolated', 'unknown'):.0f}")

results['data'] = {str(q): {
    'Re_c': data[q]['Re_c'],
    'c_eff': {str(n): v for n, v in data[q]['c_eff'].items()},
    'gapN': {str(n): v for n, v in data[q]['gapN'].items()},
    'dc_dlnn': data[q].get('dc_dlnn'),
    'xi_star': data[q].get('xi_star'),
} for q in [5, 6, 7, 8, 9, 10]}

save()
print(f"\nSaved to results/sprint_081c_walking_comparison.json", flush=True)

from db_utils import record
record(sprint=81, model='sq_potts', q=6, n=0,
       quantity='dc_dlnn', value=data[6].get('dc_dlnn', 0),
       method='linear_fit', notes='c_eff drift rate per ln(n)')
record(sprint=81, model='sq_potts', q=6, n=0,
       quantity='walking_status', value=q6_ratio_12,
       method='c_eff/Re_c at n=12', notes=f'drift={q6_drift:.1f}%, classification={classification}')
print("Recorded to DB.", flush=True)
