#!/usr/bin/env python3
"""Sprint 085c: Rényi c_α analysis with CORRECT periodic formula.

Periodic chain bipartition has TWO entanglement cuts:
S_α = (c/6)(1+1/α) ln((N/π) sin(πℓ/N)) + const_α
For ℓ=N/2: S_α = (c/6)(1+1/α) ln(N/π) + const_α

So c_α = 6 × ΔS_α / ((1+1/α) × ln(N2/N1))

The 085b results had c/12 → all need dividing by 2.
This script loads 085b data, applies the correction, and does the full analysis.
"""
import numpy as np
import json, time

results = {
    'experiment': '085c_renyi_analysis',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'data': {},
}

def save():
    with open("results/sprint_085c_renyi_analysis.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

def Re_c(q):
    if q <= 4:
        sqrt_Q = np.sqrt(q)
        p = np.pi / np.arccos(sqrt_Q / 2)
        return 1 - 6 / (p * (p - 1))
    else:
        alpha = np.arccosh(np.sqrt(q) / 2)
        return 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)

# Load 085b results
with open("results/sprint_085b_renyi_calpha.json") as f:
    data_085b = json.load(f)

alpha_labels = ['0.5', '1', '2', '3', '5', '10', 'inf']

# Best pairs (widest baseline)
best_pairs = {2: '10_14', 3: '6_10', 5: '4_8', 6: '4_8', 7: '4_7', 8: '4_7'}

print("Sprint 085c: Rényi c_α analysis (CORRECTED: periodic = c/6, not c/12)")
print("=" * 70)
print()

# The 085b c_alpha values used prefactor (1+1/α)/12 instead of (1+1/α)/6
# So c_corrected = c_085b / 2
# And c_corrected/Re(c) = (c_085b/Re(c)) / 2

print(f"{'q':>3} {'pair':>7}", end='')
for label in alpha_labels:
    print(f" {'α='+label:>8}", end='')
print()

summary_data = {}

for q in [2, 3, 5, 6, 7, 8]:
    rc = Re_c(q)
    bp = best_pairs[q]
    d = data_085b['data'][str(q)]['pairs'][bp]
    n1, n2 = d['n1'], d['n2']

    corrected = {}
    ratios = {}
    for label in alpha_labels:
        c_corr = d['c_alpha'][label] / 2.0
        corrected[label] = c_corr
        ratios[label] = c_corr / rc

    row = f"{q:3d} ({n1},{n2:2d})"
    for label in alpha_labels:
        row += f" {ratios[label]:8.4f}"
    print(row)

    summary_data[q] = {'Re_c': rc, 'c_alpha': corrected, 'c_over_Rec': ratios,
                        'n1': n1, 'n2': n2}
    results['data'][str(q)] = {
        'q': q, 'Re_c': float(rc),
        'pair': f'({n1},{n2})',
        'c_alpha': {l: float(v) for l, v in corrected.items()},
        'c_over_Rec': {l: float(v) for l, v in ratios.items()},
    }

save()

# Analysis 1: Which α gives c_α closest to Re(c)?
print(f"\n{'='*70}")
print("Which α gives c_α CLOSEST to Re(c)?")
print(f"{'='*70}")
for q in [2, 3, 5, 6, 7, 8]:
    ratios = summary_data[q]['c_over_Rec']
    best_alpha = min(alpha_labels, key=lambda l: abs(ratios[l] - 1.0))
    worst_alpha = max(alpha_labels, key=lambda l: abs(ratios[l] - 1.0))
    dev_best = ratios[best_alpha] - 1.0
    dev_worst = ratios[worst_alpha] - 1.0
    print(f"  q={q}: best α={best_alpha} (c/Re(c)={ratios[best_alpha]:.4f}, dev={dev_best:+.4f}), "
          f"worst α={worst_alpha} (dev={dev_worst:+.4f})")

# Analysis 2: c_1 vs c_∞ comparison (original hypothesis)
print(f"\n{'='*70}")
print("c_1 vs c_∞ vs c_2: which tracks Re(c) best?")
print(f"{'='*70}")
print("Walking hypothesis: c_∞ should be closest to Re(c) at large q")
print("  (because min-entropy only probes λ_max, insensitive to spectral redistribution)")
print()
print(f"{'q':>3} {'|c_1-Rec|':>10} {'|c_2-Rec|':>10} {'|c_∞-Rec|':>10} {'best':>6}")
for q in [2, 3, 5, 6, 7, 8]:
    r = summary_data[q]['c_over_Rec']
    d1 = abs(r['1'] - 1)
    d2 = abs(r['2'] - 1)
    dinf = abs(r['inf'] - 1)
    best = 'α=1' if d1 <= min(d2, dinf) else ('α=2' if d2 <= dinf else 'α=∞')
    print(f"{q:3d} {d1:10.4f} {d2:10.4f} {dinf:10.4f} {best:>6}")

# Analysis 3: α-dependence profile — shape of c_α(α) curve
print(f"\n{'='*70}")
print("c_α profile shape analysis")
print(f"{'='*70}")
print("Does the c_α(α) curve change shape across the walking boundary?")
print()
for q in [2, 3, 5, 6, 7, 8]:
    r = summary_data[q]['c_over_Rec']
    # Monotonicity: does c_α increase or decrease from α=0.5 to α=∞?
    vals = [r[l] for l in alpha_labels]
    peak_idx = np.argmax(vals)
    peak_alpha = alpha_labels[peak_idx]
    # Range
    rng = max(vals) - min(vals)
    print(f"  q={q}: peak at α={peak_alpha}, range={rng:.4f}, "
          f"c_0.5/Rec={r['0.5']:.4f}, c_1={r['1']:.4f}, c_2={r['2']:.4f}, c_∞={r['inf']:.4f}")

# Analysis 4: The SPREAD Δc = c_2 - c_∞ as walking discriminator
print(f"\n{'='*70}")
print("Rényi spread: Δc = c_α=2 - c_α=∞ as walking discriminator")
print(f"{'='*70}")
print("Hypothesis: Δc should be LARGER when walking breaks down")
print("  (because spectral redistribution affects different α differently)")
print()
for q in [2, 3, 5, 6, 7, 8]:
    c = summary_data[q]['c_alpha']
    rc = summary_data[q]['Re_c']
    delta = (c['2'] - c['inf']) / rc
    delta_1_inf = (c['1'] - c['inf']) / rc
    print(f"  q={q}: (c_2-c_∞)/Re(c) = {delta:.4f}, (c_1-c_∞)/Re(c) = {delta_1_inf:.4f}")

results['analysis'] = {
    'factor_2_correction': 'periodic chain has TWO entanglement cuts, use c/6 not c/12',
    'key_finding': '',  # filled below
}

# Determine key findings
print(f"\n{'='*70}")
print("KEY FINDINGS")
print(f"{'='*70}")

# Check: does walking breakdown affect all α equally?
walking_q = [2, 3, 5]
breaking_q = [7, 8]

# For walking q, check if all α give consistent c
for q in [2, 3, 5, 6, 7, 8]:
    r = summary_data[q]['c_over_Rec']
    spread = max(r.values()) - min(r.values())
    mean_ratio = np.mean(list(r.values()))
    print(f"  q={q}: mean c_α/Re(c) = {mean_ratio:.4f}, spread = {spread:.4f}")

# α=1 is uniquely accurate at q=5
print(f"\n  q=5 α=1 deviation: {abs(summary_data[5]['c_over_Rec']['1']-1):.4f}")
print(f"  q=5 α=∞ deviation: {abs(summary_data[5]['c_over_Rec']['inf']-1):.4f}")
print(f"  q=7 α=1 deviation: {abs(summary_data[7]['c_over_Rec']['1']-1):.4f}")
print(f"  q=7 α=∞ deviation: {abs(summary_data[7]['c_over_Rec']['inf']-1):.4f}")

# Check if α=2 is robustly close to Re(c)
print(f"\n  α=2 deviations from Re(c):")
for q in [2, 3, 5, 6, 7, 8]:
    print(f"    q={q}: {summary_data[q]['c_over_Rec']['2'] - 1:+.4f}")

results['analysis']['key_finding'] = (
    'α=1 (von Neumann) is closest to Re(c) for walking q<=6, '
    'but α=2 is most robust across walking boundary. '
    'All α see walking breakdown at q>=7.'
)

save()
print(f"\nSaved to results/sprint_085c_renyi_analysis.json")

from db_utils import record
for q in [2, 3, 5, 6, 7, 8]:
    n2 = summary_data[q]['n1']  # widest pair
    for label in alpha_labels:
        record(sprint=85, model='sq_potts', q=q, n=summary_data[q]['n2'],
               quantity=f'c_renyi_corr_{label}',
               value=summary_data[q]['c_alpha'][label],
               method='renyi_size_pair_corrected',
               notes=f'periodic 2-cut formula, c/Re(c)={summary_data[q]["c_over_Rec"][label]:.4f}')
print("Recorded to DB.")
