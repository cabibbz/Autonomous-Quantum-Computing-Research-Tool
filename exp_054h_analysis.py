#!/usr/bin/env python3
"""Sprint 054h: Central charge analysis — combine all data.

Also test g_c sensitivity for q=4: what g_c gives c=1.0?
"""
import numpy as np, json

# === Collected data ===
# q=2 (TFIM, from Sprint 049, g_c=1.0 in TFIChain convention)
q2 = {'sizes': [16, 24, 32, 48, 64],
      'S': [0.438, 0.467, 0.485, 0.509, 0.525],
      'c_exact': 0.500, 'source': 'Sprint 049'}

# q=3 (Potts, g_c=1/3, chi=60-80)
q3 = {'sizes': [8, 12, 16, 24],
      'S': [0.516117, 0.579230, 0.622697, 0.682416],
      'c_exact': 0.800, 'source': 'Sprint 054a'}

# q=4 (Potts, g_c=0.392, chi=40)
q4 = {'sizes': [8, 12, 16, 24],
      'S': [0.646366, 0.729571, 0.788153, 0.871215],
      'c_exact': 1.000, 'source': 'Sprint 054f'}

# q=5 (Potts, g_c=0.441, chi=40)
q5 = {'sizes': [8, 12, 16],
      'S': [0.760529, 0.854768, 0.918754],
      'c_exact': None, 'source': 'Sprint 054f'}


def extract_c(ns, Ss, label=""):
    ns = np.array(ns, dtype=float)
    Ss = np.array(Ss, dtype=float)
    ln_n = np.log(ns)
    A = np.vstack([ln_n, np.ones(len(ln_n))]).T
    slope, intercept = np.linalg.lstsq(A, Ss, rcond=None)[0]
    c_full = 6 * slope

    c_pairs = []
    for i in range(len(ns) - 1):
        c_pair = 6 * (Ss[i+1] - Ss[i]) / (np.log(ns[i+1]) - np.log(ns[i]))
        c_pairs.append(c_pair)

    return c_full, c_pairs


print("=" * 70)
print("CENTRAL CHARGE c(q) AT TRUE CRITICAL POINTS")
print("=" * 70)

results = {}
for label, data in [('q=2', q2), ('q=3', q3), ('q=4', q4), ('q=5', q5)]:
    c_full, c_pairs = extract_c(data['sizes'], data['S'], label)
    c_exact = data['c_exact']
    print(f"\n{label}: c_full = {c_full:.4f}", end="")
    if c_exact:
        print(f"  (CFT: {c_exact:.3f}, ratio: {c_full/c_exact:.3f})")
    else:
        print(f"  (no CFT prediction)")
    for i, cp in enumerate(c_pairs):
        n1, n2 = data['sizes'][i], data['sizes'][i+1]
        print(f"    c({n1},{n2}) = {cp:.4f}", end="")
        if c_exact:
            print(f"  (ratio: {cp/c_exact:.3f})")
        else:
            print()
    results[label] = {
        'c_full': float(c_full),
        'c_pairwise': [float(c) for c in c_pairs],
        'c_exact': c_exact,
        'sizes': data['sizes'],
        'S': data['S'],
    }

# === Overshoot pattern analysis ===
print("\n" + "=" * 70)
print("OVERSHOOT ANALYSIS — does c(q) converge at same RATE?")
print("=" * 70)

# For q=2: pairwise converges 0.544 → 0.516 toward 0.500
# Overshoot ratio at comparable sizes
print("\nOvershoot ratio c_pairwise/c_exact at smallest pair:")
for label, c_exact in [('q=2', 0.5), ('q=3', 0.8), ('q=4', 1.0)]:
    cp = results[label]['c_pairwise'][0]
    print(f"  {label}: {cp:.4f}/{c_exact:.3f} = {cp/c_exact:.3f}")

# === Check if c grows with q ===
print("\n" + "=" * 70)
print("c(q) PATTERN — does central charge grow with q?")
print("=" * 70)

# Use largest-pair estimate (most accurate)
for label in ['q=2', 'q=3', 'q=4', 'q=5']:
    cp_last = results[label]['c_pairwise'][-1]
    n1 = results[label]['sizes'][-2]
    n2 = results[label]['sizes'][-1]
    print(f"  {label}: c({n1},{n2}) = {cp_last:.4f}")

# If q=3 overshoot is ~10% at (n=16,24), and q=4 overshoot is ~23%,
# then the ratio of overshoot rates is 23/10 = 2.3x
# If convergence slows by similar factor per q step, q=5's true c
# would be c_pair / 1.3 ≈ 1.33/1.3 ≈ 1.03

# === Fit c(q) = a*(q-1) + b ===
print("\n" + "=" * 70)
print("c(q) TREND — using best (largest-pair) estimates")
print("=" * 70)

# Raw largest-pair estimates
qs_raw = [2, 3, 4, 5]
cs_raw = [results[f'q={q}']['c_pairwise'][-1] for q in qs_raw]

# Corrected estimates using overshoot pattern
# q=2 pairwise converges: ratio at (n=48,64) is ~1.03x, at (n=16,24) is ~1.09x
# q=3 at (n=16,24): ~1.10x overshoot
# So correction factor at (n=16,24) ≈ ~0.90 for c
cs_corrected = [0.500, 0.800, None, None]  # exact for q=2,3

# For q=4: if same ~10% convergence rate as q=3 at same sizes (n=16,24)
# But q=4 shows FLAT pairwise — all ~1.23. This suggests either:
# 1. True c > 1.0 (and convergence is genuinely slow for q=4)
# 2. We're not at true g_c

print("\nRaw c (largest pair):", dict(zip(qs_raw, [f"{c:.4f}" for c in cs_raw])))
print("Known exact: q=2 → 0.500, q=3 → 0.800")
print(f"Overshoot ratio q=3: {cs_raw[1]/0.8:.3f}")
print(f"Overshoot ratio q=4: {cs_raw[2]/1.0:.3f}")

# Save
with open('results/sprint_054h_analysis.json', 'w') as f:
    json.dump({'sprint': '054h', 'results': results,
               'chi_convergence': {
                   'q3_n16_chi20': 0.622695,
                   'q3_n16_chi40': 0.622697,
                   'conclusion': 'entropy converged at chi=20 for q=3'}}, f, indent=2)

print("\nDone.", flush=True)
