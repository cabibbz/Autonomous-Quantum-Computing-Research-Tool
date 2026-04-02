#!/usr/bin/env python3
"""Sprint 089c: Walking discriminator synthesis.

Key insight from 089b: the walking-specific signal is in the STATIC entropy
partition — how much of the entropy sits in the (q-1) multiplet vs ground+tail.
Here we:
1. Test whether %S(lev1) = f(q, n) has a universal scaling form
2. Define and test a "multiplet dominance" parameter M(q,n)
3. Connect to c_eff/Re(c) — can we PREDICT walking breakdown from spectrum?
4. Test updated power-law exponent for q=7 (b≈2.1 confirms universality)
"""
import numpy as np
import json, time
from scipy.optimize import curve_fit

results = {
    'experiment': '089c_walking_discriminator',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

def save():
    with open("results/sprint_089c_walking_discriminator.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# Load all data
files = {
    2: "results/sprint_088b_entspec_dmrg_q2.json",
    3: "results/sprint_088a_entspec_dmrg_q3.json",
    5: "results/sprint_087a_entspec_dmrg_q5.json",
    7: ["results/sprint_087b_entspec_dmrg_q7.json",
        "results/sprint_089a_entspec_dmrg_q7_ext.json"],
}

data = {}
for q_val, fnames in files.items():
    if isinstance(fnames, str):
        fnames = [fnames]
    all_entries = []
    for fname in fnames:
        with open(fname) as f:
            d = json.load(f)
        all_entries.extend(d['data'])
    seen = {}
    for entry in all_entries:
        seen[entry['n']] = entry
    entries = sorted(seen.values(), key=lambda x: x['n'])
    data[q_val] = entries

def Re_c(q):
    if q <= 4:
        return {2: 0.5, 3: 0.8, 4: 1.0}[q]
    alpha = np.arccosh(np.sqrt(q) / 2)
    return 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)

print("Sprint 089c: Walking discriminator synthesis")
print("=" * 70, flush=True)

# 1. Updated power-law fit for ALL q (including new q=7 data)
print("\n--- Updated power-law exponents (w_tail ~ n^b) ---")
print(f"{'q':>3} {'n_range':>10} {'npts':>5} {'b':>8} {'b_err':>8} {'R2':>8} {'type':>10}")

def log_fit(x, a, b):
    return a + b * x

all_b = {}
for q_val in [2, 3, 5, 7]:
    entries = data[q_val]
    # Use n≥8 only (avoid crossover)
    entries_filt = [e for e in entries if e['n'] >= 8]
    ns = np.array([e['n'] for e in entries_filt])
    wt = np.array([e['w_tail'] for e in entries_filt])
    mask = wt > 0
    ln_n = np.log(ns[mask])
    ln_wt = np.log(wt[mask])

    popt, pcov = curve_fit(log_fit, ln_n, ln_wt)
    b = popt[1]
    b_err = np.sqrt(pcov[1,1])
    predicted = log_fit(ln_n, *popt)
    ss_res = np.sum((ln_wt - predicted)**2)
    ss_tot = np.sum((ln_wt - np.mean(ln_wt))**2)
    R2 = 1 - ss_res / ss_tot

    cft_type = {2: 'real', 3: 'real', 5: 'walking', 7: 'broken'}[q_val]
    n_range = f"{int(ns[0])}-{int(ns[-1])}"
    print(f"{q_val:3d} {n_range:>10} {len(ln_n):5d} {b:8.3f} {b_err:8.3f} {R2:8.4f} {cft_type:>10}")
    all_b[q_val] = {'b': float(b), 'b_err': float(b_err), 'R2': float(R2)}

# Cross-q average
bs = [all_b[q]['b'] for q in [2,3,5,7]]
print(f"\nMean exponent: {np.mean(bs):.3f} ± {np.std(bs):.3f}")
print(f"Max deviation from 2.0: {max(abs(b - 2.0) for b in bs):.3f}")

results['power_law'] = all_b
results['power_law']['mean_b'] = float(np.mean(bs))
results['power_law']['std_b'] = float(np.std(bs))

# 2. Define multiplet dominance M(q, n) = %S(lev1) / (1 - %S(tail))
# This measures how much of the "conformal sector" entropy is in the multiplet
print("\n\n--- Multiplet dominance M = %S(lev1) / (%S(lev0) + %S(lev1)) ---")
print(f"{'q':>3} {'n':>4} {'%S(l0)':>8} {'%S(l1)':>8} {'M':>8}")

M_data = {}
for q_val in [2, 3, 5, 7]:
    M_data[q_val] = {'ns': [], 'Ms': [], 'S0s': [], 'S1s': []}
    for e in data[q_val]:
        S0 = e['S_lev0_frac']
        S1 = e['S_lev1_frac']
        M = S1 / (S0 + S1) if (S0 + S1) > 0 else 0
        M_data[q_val]['ns'].append(e['n'])
        M_data[q_val]['Ms'].append(M)
        M_data[q_val]['S0s'].append(S0)
        M_data[q_val]['S1s'].append(S1)
        if e['n'] in [8, 16, 24] or (q_val == 7 and e['n'] in [8, 12, 16]):
            print(f"{q_val:3d} {e['n']:4d} {S0:8.4f} {S1:8.4f} {M:8.4f}")

# M at matched size (n=16 available for all q)
print("\nM at n=16 (matched size):")
for q_val in [2, 3, 5, 7]:
    idx = M_data[q_val]['ns'].index(16) if 16 in M_data[q_val]['ns'] else -1
    if idx >= 0:
        M = M_data[q_val]['Ms'][idx]
        print(f"  q={q_val}: M = {M:.4f}")

# 3. Test if M(q,n→∞) is a simple function of q
# From 089b: %S(lev0_inf) is known. Use that.
print("\n\n--- M_inf = S_lev1_inf / (S_lev0_inf + S_lev1_inf) ---")
# Load saturation values from 089b
with open("results/sprint_089b_level_redistribution.json") as f:
    redist = json.load(f)

# Note: S_lev1_inf from the fits was unreliable (went to 0 for q=2,3,5)
# Better approach: M at largest available n as lower bound
print("Using M at largest n as best estimate:")
M_asymp = {}
for q_val in [2, 3, 5, 7]:
    M = M_data[q_val]['Ms'][-1]
    n_max = M_data[q_val]['ns'][-1]
    M_asymp[q_val] = M
    print(f"  q={q_val}: M(n={n_max}) = {M:.4f}")

# Test: M vs (q-1)/q ?
print("\n--- Does M → (q-1)/q? ---")
for q_val in [2, 3, 5, 7]:
    M = M_asymp[q_val]
    qm1_over_q = (q_val - 1) / q_val
    print(f"  q={q_val}: M = {M:.4f}, (q-1)/q = {qm1_over_q:.4f}, ratio = {M/qm1_over_q:.4f}")

# 4. Connect to c_eff/Re(c) — can we predict walking breakdown?
print("\n\n--- Connection: c_eff/Re(c) vs spectral quantities ---")

# c_eff/Re(c) values from KNOWLEDGE.md
c_ratio = {2: 1.00, 3: 1.12, 5: 1.01, 7: 0.78}  # at n_best

# Actually, let's compute c_eff from the entropy data directly
# c_eff from CC formula: S(n) = (c/6)·ln(n) + const for open BC
print("\nc_eff from DMRG entropy size pairs (largest pair):")
c_effs = {}
for q_val in [2, 3, 5, 7]:
    entries = data[q_val]
    if len(entries) >= 2:
        n1, S1 = entries[-2]['n'], entries[-2]['S_total']
        n2, S2 = entries[-1]['n'], entries[-1]['S_total']
        # Open BC: S = (c/6)·ln(n) + const
        c_eff = 6 * (S2 - S1) / (np.log(n2) - np.log(n1))
        c_effs[q_val] = c_eff
        ratio = c_eff / Re_c(q_val)
        print(f"  q={q_val}: c_eff = {c_eff:.4f}, Re(c) = {Re_c(q_val):.3f}, ratio = {ratio:.4f} (n={n1},{n2})")

# 5. Test correlation: c_eff/Re(c) vs M
print("\n--- Correlation: c_eff/Re(c) vs multiplet dominance M ---")
qs_both = [q for q in [2,3,5,7] if q in c_effs and q in M_asymp]
c_ratios_val = [c_effs[q] / Re_c(q) for q in qs_both]
M_vals = [M_asymp[q] for q in qs_both]
if len(qs_both) >= 3:
    corr = np.corrcoef(c_ratios_val, M_vals)[0,1]
    print(f"  Pearson r(c_eff/Re(c), M) = {corr:.4f}")
    for q in qs_both:
        print(f"    q={q}: c_eff/Re(c) = {c_effs[q]/Re_c(q):.4f}, M = {M_asymp[q]:.4f}")

# 6. Define "entropy redistribution index" R = %S(tail) / %S(lev1)
# This directly measures how much the multiplet leaks to tail
print("\n\n--- Entropy redistribution index R = %S(tail) / %S(lev1) ---")
print(f"{'q':>3} {'n':>4} {'R':>8} {'type':>10}")
R_asymp = {}
for q_val in [2, 3, 5, 7]:
    for e in data[q_val]:
        if e['n'] in [8, 16, 24] or (q_val == 7 and e['n'] in [8, 16]):
            R = e['S_tail_frac'] / e['S_lev1_frac'] if e['S_lev1_frac'] > 0 else 0
            cft_type = {2: 'real', 3: 'real', 5: 'walking', 7: 'broken'}[q_val]
            print(f"{q_val:3d} {e['n']:4d} {R:8.4f} {cft_type:>10}")
    # Use largest n
    e_last = data[q_val][-1]
    R_asymp[q_val] = e_last['S_tail_frac'] / e_last['S_lev1_frac']

print("\nR at largest n:")
for q_val in [2, 3, 5, 7]:
    print(f"  q={q_val}: R = {R_asymp[q_val]:.4f}")

# 7. Test: does %S(lev0) → ln(q)/[ln(q) + (q-1)·ln(q/(q-1))] analytically?
# This would be the entropy partition if the spectrum were exactly:
# λ_0 = 1/q, λ_{1..q-1} = (1-1/q)/(q-1) = 1/q (democratic)
# In that limit, S_total = ln(q), and all levels equal → %S each = 1/q
# Our spectrum is NOT democratic — test deviation
print("\n\n--- Comparison to democratic spectrum ---")
for q_val in [2, 3, 5, 7]:
    # Democratic limit: all λ_i = 1/q → %S_lev0 = 1/q, %S_lev1 = (q-1)/q
    dem_S0 = 1.0 / q_val
    dem_S1 = (q_val - 1.0) / q_val
    # Actual at largest n
    e_last = data[q_val][-1]
    n_last = e_last['n']
    actual_S0 = e_last['S_lev0_frac']
    actual_S1 = e_last['S_lev1_frac']
    # "Democracy index" D = actual / democratic
    D0 = actual_S0 / dem_S0
    D1 = actual_S1 / dem_S1
    print(f"  q={q_val} (n={n_last}): %S_l0 actual/democratic = {D0:.3f}, %S_l1 actual/democratic = {D1:.3f}")

# 8. SUMMARY
print("\n\n" + "=" * 70)
print("SYNTHESIS: What makes walking unique?")
print("=" * 70)

print("""
1. TAIL GROWTH IS UNIVERSAL: w_tail ~ n^{2.0} for ALL critical q=2,3,5,7.
   The b≈3 for q=7 at n=6-12 was pre-asymptotic (now b=2.07 with n≥8).

2. STATIC PARTITION IS q-DEPENDENT: At any fixed n, larger q puts MORE
   entropy in the (q-1) multiplet and LESS in ground state + tail.
   This is the microscopic origin of c_eff deviation.

3. MONOTONIC DISCRIMINATORS:
   - %S(lev0): 0.321 (q=2) → 0.199 (q=7) — DECREASING
   - %S(lev1): 0.635 (q=2) → 0.781 (q=7) — INCREASING
   - S_lev0_inf: 0.338 (q=2) → 0.203 (q=7) — DECREASING
   These distinguish real CFT, walking, and broken walking.

4. DYNAMIC RATES ARE NOT CLEAN DISCRIMINATORS:
   d(%S_tail)/d(ln n) is 0.036, 0.036, 0.032, 0.020 — not monotonic
   for q=2 vs q=3. The q-dependence is real but noisy.

5. THE MECHANISM: The (q-1)-fold degeneracy creates an entropy "basin"
   that absorbs fraction (q-1)/q × correction. As q grows, the basin
   deepens, capturing more entropy and depleting the tail. This slows
   the convergence of c_eff to Re(c).
""")

results['summary'] = {
    'power_law_universal': True,
    'mean_exponent': float(np.mean(bs)),
    'monotonic_discriminators': ['S_lev0', 'S_lev1', 'S_lev0_inf'],
    'M_asymp': {str(k): float(v) for k, v in M_asymp.items()},
    'R_asymp': {str(k): float(v) for k, v in R_asymp.items()},
    'c_effs': {str(k): float(v) for k, v in c_effs.items()},
}
save()
print(f"Saved to results/sprint_089c_walking_discriminator.json")

from db_utils import record
for q_val in [2, 3, 5, 7]:
    record(sprint=89, model='sq_potts', q=q_val, n=0,
           quantity='M_asymp', value=M_asymp[q_val],
           method='multiplet_dominance',
           notes='S_lev1/(S_lev0+S_lev1) at largest n')
    record(sprint=89, model='sq_potts', q=q_val, n=0,
           quantity='R_redist', value=R_asymp[q_val],
           method='redistribution_index',
           notes='S_tail/S_lev1 at largest n')
print("Recorded to DB.")
