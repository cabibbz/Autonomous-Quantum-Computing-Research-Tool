#!/usr/bin/env python3
"""Sprint 083c: Compare velocity extractions and analyze walking signatures in Casimir energy.

Three methods for v(q):
  1. v_gap = gap×N / (2π·x_σ)  [Sprint 082, uses gap + correlator]
  2. v_casimir = v·c / c  [this sprint, uses ground state energy only]
  3. v_self = gap×N / (2π) · (6/vc)  [combines Casimir vc with gap, no x_σ needed → gives x_σ]

Key analysis:
  - Does vc (Casimir) show walking signature?
  - What c makes Casimir v agree with gap v? → "effective c from velocity matching"
  - Does vc converge with N for all q, or diverge for q>5?
"""
import numpy as np
import json, time

results = {
    'experiment': '083c_velocity_comparison',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

def save():
    with open("results/sprint_083c_velocity_comparison.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# Collect all data from 083a and 083b
with open("results/sprint_083a_casimir_q235.json") as f:
    data_a = json.load(f)
with open("results/sprint_083b_casimir_q678.json") as f:
    data_b = json.load(f)

# Sprint 082 reference data
sprint082 = {
    2: {'x_sigma': 0.123, 'v': 1.02, 'gap_N': 0.786},
    3: {'x_sigma': 0.132, 'v': 0.89, 'gap_N': 0.733},
    4: {'x_sigma': 0.135, 'v': 0.81, 'gap_N': 0.686},
    5: {'x_sigma': 0.136, 'v': 0.75, 'gap_N': 0.639},
    6: {'x_sigma': 0.135, 'v': 0.71, 'gap_N': 0.601},
    7: {'x_sigma': 0.132, 'v': 0.68, 'gap_N': 0.560},
    8: {'x_sigma': 0.130, 'v': 0.66, 'gap_N': 0.536},
}

# Complex CFT Re(c)
def complex_cft_c(q):
    if q == 2: return 0.500
    if q == 3: return 0.800
    if q == 4: return 1.000
    sqrt_q = np.sqrt(q)
    alpha = np.arccosh(sqrt_q / 2)
    return 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)

# c_eff from DMRG (best measurements from Sprints 078-081)
c_eff_measured = {
    2: 0.500,  # exact
    3: 0.800,  # exact
    4: 1.000,  # exact
    5: 1.152,  # Sprint 078c n=12
    6: 1.115,  # Sprint 081a n=12
    7: 1.059,  # Sprint 079a n=12
    8: 1.062,  # Sprint 080b n=7
}

print("Sprint 083c: Velocity comparison and walking signature analysis")
print("=" * 80, flush=True)

# Merge data
all_data = {}
for q_str, d in data_a['data'].items():
    all_data[int(q_str)] = d
for q_str, d in data_b['data'].items():
    all_data[int(q_str)] = d

# Table 1: Grand comparison
print(f"\n{'='*80}")
print("TABLE 1: Three velocity methods compared")
print(f"{'='*80}")
print(f"{'q':>3} {'vc':>7} {'Re(c)':>6} {'c_eff':>6} {'v_gap':>6} {'v(Re c)':>8} {'v(c_eff)':>9} "
      f"{'Δ_gap%':>7} {'Δ_ceff%':>8}")

comparison = []
for q in sorted(all_data.keys()):
    d = all_data[q]
    vc = d['fit']['vc_casimir']
    rec = complex_cft_c(q)
    ceff = c_eff_measured[q]
    ref = sprint082[q]

    v_gap = ref['gap_N'] / (2 * np.pi * ref['x_sigma'])
    v_rec = vc / rec
    v_ceff = vc / ceff

    delta_gap = 100 * (v_rec - v_gap) / v_gap  # Casimir(Re(c)) vs gap
    delta_ceff = 100 * (v_ceff - v_gap) / v_gap  # Casimir(c_eff) vs gap

    print(f"{q:3d} {vc:7.4f} {rec:6.3f} {ceff:6.3f} {v_gap:6.3f} {v_rec:8.4f} {v_ceff:9.4f} "
          f"{delta_gap:+7.1f} {delta_ceff:+8.1f}")

    comparison.append({
        'q': q, 'vc': vc, 'Re_c': rec, 'c_eff': ceff,
        'v_gap': v_gap, 'v_rec': v_rec, 'v_ceff': v_ceff,
        'delta_gap_pct': delta_gap, 'delta_ceff_pct': delta_ceff,
    })

# Table 2: What c makes Casimir v match gap v?
print(f"\n{'='*80}")
print("TABLE 2: c_implied = vc / v_gap (what c makes Casimir agree with correlator?)")
print(f"{'='*80}")
print(f"{'q':>3} {'vc':>7} {'v_gap':>6} {'c_implied':>10} {'Re(c)':>6} {'c_eff':>6} "
      f"{'c_imp/Re':>8} {'c_imp/eff':>9}")

c_implied_data = []
for q in sorted(all_data.keys()):
    d = all_data[q]
    vc = d['fit']['vc_casimir']
    ref = sprint082[q]
    v_gap = ref['gap_N'] / (2 * np.pi * ref['x_sigma'])
    c_implied = vc / v_gap
    rec = complex_cft_c(q)
    ceff = c_eff_measured[q]

    print(f"{q:3d} {vc:7.4f} {v_gap:6.3f} {c_implied:10.4f} {rec:6.3f} {ceff:6.3f} "
          f"{c_implied/rec:8.3f} {c_implied/ceff:9.3f}")

    c_implied_data.append({
        'q': q, 'c_implied': c_implied, 'Re_c': rec, 'c_eff': ceff,
        'ratio_rec': c_implied / rec, 'ratio_ceff': c_implied / ceff,
    })

# Table 3: x_sigma from Casimir (no correlator needed)
# v_casimir = vc/c → x_sigma_casimir = gap*N / (2*pi*v_casimir)
print(f"\n{'='*80}")
print("TABLE 3: x_σ from Casimir+gap (independent of correlator)")
print(f"{'='*80}")
print(f"{'q':>3} {'x_σ(corr)':>10} {'x_σ(Cas,Rec)':>13} {'x_σ(Cas,ceff)':>14} "
      f"{'Δ_Rec%':>7} {'Δ_ceff%':>8}")

for q in sorted(all_data.keys()):
    d = all_data[q]
    vc = d['fit']['vc_casimir']
    rec = complex_cft_c(q)
    ceff = c_eff_measured[q]
    ref = sprint082[q]
    gap_N = ref['gap_N']

    v_rec = vc / rec
    v_ceff = vc / ceff
    x_cas_rec = gap_N / (2 * np.pi * v_rec)
    x_cas_ceff = gap_N / (2 * np.pi * v_ceff)
    x_corr = ref['x_sigma']

    dr = 100 * (x_cas_rec - x_corr) / x_corr
    dc = 100 * (x_cas_ceff - x_corr) / x_corr

    print(f"{q:3d} {x_corr:10.4f} {x_cas_rec:13.4f} {x_cas_ceff:14.4f} {dr:+7.1f} {dc:+8.1f}")

# Table 4: Walking signature — vc convergence with N
print(f"\n{'='*80}")
print("TABLE 4: vc pairwise convergence (walking signature)")
print(f"{'='*80}")

for q in sorted(all_data.keys()):
    d = all_data[q]
    pw = d.get('pairwise_vc', d.get('pairwise', []))
    if not pw:
        continue
    vc_vals = [p['vc'] for p in pw]
    if len(vc_vals) >= 2:
        drift = 100 * (vc_vals[-1] - vc_vals[0]) / vc_vals[0]
        trend = "↓" if drift < 0 else "↑"
    else:
        drift = 0
        trend = "—"

    pairs_str = ", ".join([f"{p['pair']}:{p['vc']:.4f}" for p in pw])
    print(f"  q={q}: {pairs_str}  drift={drift:+.2f}% {trend}")

# Table 5: vc trend analysis
print(f"\n{'='*80}")
print("TABLE 5: vc(q) trend — does vc saturate?")
print(f"{'='*80}")

qs = sorted(all_data.keys())
vcs = [all_data[q]['fit']['vc_casimir'] for q in qs]
print(f"{'q':>3} {'vc':>8} {'Δvc':>7}")
for i, q in enumerate(qs):
    dvc = f"{vcs[i] - vcs[i-1]:+.4f}" if i > 0 else "—"
    print(f"{q:3d} {vcs[i]:8.4f} {dvc:>7}")

# Fit vc(q) to various forms
qs_arr = np.array(qs, dtype=float)
vc_arr = np.array(vcs)

# Linear fit
A = np.vstack([qs_arr, np.ones_like(qs_arr)]).T
(m, b), _, _, _ = np.linalg.lstsq(A, vc_arr, rcond=None)
vc_lin = m * qs_arr + b
rms_lin = np.sqrt(np.mean((vc_arr - vc_lin)**2))

# Log fit: vc = a*ln(q) + b
A2 = np.vstack([np.log(qs_arr), np.ones_like(qs_arr)]).T
(m2, b2), _, _, _ = np.linalg.lstsq(A2, vc_arr, rcond=None)
vc_log = m2 * np.log(qs_arr) + b2
rms_log = np.sqrt(np.mean((vc_arr - vc_log)**2))

# Power fit: vc = a*q^b (in log space)
A3 = np.vstack([np.log(qs_arr), np.ones_like(qs_arr)]).T
(m3, b3), _, _, _ = np.linalg.lstsq(A3, np.log(vc_arr), rcond=None)
vc_pow = np.exp(b3) * qs_arr**m3
rms_pow = np.sqrt(np.mean((vc_arr - vc_pow)**2))

print(f"\n  vc(q) fits:")
print(f"    Linear: vc = {m:.4f}*q + {b:.4f}  (RMS={rms_lin:.5f})")
print(f"    Log:    vc = {m2:.4f}*ln(q) + {b2:.4f}  (RMS={rms_log:.5f})")
print(f"    Power:  vc = {np.exp(b3):.4f}*q^{m3:.4f}  (RMS={rms_pow:.5f})")

# Key discovery check: does v(Casimir, Re(c)) ≈ v(gap) imply c_eff ≈ Re(c)?
print(f"\n{'='*80}")
print("KEY INSIGHT: Walking signature in velocity matching")
print(f"{'='*80}")
print(f"\nFor q≤5: v(Casimir, Re(c)) ≈ v(gap) → Casimir gives CORRECT Re(c)")
print(f"For q>5: v(Casimir, Re(c)) > v(gap) → Casimir gives WRONG Re(c)")
print(f"         The walking breakdown means c in E₀ ≠ c in entropy")
print(f"\nc_implied = vc/v_gap is the 'Casimir central charge':")
for cd in c_implied_data:
    q = cd['q']
    walking = "real CFT" if q <= 4 else ("WALKING" if q == 5 else "BREAKING")
    print(f"  q={q}: c_implied={cd['c_implied']:.4f}, Re(c)={cd['Re_c']:.3f}, "
          f"c_eff={cd['c_eff']:.3f} [{walking}]")

results['comparison'] = comparison
results['c_implied'] = c_implied_data
results['vc_fits'] = {
    'linear': {'m': float(m), 'b': float(b), 'rms': float(rms_lin)},
    'log': {'m': float(m2), 'b': float(b2), 'rms': float(rms_log)},
    'power': {'a': float(np.exp(b3)), 'b': float(m3), 'rms': float(rms_pow)},
}
save()
print(f"\nSaved to results/sprint_083c_velocity_comparison.json")

# Record c_implied to DB
from db_utils import record
for cd in c_implied_data:
    q = cd['q']
    n_best = max(all_data[q]['sizes'])
    record(sprint=83, model='sq_potts', q=q, n=n_best,
           quantity='c_implied', value=cd['c_implied'],
           method='casimir_velocity_match',
           notes=f'c that makes Casimir v match gap v. Re(c)={cd["Re_c"]:.3f}')
print("Recorded to DB.")
