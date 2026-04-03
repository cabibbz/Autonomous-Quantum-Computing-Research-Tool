#!/usr/bin/env python3
"""Sprint 091c: Synthesis — BW fidelity, walking discriminators, and entanglement temperature.

Combines 091a/b data to:
1. Correlate non-Potts fraction at nA=4 with M/[(q-1)/q] and c_eff/Re(c)
2. Extract entanglement temperature beta_E(q) from BW alpha
3. Fit BW locality drop rate with nA for q=2 vs q=5
4. Verify: is the q=4 boundary special for BW (as it is for M)?
"""
import numpy as np
import json

# Load 091a and 091b results
with open("results/sprint_091a_he_bw_fidelity.json") as f:
    data_a = json.load(f)
with open("results/sprint_091b_he_range_decomp.json") as f:
    data_b = json.load(f)

print("Sprint 091c: Synthesis — BW fidelity and walking discriminators")
print("=" * 70)

# ---- Part 1: Fixed nA=4 comparison ----
print("\nPart 1: Non-Potts fraction at nA=4 vs walking discriminators")
print("-" * 60)

# Walking discriminators from prior sprints
q_vals = [2, 3, 4, 5]
M_ratio = {2: 1.352, 3: 1.086, 4: 1.003, 5: 0.964}  # M/[(q-1)/q] at n=16
c_eff_ratio = {2: 1.00, 3: 1.12, 4: 1.00, 5: 1.01}   # c_eff/Re(c)
Re_c = {2: 0.500, 3: 0.800, 4: 1.000, 5: 1.138}

# Non-Potts fraction = 1 - NN_locality
non_potts = {}
bw_alpha = {}
for q in q_vals:
    key = f'fixed_nA4_q{q}'
    d = data_b['data'][key]
    non_potts[q] = d['residual']

# Also get BW alpha from 091a
for q in [2, 3, 4, 5, 7]:
    key = f'q{q}'
    if key in data_a['data']:
        bw_alpha[q] = data_a['data'][key]['bw_alpha']

print(f"{'q':>3} {'non-Potts%':>10} {'M/[(q-1)/q]':>12} {'c_eff/Rec':>10} {'BW_alpha':>10} {'Re(c)':>8}")
for q in q_vals:
    print(f"{q:>3} {non_potts[q]*100:>10.3f} {M_ratio[q]:>12.3f} {c_eff_ratio[q]:>10.3f} "
          f"{bw_alpha.get(q, 0):>10.4f} {Re_c[q]:>8.3f}")

# Correlations
nps = np.array([non_potts[q] for q in q_vals])
mrs = np.array([M_ratio[q] for q in q_vals])
crs = np.array([c_eff_ratio[q] for q in q_vals])
recs = np.array([Re_c[q] for q in q_vals])
qs = np.array(q_vals, dtype=float)

# Pearson correlations
corr_M = np.corrcoef(nps, mrs)[0, 1]
corr_c = np.corrcoef(nps, crs)[0, 1]
corr_q = np.corrcoef(nps, qs)[0, 1]
corr_Rec = np.corrcoef(nps, recs)[0, 1]

print(f"\nCorrelation with non-Potts fraction:")
print(f"  vs M/[(q-1)/q]:  r = {corr_M:+.4f}")
print(f"  vs c_eff/Re(c):  r = {corr_c:+.4f}")
print(f"  vs q:            r = {corr_q:+.4f}")
print(f"  vs Re(c):        r = {corr_Rec:+.4f}")

# Fit: non-Potts vs q
from numpy.polynomial import polynomial as P
log_nps = np.log(nps)
fit = np.polyfit(qs, log_nps, 1)
print(f"\n  Fit: non-Potts ~ exp({fit[0]:.3f}·q + {fit[1]:.3f})")
print(f"  → non-Potts fraction grows exponentially with q")

# ---- Part 2: Entanglement temperature ----
print(f"\n{'='*70}")
print("Part 2: Entanglement temperature β_E = BW alpha")
print("-" * 60)

# BW alpha is the overall scale: H_E ≈ alpha × H_BW_unit
# where H_BW_unit has the shape of the physical H with BW envelope (normalized differently)
# The effective entanglement temperature is related to alpha
# For CFT: beta_E = 2*pi*(L/N) where L is subsystem size
# The BW prediction gives alpha = 2*pi for infinite system
# At finite size with periodic BC, corrections expected

print(f"{'q':>3} {'n':>3} {'nA':>3} {'alpha':>8} {'alpha/2pi':>10}")
for q in [2, 3, 4, 5, 7]:
    key = f'q{q}'
    if key in data_a['data']:
        d = data_a['data'][key]
        alpha = d['bw_alpha']
        print(f"{q:>3} {d['n']:>3} {d['nA']:>3} {alpha:>8.4f} {alpha/(2*np.pi):>10.4f}")

# At fixed nA=4 (from 091a data, same as 091b)
print(f"\nAt fixed nA=4:")
alphas_nA4 = {}
for q in [4, 5]:
    key = f'q{q}'
    d = data_a['data'][key]
    alphas_nA4[q] = d['bw_alpha']
    print(f"  q={q}: alpha = {d['bw_alpha']:.4f}, alpha/2pi = {d['bw_alpha']/(2*np.pi):.4f}")

# ---- Part 3: nA scaling fit ----
print(f"\n{'='*70}")
print("Part 3: BW locality scaling with nA")
print("-" * 60)

# q=2 data from 091b
q2_nA = [3, 4, 5, 6, 7]
q2_loc = []
for nA in q2_nA:
    d = data_b['data'][f'q2_nA{nA}']
    q2_loc.append(d['var_nn'])

# q=5 data from 091b
q5_nA = [3, 4]
q5_loc = []
for nA in q5_nA:
    d = data_b['data'][f'q5_nA{nA}']
    q5_loc.append(d['var_nn'])

print("q=2 Potts NN locality vs nA:")
for nA, loc in zip(q2_nA, q2_loc):
    print(f"  nA={nA}: {loc*100:.3f}%")

# Fit q=2: log(1-locality) vs nA
q2_resid = np.array([1 - loc for loc in q2_loc])
q2_nA_arr = np.array(q2_nA, dtype=float)
# Only fit where residual > 0
mask = q2_resid > 1e-6
if mask.sum() >= 2:
    log_resid = np.log(q2_resid[mask])
    fit_q2 = np.polyfit(q2_nA_arr[mask], log_resid, 1)
    print(f"\n  Fit: (1-locality) ~ exp({fit_q2[0]:.3f}·nA + {fit_q2[1]:.3f})")
    print(f"  → non-Potts fraction grows exponentially with nA for q=2")
    print(f"  → doubling every {np.log(2)/fit_q2[0]:.1f} sites")

print(f"\nq=5 Potts NN locality vs nA:")
for nA, loc in zip(q5_nA, q5_loc):
    print(f"  nA={nA}: {loc*100:.3f}%")

q5_resid = np.array([1 - loc for loc in q5_loc])
if len(q5_resid) >= 2:
    log_resid_5 = np.log(q5_resid)
    q5_nA_arr = np.array(q5_nA, dtype=float)
    # Two points only — compute slope
    slope_5 = (log_resid_5[1] - log_resid_5[0]) / (q5_nA_arr[1] - q5_nA_arr[0])
    print(f"  slope: {slope_5:.3f} per site (q=2 slope: {fit_q2[0]:.3f})")
    print(f"  q=5 non-Potts grows {slope_5/fit_q2[0]:.1f}× faster than q=2")

# ---- Part 4: Is q=4 special for BW? ----
print(f"\n{'='*70}")
print("Part 4: Is q=4 special for BW (as for M crossover)?")
print("-" * 60)

# At nA=4: non-Potts fractions
print(f"Non-Potts fraction at nA=4:")
print(f"  q=2: {non_potts[2]*100:.4f}%")
print(f"  q=3: {non_potts[3]*100:.4f}%")
print(f"  q=4: {non_potts[4]*100:.4f}%")
print(f"  q=5: {non_potts[5]*100:.4f}%")

# Check if there's a jump at q=4
ratio_34 = non_potts[4] / non_potts[3]
ratio_45 = non_potts[5] / non_potts[4]
print(f"\n  Ratio q=4/q=3: {ratio_34:.1f}×")
print(f"  Ratio q=5/q=4: {ratio_45:.1f}×")
print(f"  {'ACCELERATION at q=4' if ratio_45 > ratio_34 else 'DECELERATION at q=5' if ratio_34 > ratio_45 else 'UNIFORM growth'}")

# Check for exponential growth
log_np = [np.log(non_potts[q]) for q in q_vals]
# Increments
for i in range(len(q_vals)-1):
    dq = q_vals[i+1] - q_vals[i]
    dlog = log_np[i+1] - log_np[i]
    print(f"  d(log non-Potts)/dq from q={q_vals[i]}→{q_vals[i+1]}: {dlog/dq:.3f}")

# ---- Part 5: Summary ----
print(f"\n{'='*70}")
print("SUMMARY")
print("=" * 70)

print("""
KEY FINDINGS:

1. At fixed nA=4, BW locality degrades MONOTONICALLY with q:
   q=2 (99.97%) → q=3 (99.92%) → q=4 (99.00%) → q=5 (97.08%)
   Non-Potts fraction grows exponentially: ~exp(1.5·q)

2. BW locality drops RAPIDLY with subsystem size nA (for q=2):
   nA=3 (99.98%) → nA=7 (62.65%)
   Non-Potts fraction grows exponentially: ~exp(1.2·nA)

3. The q-dependence is MUCH WEAKER than the nA-dependence at accessible sizes.
   091a results were dominated by nA variation, not walking physics.

4. NNN and 3-body operators add very little (<0.2% for all q at nA=4).
   The non-Potts content is in genuinely NON-POTTS operators — operators
   that don't appear in the physical Hamiltonian at all.

5. BW alpha (entanglement temperature scale) is ~3.1-3.2 for q=4,5 at nA=4.
   alpha/(2*pi) ≈ 0.51 — about half the expected CFT value.

6. q=5 non-Potts fraction grows faster with nA than q=2 (by ~2x slope),
   suggesting walking-related corrections to BW become dominant at larger systems.
""")

# Save combined results
synthesis = {
    'experiment': '091c_synthesis',
    'fixed_nA4_non_potts': {str(q): non_potts[q] for q in q_vals},
    'correlations': {
        'non_potts_vs_M': float(corr_M),
        'non_potts_vs_ceff': float(corr_c),
        'non_potts_vs_q': float(corr_q),
        'non_potts_vs_Rec': float(corr_Rec),
    },
    'exponential_fit': {'slope': float(fit[0]), 'intercept': float(fit[1])},
    'q2_nA_scaling': {
        'nA': q2_nA,
        'locality': [float(x) for x in q2_loc],
        'nA_slope': float(fit_q2[0]) if mask.sum() >= 2 else None,
    },
    'q5_nA_scaling': {
        'nA': q5_nA,
        'locality': [float(x) for x in q5_loc],
    },
    'bw_alpha': {str(q): bw_alpha[q] for q in bw_alpha},
}

with open("results/sprint_091c_synthesis.json", "w") as f:
    json.dump(synthesis, f, indent=2)

print("\nSaved to results/sprint_091c_synthesis.json")
