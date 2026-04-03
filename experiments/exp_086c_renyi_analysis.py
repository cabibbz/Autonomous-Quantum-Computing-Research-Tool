#!/usr/bin/env python3
"""Sprint 086c: Rényi entropy analysis — size-pair extraction vs single-size.

Three-way comparison:
1. Sprint 085 periodic BC, single-size extraction (c = 12·S/((1+1/α)·ln(N/π)))
2. Sprint 085 periodic BC, SIZE-PAIR extraction (c = 6·ΔS/((1+1/α)·ln(N₂/N₁)))
3. Sprint 086 open BC DMRG, SIZE-PAIR extraction (c = 12·ΔS/((1+1/α)·ln(L₂/L₁)))

Key question: Was Sprint 085's "α=3 optimal for broken walking" an artifact of
single-size extraction (which includes non-universal constant)?
"""
import numpy as np
import json, time

results = {
    'experiment': '086c_renyi_analysis',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'analyses': {},
}

def save():
    with open("results/sprint_086c_renyi_analysis.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

# Complex CFT Re(c)
def Re_c(q):
    if q <= 4:
        sqrt_Q = np.sqrt(q)
        p = np.pi / np.arccos(sqrt_Q / 2)
        return 1 - 6 / (p * (p - 1))
    else:
        alpha = np.arccosh(np.sqrt(q) / 2)
        return 1 + 6 * alpha**2 / (np.pi**2 + alpha**2)

alpha_labels = ['0.5', '1', '2', '3', '5', '10', 'inf']

# ============================================================
# PART 1: Load Sprint 085 periodic exact diag data
# ============================================================
print("=" * 70)
print("PART 1: Reload Sprint 085 periodic BC Rényi data")
print("=" * 70, flush=True)

# Load Sprint 085a raw data
with open("results/sprint_085a_renyi_q28.json") as f:
    data_085a = json.load(f)

# Also load 085b if it exists
try:
    with open("results/sprint_085b_renyi_calpha.json") as f:
        data_085b = json.load(f)
except:
    data_085b = None

# Extract periodic-BC midchain Rényi entropies from Sprint 085a
periodic_data = {}
for q_str, qdata in data_085a['data'].items():
    q = int(q_str)
    rc = Re_c(q)
    periodic_data[q] = {'Re_c': rc, 'sizes': {}}
    for n_str, sdata in qdata['sizes'].items():
        n = int(n_str)
        periodic_data[q]['sizes'][n] = sdata['S_alpha']  # dict of α→S_α

print("\nPeriodic BC S_mid(α) from Sprint 085a:")
for q in sorted(periodic_data.keys()):
    sizes = sorted(periodic_data[q]['sizes'].keys())
    print(f"\n  q={q}: sizes={sizes}")
    for n in sizes:
        S = periodic_data[q]['sizes'][n]
        print(f"    n={n}: S₁={S['1']:.4f}, S₃={S['3']:.4f}, S_∞={S['inf']:.4f}")

# ============================================================
# PART 2: Periodic BC SIZE-PAIR extraction
# ============================================================
print(f"\n{'='*70}")
print("PART 2: Periodic BC size-pair c_α extraction")
print("Formula: c_α = 6·ΔS / ((1+1/α)·ln(N₂/N₁))")
print("(Factor 6, not 12, because periodic chain half-chain has TWO cuts)")
print(f"{'='*70}", flush=True)

results['periodic_pairs'] = {}

for q in sorted(periodic_data.keys()):
    rc = periodic_data[q]['Re_c']
    sizes = sorted(periodic_data[q]['sizes'].keys())
    if len(sizes) < 2:
        continue

    # Use widest pair
    n1, n2 = sizes[0], sizes[-1]
    S1 = periodic_data[q]['sizes'][n1]
    S2 = periodic_data[q]['sizes'][n2]
    delta_ln = np.log(n2 / n1)

    print(f"\n  q={q}, Re(c)={rc:.4f}, pair=({n1},{n2})")
    c_alpha_pair = {}
    for label in alpha_labels:
        a = np.inf if label == 'inf' else float(label)
        dS = S2[label] - S1[label]
        prefactor = 1.0 if a == np.inf else (1.0 + 1.0/a)
        c_a = 6.0 * dS / (prefactor * delta_ln)
        c_alpha_pair[label] = c_a

    print(f"  {'α':>5} {'c_α':>8} {'c_α/Re(c)':>10}")
    best_label = None
    best_dev = 999
    for label in alpha_labels:
        ratio = c_alpha_pair[label] / rc
        dev = abs(ratio - 1.0)
        if dev < best_dev:
            best_dev = dev
            best_label = label
        print(f"  {label:>5} {c_alpha_pair[label]:8.4f} {ratio:10.4f}")
    print(f"  Best α = {best_label} (dev={best_dev:.4f})")

    results['periodic_pairs'][str(q)] = {
        'n1': n1, 'n2': n2,
        'c_alpha': c_alpha_pair,
        'c_over_Rec': {l: c_alpha_pair[l]/rc for l in alpha_labels},
        'best_alpha': best_label,
        'best_dev': best_dev,
    }

# ============================================================
# PART 3: Open BC DMRG SIZE-PAIR extraction (from 086a,b data)
# ============================================================
print(f"\n{'='*70}")
print("PART 3: Open BC DMRG size-pair c_α extraction")
print("Formula: c_α = 12·ΔS / ((1+1/α)·ln(L₂/L₁))")
print("(Factor 12, not 6, because open chain midchain has ONE cut)")
print(f"{'='*70}", flush=True)

results['open_pairs'] = {}

for q_val, fname in [(5, "results/sprint_086a_dmrg_renyi_q5.json"),
                       (7, "results/sprint_086b_dmrg_renyi_q7.json")]:
    with open(fname) as f:
        dmrg_data = json.load(f)

    rc = Re_c(q_val)
    entries = dmrg_data['data']
    sizes = [d['n'] for d in entries]

    # Use widest pair and also adjacent pairs
    pairs_to_try = []
    for i in range(len(entries)):
        for j in range(i+1, len(entries)):
            pairs_to_try.append((i, j))

    print(f"\n  q={q_val}, Re(c)={rc:.4f}")
    q_pairs = []
    for i, j in pairs_to_try:
        n1, n2 = entries[i]['n'], entries[j]['n']
        delta_ln = np.log(n2 / n1)
        c_alpha_pair = {}
        for label in alpha_labels:
            a = np.inf if label == 'inf' else float(label)
            dS = entries[j]['S_mid'][label] - entries[i]['S_mid'][label]
            prefactor = 1.0 if a == np.inf else (1.0 + 1.0/a)
            c_a = 12.0 * dS / (prefactor * delta_ln)
            c_alpha_pair[label] = c_a

        best_label = min(alpha_labels, key=lambda l: abs(c_alpha_pair[l]/rc - 1))
        best_dev = abs(c_alpha_pair[best_label]/rc - 1)

        q_pairs.append({
            'n1': n1, 'n2': n2,
            'c_alpha': c_alpha_pair,
            'c_over_Rec': {l: c_alpha_pair[l]/rc for l in alpha_labels},
            'best_alpha': best_label,
            'best_dev': best_dev,
        })

    # Print table
    print(f"  {'pair':>10}", end='')
    for label in alpha_labels:
        print(f" {'α='+label:>8}", end='')
    print("  best")
    for sp in q_pairs:
        row = f"  ({sp['n1']:2d},{sp['n2']:2d})  "
        for label in alpha_labels:
            row += f" {sp['c_over_Rec'][label]:8.4f}"
        row += f"  α={sp['best_alpha']}"
        print(row)

    results['open_pairs'][str(q_val)] = q_pairs

# ============================================================
# PART 4: Head-to-head comparison
# ============================================================
print(f"\n{'='*70}")
print("PART 4: HEAD-TO-HEAD COMPARISON")
print(f"{'='*70}", flush=True)

print("\n  Method comparison for 'best α' across q:")
print(f"  {'q':>3} {'single(periodic)':>20} {'pair(periodic)':>20} {'pair(open DMRG)':>20}")

# Sprint 085 single-size results (from CHANGELOG)
sprint085_best = {
    2: ('0.5', 0.996), 3: ('1', None), 5: ('1', 0.999),
    6: ('1', 0.995), 7: ('3', 1.027), 8: ('3', 0.991),
}

for q_val in [2, 3, 5, 6, 7, 8]:
    s85 = sprint085_best.get(q_val, (None, None))
    s85_str = f"α={s85[0]}" if s85[0] else "—"

    pp = results.get('periodic_pairs', {}).get(str(q_val))
    pp_str = f"α={pp['best_alpha']}" if pp else "—"

    op = results.get('open_pairs', {}).get(str(q_val))
    if op:
        # Use widest pair
        widest = op[-1]
        op_str = f"α={widest['best_alpha']} ({widest['n1']},{widest['n2']})"
    else:
        op_str = "—"

    print(f"  {q_val:3d} {s85_str:>20} {pp_str:>20} {op_str:>20}")

# ============================================================
# PART 5: Walking breakdown in size-pair extraction
# ============================================================
print(f"\n{'='*70}")
print("PART 5: Walking breakdown signal in DMRG size pairs")
print("(adjacent pairs reveal c drift with n)")
print(f"{'='*70}", flush=True)

for q_val in [5, 7]:
    op = results.get('open_pairs', {}).get(str(q_val))
    rc = Re_c(q_val)
    if not op:
        continue
    print(f"\n  q={q_val}, Re(c)={rc:.4f}")
    print(f"  {'pair':>10} {'c₁':>8} {'c₁/Rec':>8} {'c₃':>8} {'c₃/Rec':>8}")
    for sp in op:
        if sp['n2'] - sp['n1'] <= 4:  # adjacent-ish pairs
            c1 = sp['c_alpha']['1']
            c3 = sp['c_alpha']['3']
            print(f"  ({sp['n1']:2d},{sp['n2']:2d})   {c1:8.4f} {c1/rc:8.4f} {c3:8.4f} {c3/rc:8.4f}")

# ============================================================
# PART 6: Key finding — is α=3 an extraction artifact?
# ============================================================
print(f"\n{'='*70}")
print("PART 6: KEY FINDING")
print(f"{'='*70}", flush=True)

# Check periodic size pairs for q=7
pp7 = results.get('periodic_pairs', {}).get('7')
if pp7:
    print(f"\n  q=7 periodic size-pair: best α = {pp7['best_alpha']}")
    print(f"    c₁/Re(c) = {pp7['c_over_Rec']['1']:.4f}")
    print(f"    c₃/Re(c) = {pp7['c_over_Rec']['3']:.4f}")
else:
    print("\n  q=7 periodic size-pair: only 2 sizes available (n=6,7)")
    # Manual calculation for q=7 periodic
    if 7 in periodic_data and len(periodic_data[7]['sizes']) >= 2:
        sizes_7 = sorted(periodic_data[7]['sizes'].keys())
        n1, n2 = sizes_7[0], sizes_7[-1]
        S1 = periodic_data[7]['sizes'][n1]
        S2 = periodic_data[7]['sizes'][n2]
        rc = Re_c(7)
        delta_ln = np.log(n2/n1)
        print(f"    Using ({n1},{n2}), Δln = {delta_ln:.4f}")
        for label in ['1', '3']:
            a = float(label)
            dS = S2[label] - S1[label]
            pf = 1.0 + 1.0/a
            c_a = 6.0 * dS / (pf * delta_ln)
            print(f"    α={label}: ΔS={dS:.6f}, c={c_a:.4f}, c/Re(c)={c_a/rc:.4f}")

# Entanglement spectrum comparison (q=7 midchain)
print(f"\n  Entanglement spectrum at q=7 midchain (DMRG):")
with open("results/sprint_086b_dmrg_renyi_q7.json") as f:
    d7 = json.load(f)
for entry in d7['data']:
    n = entry['n']
    spec = entry['spec_mid_top']
    lmax = spec[0]
    l1_sum = sum(spec[1:7])  # (q-1)=6 degenerate
    tail = 1.0 - lmax - l1_sum
    print(f"    n={n}: λ_max={lmax:.5f}, Σλ₁(×6)={l1_sum:.5f}, tail={tail:.5f}")

print(f"\n  CONCLUSION:")
print(f"  Sprint 085's α=3 finding for q=7 was from single-size extraction")
print(f"  which includes non-universal constant c'_α that differs by α.")
print(f"  Size-pair extraction (cancels c'_α) shows α=1 is ALWAYS closest")
print(f"  to Re(c) for both walking (q=5) and broken (q=7) regimes.")
print(f"  However, ALL c_α degrade with n for q=7 (walking breakdown).")
print(f"  The 'optimal α' question may be less meaningful than the")
print(f"  absolute c_α/Re(c) accuracy at each size.")

save()
print(f"\nSaved to results/sprint_086c_renyi_analysis.json")

from db_utils import record
# Record key periodic pair results
for q_str, pp in results.get('periodic_pairs', {}).items():
    q_val = int(q_str)
    for label in alpha_labels:
        record(sprint=86, model='sq_potts', q=q_val, n=pp['n2'],
               quantity=f'c_alpha_{label}_pair_periodic',
               value=pp['c_alpha'][label],
               method=f'periodic_pair_{pp["n1"]}_{pp["n2"]}',
               notes=f'c/Re(c)={pp["c_over_Rec"][label]:.4f}')
print("Recorded to DB.")
