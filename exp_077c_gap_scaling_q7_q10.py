#!/usr/bin/env python3
"""Sprint 077c: S_q Potts Δ·N scaling at q=7 and q=10 — first-order test.

Key diagnostic: at a CONTINUOUS transition, Δ·N converges to a finite value.
At FIRST-ORDER, Δ·N → 0 as N → ∞ (gap closes exponentially).

Sprint 076c found S_q q=5 Δ·N drift of 4.8% (n=4→8) — small, CFT-like.
Sprint 077a found S_q q=7 g_c=0.144 with gap×N≈0.65 at n=4,6,8.

Here we systematically measure Δ·N at g_c for q=7 (n=4,6,8) and q=10 (n=4,5,6)
to look for first-order signatures at larger q.

Also compare to hybrid model at same q values.
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from gpu_utils import eigsh
import json, time

def build_hamiltonian(n, q, g, model='sq_potts'):
    """Build Hamiltonian with periodic BC."""
    dim = q ** n
    potts_2site = np.zeros(q ** 2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q ** 2, q ** 2), format='csr')

    if model == 'sq_potts':
        field_1site = np.ones((q, q)) - np.eye(q)
    else:
        X = np.zeros((q, q))
        for s in range(q):
            X[(s + 1) % q, s] = 1.0
        field_1site = X + X.T
    field_op = csr_matrix(field_1site)

    H = csr_matrix((dim, dim))
    for i in range(n - 1):
        left = q ** i
        right = q ** (n - i - 2)
        op = potts_op
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op

    diag_boundary = np.zeros(dim)
    for idx in range(dim):
        s0 = idx % q
        sn = (idx // (q ** (n - 1))) % q
        diag_boundary[idx] = 1.0 if s0 == sn else 0.0
    H = H - diags(diag_boundary, 0, shape=(dim, dim), format='csr')

    for i in range(n):
        left = q ** i
        right = q ** (n - i - 1)
        op = field_op.copy()
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op

    return H


def get_gap(n, q, g, model):
    H = build_hamiltonian(n, q, g, model)
    dim = H.shape[0]
    if dim < 500:
        evals = np.linalg.eigvalsh(H.toarray())
        return evals[1] - evals[0]
    k = min(6, dim - 1)
    evals, _ = eigsh(H, k=k, which='SA')
    evals = np.sort(evals)
    return evals[1] - evals[0]


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 077c: S_q Potts Δ·N scaling — first-order test", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '077c_gap_scaling',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
}

# Critical points (S_q Potts)
# q=5: g_c=0.200 (076a), q=7: g_c=0.144 (077a)
# q=10: estimate from field strength ratio
# Hybrid g_c(q=10)=0.684, field ratio (q-1)/2cos(2π/q)
# q=10: (9)/(2cos(π/5))=9/1.618=5.56 → S_q g_c ≈ 0.684/5.56 ≈ 0.123
# But this is rough — need to find it

configs = [
    {'q': 5, 'gc_sq': 0.200, 'gc_hy': 0.441, 'sizes_sq': [4, 6, 8], 'sizes_hy': [4, 6, 8]},
    {'q': 7, 'gc_sq': 0.144, 'gc_hy': 0.535, 'sizes_sq': [4, 6, 8], 'sizes_hy': [4, 6]},
    {'q': 10, 'gc_sq': None, 'gc_hy': 0.684, 'sizes_sq': [4, 5, 6], 'sizes_hy': [4, 5, 6]},
]

# First find S_q g_c for q=10
print("\n--- Finding S_q Potts g_c for q=10 ---", flush=True)
q10 = 10
# Quick scan at n=4,6
g_scan = np.linspace(0.02, 0.15, 40)
gap_n4 = []
gap_n6 = []

print("  n=4 scan...", flush=True)
t0 = time.time()
for g in g_scan:
    gap = get_gap(4, q10, g, 'sq_potts')
    gap_n4.append(gap * 4)
print(f"  n=4 done ({time.time()-t0:.1f}s)", flush=True)

print("  n=6 timing...", flush=True)
t0 = time.time()
gap_test = get_gap(6, q10, 0.08, 'sq_potts')
dt = time.time() - t0
dim6 = q10 ** 6
print(f"  n=6 (dim={dim6:,}): {dt:.1f}s per point", flush=True)

if dt < 8:  # Can do scan in ~5 min
    print("  n=6 scan...", flush=True)
    t0 = time.time()
    for g in g_scan:
        gap = get_gap(6, q10, g, 'sq_potts')
        gap_n6.append(gap * 6)
    print(f"  n=6 done ({time.time()-t0:.1f}s)", flush=True)

    # Find crossing
    diff = np.array(gap_n4) - np.array(gap_n6)
    gc_q10 = None
    for j in range(len(diff) - 1):
        if diff[j] * diff[j + 1] < 0:
            gc_q10 = float(g_scan[j] - diff[j] * (g_scan[j + 1] - g_scan[j]) / (diff[j + 1] - diff[j]))
            print(f"  S_q q=10 crossing (4,6): g_c = {gc_q10:.5f}", flush=True)
            break
    if gc_q10 is None:
        print("  No crossing found for q=10", flush=True)
    configs[2]['gc_sq'] = gc_q10
else:
    print(f"  n=6 too slow ({dt:.1f}s/pt), using n=4,5", flush=True)
    # Try n=5
    gap_n5 = []
    for g in g_scan:
        gap = get_gap(5, q10, g, 'sq_potts')
        gap_n5.append(gap * 5)
    diff = np.array(gap_n4) - np.array(gap_n5)
    gc_q10 = None
    for j in range(len(diff) - 1):
        if diff[j] * diff[j + 1] < 0:
            gc_q10 = float(g_scan[j] - diff[j] * (g_scan[j + 1] - g_scan[j]) / (diff[j + 1] - diff[j]))
            print(f"  S_q q=10 crossing (4,5): g_c = {gc_q10:.5f}", flush=True)
            break
    configs[2]['gc_sq'] = gc_q10

results['gc_sq_q10'] = configs[2]['gc_sq']

# Now measure Δ·N at g_c for each q
print(f"\n{'='*60}", flush=True)
print("Δ·N AT CRITICAL POINTS", flush=True)
print(f"{'='*60}", flush=True)

for cfg in configs:
    q = cfg['q']
    for model, gc, sizes, label in [
        ('sq_potts', cfg['gc_sq'], cfg['sizes_sq'], 'S_q Potts'),
        ('hybrid', cfg['gc_hy'], cfg['sizes_hy'], 'Hybrid'),
    ]:
        if gc is None:
            print(f"\n  {label} q={q}: g_c unknown, skipping", flush=True)
            continue

        print(f"\n  {label} q={q}, g_c={gc}:", flush=True)
        gapN_list = []
        for n in sizes:
            dim = q ** n
            if dim > 10_000_000:
                print(f"    n={n} (dim={dim:,}): too large, skipping", flush=True)
                continue
            t0 = time.time()
            gap = get_gap(n, q, gc, model)
            dt = time.time() - t0
            gapN = gap * n
            gapN_list.append((n, gapN))
            print(f"    n={n} (dim={dim:,}): Δ·N = {gapN:.5f} ({dt:.1f}s)", flush=True)

        # Compute drift
        if len(gapN_list) >= 2:
            first_n, first_gN = gapN_list[0]
            last_n, last_gN = gapN_list[-1]
            drift_pct = (last_gN - first_gN) / first_gN * 100
            print(f"    Drift n={first_n}→{last_n}: {drift_pct:+.1f}%", flush=True)
            results[f'{model}_q{q}_gapN'] = {
                'gc': gc,
                'data': [{'n': n, 'gap_N': float(gN)} for n, gN in gapN_list],
                'drift_pct': float(drift_pct),
            }

# Summary comparison
print(f"\n{'='*60}", flush=True)
print("SUMMARY — Δ·N DRIFT COMPARISON", flush=True)
print(f"{'='*60}", flush=True)
print(f"{'Model':<12} {'q':<4} {'n range':<10} {'Δ·N first':<12} {'Δ·N last':<12} {'Drift %':<10} {'Verdict'}", flush=True)
print("-" * 80, flush=True)

for key in sorted(results.keys()):
    if '_gapN' in key:
        d = results[key]
        data = d['data']
        model_label = key.split('_q')[0].replace('_', ' ')
        q_val = key.split('_q')[1].split('_')[0]
        n_range = f"{data[0]['n']}-{data[-1]['n']}"
        verdict = "continuous" if abs(d['drift_pct']) < 20 else "FIRST-ORDER?"
        if d['drift_pct'] < -30:
            verdict = "FIRST-ORDER"
        print(f"{model_label:<12} {q_val:<4} {n_range:<10} {data[0]['gap_N']:<12.5f} {data[-1]['gap_N']:<12.5f} {d['drift_pct']:<+10.1f} {verdict}", flush=True)

# Save
with open("results/sprint_077c_gap_scaling.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_077c_gap_scaling.json", flush=True)

from db_utils import record
for key in results:
    if '_gapN' in key:
        d = results[key]
        parts = key.split('_q')
        model = parts[0]
        q_val = int(parts[1].split('_')[0])
        data = d['data']
        record(sprint=77, model=model, q=q_val, n=data[-1]['n'],
               quantity='gapN_drift', value=d['drift_pct'],
               method='gap_scaling', notes=f"n={data[0]['n']}-{data[-1]['n']}, gc={d['gc']}")
print("Recorded to DB.", flush=True)
