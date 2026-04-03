#!/usr/bin/env python3
"""Sprint 077a: S_q Potts model gap crossings at q=7.

Sprint 076c found 69% gap collapse at q=7 n=4→6 at estimated g_c≈0.12.
Need proper gap×N crossing scan to find true g_c and test first-order.

S_q Potts: H = -Σ δ(s_i,s_j) - g Σ_{k=1}^{q-1} X^k  (full S_q symmetry)
Field strength per unit g: q-1 = 6 for q=7.

Sizes: n=4 (dim=2401), n=6 (dim=117649), n=8 (dim=5764801 ~5.8M, GPU).
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from gpu_utils import eigsh
import json, time

q = 7


def sq_potts_hamiltonian_periodic(n, q, g):
    """Build S_q Potts Hamiltonian with periodic BC."""
    dim = q ** n
    # Potts coupling
    potts_2site = np.zeros(q ** 2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q ** 2, q ** 2), format='csr')
    # S_q field: all-ones - I
    field_1site = np.ones((q, q)) - np.eye(q)
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

    # Periodic boundary
    diag_boundary = np.zeros(dim)
    for idx in range(dim):
        s0 = idx % q
        sn = (idx // (q ** (n - 1))) % q
        diag_boundary[idx] = 1.0 if s0 == sn else 0.0
    H = H - diags(diag_boundary, 0, shape=(dim, dim), format='csr')

    # S_q field
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


def get_gap(n, q, g):
    """Get energy gap and first few eigenvalues."""
    H = sq_potts_hamiltonian_periodic(n, q, g)
    dim = H.shape[0]
    if dim < 500:
        evals = np.linalg.eigvalsh(H.toarray())
        return evals[1] - evals[0], evals[:min(8, len(evals))]
    k = min(8, dim - 1)
    evals, _ = eigsh(H, k=k, which='SA')
    evals = np.sort(evals)
    return evals[1] - evals[0], evals


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 077a: S_q Potts gap×N crossings at q=7", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '077a_sq_potts_gap_q7',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q,
}

# Timing test
print("\n--- Timing tests ---", flush=True)
feasible = []
for n_test in [4, 6, 8]:
    dim = q ** n_test
    print(f"  n={n_test} (dim={dim:,})...", end=" ", flush=True)
    t0 = time.time()
    gap, evals = get_gap(n_test, q, 0.1)
    dt = time.time() - t0
    print(f"{dt:.1f}s, gap={gap:.6f}", flush=True)
    results[f'timing_n{n_test}'] = {'dim': dim, 'time_s': float(dt), 'gap': float(gap)}
    if dt < 60:
        feasible.append(n_test)
    else:
        print(f"  Too slow for scan, skipping n={n_test}", flush=True)
        break

print(f"\nFeasible sizes for scan: {feasible}", flush=True)

# Gap×N scan
# S_q field strength = q-1 = 6, hybrid field = 2cos(2π/7) ≈ 1.247
# Hybrid g_c(q=7) = 0.535, so S_q g_c ≈ 0.535 * 1.247/6 ≈ 0.111
# Sprint 076c estimated g_c≈0.12 with 69% collapse
# Scan wider: 0.02 to 0.25
g_range = np.linspace(0.02, 0.25, 60)

print(f"\n--- S_q Potts gap×N scan at q={q} ---", flush=True)
print(f"g range: [{g_range[0]:.3f}, {g_range[-1]:.3f}], {len(g_range)} points", flush=True)

gap_data = {}
for n in feasible:
    dim = q ** n
    t0 = time.time()
    data = []
    for gi, g in enumerate(g_range):
        gap, evals = get_gap(n, q, g)
        data.append({
            'g': round(float(g), 5),
            'gap': float(gap),
            'gap_N': float(gap * n),
            'E0': float(evals[0]),
            'E1': float(evals[1]),
            'E2': float(evals[2]) if len(evals) > 2 else None,
        })
        if dim > 10000 and (gi + 1) % 10 == 0:
            elapsed = time.time() - t0
            est = elapsed * len(g_range) / (gi + 1)
            print(f"  n={n}: {gi+1}/{len(g_range)} ({elapsed:.0f}s, est {est:.0f}s)", flush=True)
    dt = time.time() - t0
    gN = [d['gap_N'] for d in data]
    gaps = [d['gap'] for d in data]
    print(f"n={n} (dim={dim:,}): {dt:.1f}s, gap×N=[{min(gN):.4f}, {max(gN):.4f}], min gap={min(gaps):.6f}", flush=True)
    gap_data[n] = data

results['gap_data'] = {str(k): v for k, v in gap_data.items()}

# Find crossings
print("\n--- Gap×N crossings ---", flush=True)
crossings = []
n_sorted = sorted(gap_data.keys())
for i in range(len(n_sorted) - 1):
    n1, n2 = n_sorted[i], n_sorted[i + 1]
    g1 = np.array([d['g'] for d in gap_data[n1]])
    v1 = np.array([d['gap_N'] for d in gap_data[n1]])
    g2 = np.array([d['g'] for d in gap_data[n2]])
    v2 = np.array([d['gap_N'] for d in gap_data[n2]])
    g_common = np.linspace(max(g1[0], g2[0]), min(g1[-1], g2[-1]), 2000)
    iv1 = np.interp(g_common, g1, v1)
    iv2 = np.interp(g_common, g2, v2)
    diff = iv1 - iv2
    found = False
    for j in range(len(diff) - 1):
        if diff[j] * diff[j + 1] < 0:
            g_cross = g_common[j] - diff[j] * (g_common[j + 1] - g_common[j]) / (diff[j + 1] - diff[j])
            gN_at = float(np.interp(g_cross, g1, v1))
            print(f"  ({n1},{n2}): g_c = {g_cross:.5f}, gap×N = {gN_at:.4f}", flush=True)
            crossings.append({
                'n_pair': [n1, n2], 'g_c': float(g_cross), 'gap_N_at': gN_at
            })
            found = True
            break
    if not found:
        print(f"  ({n1},{n2}): NO crossing found", flush=True)
        min_diff = float(min(abs(diff)))
        closest_idx = int(np.argmin(abs(diff)))
        print(f"    Closest approach: |diff|={min_diff:.6f} at g={g_common[closest_idx]:.4f}", flush=True)
        # Check sign of diff — does larger n have bigger or smaller gap×N?
        avg_diff = float(np.mean(diff))
        print(f"    Average diff (n1-n2): {avg_diff:.4f} ({'n1 larger' if avg_diff > 0 else 'n2 larger'})", flush=True)
        crossings.append({
            'n_pair': [n1, n2], 'g_c': None, 'min_diff': min_diff,
            'closest_g': float(g_common[closest_idx]), 'avg_diff': avg_diff
        })

results['crossings'] = crossings

# Check for first-order signatures: gap collapse and level crossing
print("\n--- First-order diagnostics ---", flush=True)
for n in feasible:
    data = gap_data[n]
    min_gap = min(d['gap'] for d in data)
    min_gap_g = [d['g'] for d in data if d['gap'] == min_gap][0]
    # Gap×N at minimum
    min_gapN = min_gap * n
    print(f"  n={n}: min gap = {min_gap:.6f} at g={min_gap_g:.4f}, gap×N={min_gapN:.4f}", flush=True)
    results[f'first_order_n{n}'] = {
        'min_gap': float(min_gap),
        'min_gap_g': float(min_gap_g),
        'min_gap_N': float(min_gapN),
    }

# If we have 2+ sizes, check gap×N at minimum
if len(feasible) >= 2:
    print("\n  Gap×N at min gap:", flush=True)
    for n in feasible:
        fo = results[f'first_order_n{n}']
        print(f"    n={n}: gap×N_min = {fo['min_gap_N']:.4f}", flush=True)
    gNmin_vals = [results[f'first_order_n{n}']['min_gap_N'] for n in feasible]
    if len(gNmin_vals) >= 2:
        trend = "DECREASING (first-order)" if gNmin_vals[-1] < gNmin_vals[0] else "INCREASING/STABLE (continuous)"
        print(f"  Trend: {trend}", flush=True)
        results['gap_N_min_trend'] = trend

# Summary
print(f"\n{'='*60}", flush=True)
print("SUMMARY", flush=True)
print(f"{'='*60}", flush=True)
print(f"q={q}, sizes tested: {feasible}", flush=True)
gc_vals = [c['g_c'] for c in crossings if c.get('g_c') is not None]
if gc_vals:
    print(f"S_q Potts g_c(q=7) from crossings: {gc_vals}", flush=True)
    print(f"Best estimate: {np.mean(gc_vals):.5f}", flush=True)
else:
    print("NO gap×N crossing found — possible first-order signature", flush=True)

# Save
with open("results/sprint_077a_sq_potts_gap_q7.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_077a_sq_potts_gap_q7.json", flush=True)

from db_utils import record
for c in crossings:
    if c.get('g_c') is not None:
        record(sprint=77, model='sq_potts', q=q, n=max(c['n_pair']),
               quantity='gc_gap', value=c['g_c'],
               method='gapN_crossing', notes=f"pair {c['n_pair']}")
    else:
        record(sprint=77, model='sq_potts', q=q, n=max(c['n_pair']),
               quantity='gc_gap', value=-1,
               method='gapN_crossing', notes=f"NO crossing, pair {c['n_pair']}, min_diff={c.get('min_diff',0):.6f}")
print("Recorded to DB.", flush=True)
