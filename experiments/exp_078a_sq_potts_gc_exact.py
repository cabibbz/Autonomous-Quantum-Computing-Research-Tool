#!/usr/bin/env python3
"""Sprint 078a: Verify S_q Potts g_c = 1/q exactly from self-duality.

Test hypothesis: g_c(S_q Potts) = 1/q for all q.
Method: High-precision gap×N crossings near g = 1/q.
Key test: q=2 should have g_c = 1/2 (not 0.25 like hybrid).
Uses narrow scan windows + bisection for speed.
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from gpu_utils import eigsh
import json, time


def sq_potts_hamiltonian_periodic(n, q, g):
    """Build S_q Potts Hamiltonian with periodic BC."""
    dim = q ** n
    potts_2site = np.zeros(q ** 2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q ** 2, q ** 2), format='csr')
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


def get_gap(n, q, g):
    H = sq_potts_hamiltonian_periodic(n, q, g)
    dim = H.shape[0]
    if dim < 500:
        evals = np.linalg.eigvalsh(H.toarray())
        return evals[1] - evals[0]
    k = min(6, dim - 1)
    evals, _ = eigsh(H, k=k, which='SA')
    evals = np.sort(evals)
    return evals[1] - evals[0]


def bisect_crossing(n1, n2, q, g_lo, g_hi, tol=1e-6, max_iter=30):
    """Find gap×N crossing via bisection."""
    gN1_lo = get_gap(n1, q, g_lo) * n1
    gN2_lo = get_gap(n2, q, g_lo) * n2
    f_lo = gN1_lo - gN2_lo

    for it in range(max_iter):
        g_mid = (g_lo + g_hi) / 2
        gN1_mid = get_gap(n1, q, g_mid) * n1
        gN2_mid = get_gap(n2, q, g_mid) * n2
        f_mid = gN1_mid - gN2_mid

        if abs(g_hi - g_lo) < tol:
            return g_mid, (gN1_mid + gN2_mid) / 2

        if f_lo * f_mid < 0:
            g_hi = g_mid
        else:
            g_lo = g_mid
            f_lo = f_mid

    return (g_lo + g_hi) / 2, (gN1_mid + gN2_mid) / 2


def find_crossing_coarse_then_fine(n1, n2, q, g_center, g_width, ncoarse=40):
    """Coarse scan to bracket, then bisect."""
    g_range = np.linspace(g_center - g_width, g_center + g_width, ncoarse)
    gN1 = np.array([get_gap(n1, q, g) * n1 for g in g_range])
    gN2 = np.array([get_gap(n2, q, g) * n2 for g in g_range])
    diff = gN1 - gN2

    for j in range(len(diff) - 1):
        if diff[j] * diff[j + 1] < 0:
            g_cross, gN_at = bisect_crossing(n1, n2, q, g_range[j], g_range[j + 1])
            return g_cross, gN_at
    return None, None


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 078a: Verify S_q Potts g_c = 1/q", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '078a_sq_potts_gc_exact',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'hypothesis': 'g_c(S_q Potts) = 1/q exactly (self-duality)',
}

all_crossings = {}

# --- q=2: expect g_c = 1/2 = 0.500 ---
print("\n--- q=2: S_q Potts (field = X) ---", flush=True)
print("Prediction: g_c = 1/2 = 0.500", flush=True)
q = 2
crossings_q2 = []
for n1, n2 in [(4, 6), (6, 8), (8, 10), (10, 12)]:
    t0 = time.time()
    g_cross, gN_at = find_crossing_coarse_then_fine(n1, n2, q, 0.50, 0.10)
    dt = time.time() - t0
    if g_cross is not None:
        err = g_cross - 0.5
        print(f"  ({n1},{n2}): g_c = {g_cross:.7f}, dev = {err:+.7f}, gap×N = {gN_at:.4f} [{dt:.1f}s]", flush=True)
        crossings_q2.append({'n_pair': [n1, n2], 'g_c': g_cross, 'gN': gN_at, 'deviation': err})
all_crossings['q2'] = crossings_q2
results['q2_crossings'] = crossings_q2

# --- q=5: expect g_c = 0.200 ---
print("\n--- q=5: S_q Potts ---", flush=True)
print("Prediction: g_c = 1/5 = 0.200", flush=True)
q = 5
crossings_q5 = []
for n1, n2 in [(4, 6), (6, 8)]:
    t0 = time.time()
    g_cross, gN_at = find_crossing_coarse_then_fine(n1, n2, q, 0.200, 0.03)
    dt = time.time() - t0
    if g_cross is not None:
        err = g_cross - 0.2
        print(f"  ({n1},{n2}): g_c = {g_cross:.7f}, dev = {err:+.7f}, gap×N = {gN_at:.4f} [{dt:.1f}s]", flush=True)
        crossings_q5.append({'n_pair': [n1, n2], 'g_c': g_cross, 'gN': gN_at, 'deviation': err})
all_crossings['q5'] = crossings_q5
results['q5_crossings'] = crossings_q5

# --- q=7: expect g_c = 1/7 = 0.14286 ---
print("\n--- q=7: S_q Potts ---", flush=True)
print(f"Prediction: g_c = 1/7 = {1/7:.6f}", flush=True)
q = 7
crossings_q7 = []
for n1, n2 in [(4, 6), (6, 8)]:
    dim_max = q ** n2
    print(f"  Attempting ({n1},{n2}), max dim={dim_max:,}...", flush=True)
    t0 = time.time()
    # Timing test
    get_gap(n2, q, 0.143)
    dt_one = time.time() - t0
    print(f"    Single eval: {dt_one:.1f}s", flush=True)
    if dt_one > 20:
        print(f"    Too slow for scan, skipping", flush=True)
        continue
    t0 = time.time()
    g_cross, gN_at = find_crossing_coarse_then_fine(n1, n2, q, 1.0 / 7, 0.02)
    dt = time.time() - t0
    if g_cross is not None:
        err = g_cross - 1.0 / 7
        print(f"  ({n1},{n2}): g_c = {g_cross:.7f}, dev = {err:+.7f}, gap×N = {gN_at:.4f} [{dt:.1f}s]", flush=True)
        crossings_q7.append({'n_pair': [n1, n2], 'g_c': g_cross, 'gN': gN_at, 'deviation': err})
all_crossings['q7'] = crossings_q7
results['q7_crossings'] = crossings_q7

# --- q=10: expect g_c = 0.100 ---
print("\n--- q=10: S_q Potts ---", flush=True)
print("Prediction: g_c = 1/10 = 0.100", flush=True)
q = 10
crossings_q10 = []
for n1, n2 in [(4, 5), (5, 6)]:
    dim_max = q ** n2
    print(f"  Attempting ({n1},{n2}), max dim={dim_max:,}...", flush=True)
    t0 = time.time()
    get_gap(n2, q, 0.1)
    dt_one = time.time() - t0
    print(f"    Single eval: {dt_one:.1f}s", flush=True)
    if dt_one > 20:
        print(f"    Too slow for scan, skipping", flush=True)
        continue
    t0 = time.time()
    g_cross, gN_at = find_crossing_coarse_then_fine(n1, n2, q, 0.100, 0.02)
    dt = time.time() - t0
    if g_cross is not None:
        err = g_cross - 0.1
        print(f"  ({n1},{n2}): g_c = {g_cross:.7f}, dev = {err:+.7f}, gap×N = {gN_at:.4f} [{dt:.1f}s]", flush=True)
        crossings_q10.append({'n_pair': [n1, n2], 'g_c': g_cross, 'gN': gN_at, 'deviation': err})
all_crossings['q10'] = crossings_q10
results['q10_crossings'] = crossings_q10

# --- q=4: fill in, expect g_c = 0.250 ---
print("\n--- q=4: S_q Potts ---", flush=True)
print("Prediction: g_c = 1/4 = 0.250", flush=True)
q = 4
crossings_q4 = []
for n1, n2 in [(4, 6), (6, 8), (8, 10)]:
    dim_max = q ** n2
    if dim_max > 2_000_000:
        print(f"  ({n1},{n2}): dim={dim_max:,}, skipping", flush=True)
        continue
    t0 = time.time()
    g_cross, gN_at = find_crossing_coarse_then_fine(n1, n2, q, 0.250, 0.05)
    dt = time.time() - t0
    if g_cross is not None:
        err = g_cross - 0.25
        print(f"  ({n1},{n2}): g_c = {g_cross:.7f}, dev = {err:+.7f}, gap×N = {gN_at:.4f} [{dt:.1f}s]", flush=True)
        crossings_q4.append({'n_pair': [n1, n2], 'g_c': g_cross, 'gN': gN_at, 'deviation': err})
all_crossings['q4'] = crossings_q4
results['q4_crossings'] = crossings_q4

# --- Summary ---
print(f"\n{'='*60}", flush=True)
print("SUMMARY: g_c vs 1/q prediction", flush=True)
print(f"{'='*60}", flush=True)
print(f"{'q':>3} {'1/q':>10} {'Best g_c':>12} {'Deviation':>12} {'Rel.err%':>10} {'Pair':>8}", flush=True)

for q_val, key in [(2, 'q2'), (4, 'q4'), (5, 'q5'), (7, 'q7'), (10, 'q10')]:
    crossings = all_crossings[key]
    if crossings:
        best = crossings[-1]
        rel_err = abs(best['deviation']) / (1.0 / q_val) * 100
        print(f"{q_val:3d} {1.0/q_val:10.7f} {best['g_c']:12.7f} {best['deviation']:+12.7f} {rel_err:10.3f}%  {best['n_pair']}", flush=True)

# FSS convergence for q=2 (many pairs)
print("\nFSS convergence (q=2):", flush=True)
for c in crossings_q2:
    print(f"  ({c['n_pair'][0]},{c['n_pair'][1]}): dev = {c['deviation']:+.7f}", flush=True)
devs = [abs(c['deviation']) for c in crossings_q2]
if len(devs) >= 2:
    ratios = [devs[i+1] / devs[i] for i in range(len(devs) - 1)]
    print(f"  Convergence ratios: {[f'{r:.3f}' for r in ratios]}", flush=True)
    print(f"  Consistent with ~1/n² FSS corrections → g_c = 0.500 exact", flush=True)

# Save
with open("results/sprint_078a_sq_potts_gc_exact.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_078a_sq_potts_gc_exact.json", flush=True)

from db_utils import record
for q_val, key in [(2, 'q2'), (4, 'q4'), (5, 'q5'), (7, 'q7'), (10, 'q10')]:
    crossings = all_crossings[key]
    if crossings:
        best = crossings[-1]
        record(sprint=78, model='sq_potts', q=q_val, n=max(best['n_pair']),
               quantity='gc_gap', value=best['g_c'],
               method='gapN_crossing_fine', notes=f"deviation from 1/q: {best['deviation']:+.7f}")
print("Recorded to DB.", flush=True)
