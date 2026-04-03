#!/usr/bin/env python3
"""Sprint 076a: S_q Potts model gap crossings at q=5.

Build the STANDARD S_q Potts model:
  H = -Σ δ(s_i,s_j) - g Σ_{k=1}^{q-1} X^k

This has FULL S_q symmetry (permutation group), not just Z_q.
For q>4, this should be FIRST-ORDER (Gorbenko et al. 2018).

Compare gap×N crossings to our hybrid model (known g_c=0.441).
"""
import numpy as np
from scipy.sparse import csr_matrix, kron as sp_kron, eye as sp_eye, diags
from gpu_utils import eigsh
import json, time

q = 5


def sq_potts_hamiltonian_periodic(n, q, g):
    """Build S_q Potts Hamiltonian with periodic BC.

    Coupling: -δ(s_i, s_j)  (same as hybrid)
    Field: -g Σ_{k=1}^{q-1} X^k  (full S_q symmetric, not just X+X†)
    """
    dim = q ** n

    # Potts coupling (same as hybrid)
    potts_2site = np.zeros(q ** 2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q ** 2, q ** 2), format='csr')

    # S_q symmetric field: Σ_{k=1}^{q-1} X^k
    # X is the cyclic shift: X|s⟩ = |s+1 mod q⟩
    # X^k|s⟩ = |s+k mod q⟩
    # Σ_{k=1}^{q-1} X^k = (all-ones matrix) - I  (since X^0=I)
    # This is the S_q invariant: applies all non-identity permutations equally
    field_1site = np.ones((q, q)) - np.eye(q)
    field_op = csr_matrix(field_1site)

    H = csr_matrix((dim, dim))

    # Nearest-neighbor coupling (periodic)
    for i in range(n - 1):
        left = q ** i
        right = q ** (n - i - 2)
        op = potts_op
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - op

    # Periodic boundary (site n-1 to site 0)
    diag_boundary = np.zeros(dim)
    for idx in range(dim):
        s0 = idx % q
        sn = (idx // (q ** (n - 1))) % q
        diag_boundary[idx] = 1.0 if s0 == sn else 0.0
    H = H - diags(diag_boundary, 0, shape=(dim, dim), format='csr')

    # Transverse field: -g Σ_i Σ_{k=1}^{q-1} X^k_i
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


def get_gap(n, q, g, model='sq_potts'):
    """Get energy gap for S_q Potts or hybrid model."""
    if model == 'sq_potts':
        H = sq_potts_hamiltonian_periodic(n, q, g)
    else:
        H = hybrid_hamiltonian_periodic(n, q, g)
    dim = H.shape[0]
    if dim < 500:
        evals = np.linalg.eigvalsh(H.toarray())
        return evals[1] - evals[0], evals[:min(6, len(evals))]
    k = min(6, dim - 1)
    evals, _ = eigsh(H, k=k, which='SA')
    evals = np.sort(evals)
    return evals[1] - evals[0], evals


def hybrid_hamiltonian_periodic(n, q, g):
    """Hybrid model: Potts coupling + Z_q clock field (X + X†)."""
    dim = q ** n
    potts_2site = np.zeros(q ** 2)
    for a in range(q):
        potts_2site[a * q + a] = 1.0
    potts_op = diags(potts_2site, 0, shape=(q ** 2, q ** 2), format='csr')
    X = np.zeros((q, q))
    for s in range(q):
        X[(s + 1) % q, s] = 1.0
    XpXd = csr_matrix(X + X.T)
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
        op = XpXd
        if left > 1:
            op = sp_kron(sp_eye(left), op, format='csr')
        if right > 1:
            op = sp_kron(op, sp_eye(right), format='csr')
        H = H - g * op
    return H


# ---- MAIN ----
print("=" * 60, flush=True)
print("Sprint 076a: S_q Potts gap×N crossings at q=5", flush=True)
print("=" * 60, flush=True)

results = {
    'experiment': '076a_sq_potts_gap_q5',
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    'q': q,
}

# First verify: at q=3, S_q and hybrid should be IDENTICAL
print("\n--- Verification: q=3, S_q ≡ hybrid ---", flush=True)
for n_test in [4, 6]:
    g_test = 0.333
    gap_sq, _ = get_gap(n_test, 3, g_test, 'sq_potts')
    gap_hy, _ = get_gap(n_test, 3, g_test, 'hybrid')
    match = abs(gap_sq - gap_hy) < 1e-10
    print(f"  q=3 n={n_test}: S_q gap={gap_sq:.10f}, hybrid gap={gap_hy:.10f}, match={match}", flush=True)
    results[f'verify_q3_n{n_test}'] = {'sq_gap': float(gap_sq), 'hybrid_gap': float(gap_hy), 'match': match}

# Timing test for q=5
print("\n--- Timing tests for q=5 ---", flush=True)
for n_test in [4, 6, 8]:
    dim = q ** n_test
    print(f"  n={n_test} (dim={dim:,})...", end=" ", flush=True)
    t0 = time.time()
    gap, _ = get_gap(n_test, q, 0.5, 'sq_potts')
    dt = time.time() - t0
    print(f"{dt:.1f}s", flush=True)
    results[f'timing_n{n_test}'] = {'dim': dim, 'time_s': float(dt)}
    if dt > 200:
        print(f"  TOO SLOW — skipping n={n_test} in scan", flush=True)
        break

# Determine feasible sizes
sizes = []
for n in [4, 6, 8]:
    key = f'timing_n{n}'
    if key in results and results[key]['time_s'] < 200:
        sizes.append(n)

print(f"\nFeasible sizes: {sizes}", flush=True)

# Gap×N scan for S_q Potts at q=5
# Hybrid g_c=0.441, so S_q should be in similar range
# S_q Potts field is (q-1) times stronger per unit g than hybrid's 2cos term
# So S_q g_c should be SMALLER than hybrid g_c
# Rough estimate: hybrid field strength = 2cos(2π/5) ≈ 0.618 per unit g
# S_q field strength = (q-1) = 4 per unit g
# Ratio ≈ 0.618/4 ≈ 0.15, so S_q g_c ≈ 0.441 * 0.15 ≈ 0.068
# But this is rough — scan wider

g_range = np.linspace(0.01, 0.25, 50)

print(f"\n--- S_q Potts gap×N scan at q={q} ---", flush=True)
print(f"g range: [{g_range[0]:.3f}, {g_range[-1]:.3f}], {len(g_range)} points", flush=True)

gap_data_sq = {}
for n in sizes:
    dim = q ** n
    t0 = time.time()
    data = []
    for gi, g in enumerate(g_range):
        gap, evals = get_gap(n, q, g, 'sq_potts')
        data.append({
            'g': round(float(g), 5),
            'gap': float(gap),
            'gap_N': float(gap * n),
            'E0': float(evals[0]),
        })
        if dim > 10000 and (gi + 1) % 10 == 0:
            elapsed = time.time() - t0
            est = elapsed * len(g_range) / (gi + 1)
            print(f"  n={n}: {gi+1}/{len(g_range)} ({elapsed:.0f}s, est {est:.0f}s)", flush=True)
    dt = time.time() - t0
    gN = [d['gap_N'] for d in data]
    print(f"n={n} (dim={dim:,}): {dt:.1f}s, gap×N=[{min(gN):.4f}, {max(gN):.4f}]", flush=True)
    gap_data_sq[n] = data

results['gap_data_sq'] = {str(k): v for k, v in gap_data_sq.items()}

# Find crossings
print("\n--- S_q Potts crossings ---", flush=True)
crossings_sq = []
n_sorted = sorted(gap_data_sq.keys())
for i in range(len(n_sorted) - 1):
    n1, n2 = n_sorted[i], n_sorted[i + 1]
    g1 = np.array([d['g'] for d in gap_data_sq[n1]])
    v1 = np.array([d['gap_N'] for d in gap_data_sq[n1]])
    g2 = np.array([d['g'] for d in gap_data_sq[n2]])
    v2 = np.array([d['gap_N'] for d in gap_data_sq[n2]])
    g_common = np.linspace(max(g1[0], g2[0]), min(g1[-1], g2[-1]), 1000)
    iv1 = np.interp(g_common, g1, v1)
    iv2 = np.interp(g_common, g2, v2)
    diff = iv1 - iv2
    found = False
    for j in range(len(diff) - 1):
        if diff[j] * diff[j + 1] < 0:
            g_cross = g_common[j] - diff[j] * (g_common[j + 1] - g_common[j]) / (diff[j + 1] - diff[j])
            gN_at = float(np.interp(g_cross, g1, v1))
            print(f"  ({n1},{n2}): g_c = {g_cross:.5f}, gap×N = {gN_at:.4f}", flush=True)
            crossings_sq.append({
                'n_pair': [n1, n2], 'g_c': float(g_cross), 'gap_N_at': gN_at
            })
            found = True
            break
    if not found:
        print(f"  ({n1},{n2}): NO crossing", flush=True)
        # Check if curves approach each other
        min_diff = min(abs(diff))
        closest_idx = np.argmin(abs(diff))
        print(f"    Closest approach: |diff|={min_diff:.6f} at g={g_common[closest_idx]:.4f}", flush=True)
        crossings_sq.append({
            'n_pair': [n1, n2], 'g_c': None, 'min_diff': float(min_diff),
            'closest_g': float(g_common[closest_idx])
        })

results['crossings_sq'] = crossings_sq

# Also scan hybrid for comparison
print(f"\n--- Hybrid gap×N scan at q={q} (for comparison) ---", flush=True)
g_range_hy = np.linspace(0.2, 0.7, 50)
gap_data_hy = {}
for n in sizes:
    dim = q ** n
    t0 = time.time()
    data = []
    for gi, g in enumerate(g_range_hy):
        gap, evals = get_gap(n, q, g, 'hybrid')
        data.append({
            'g': round(float(g), 5),
            'gap': float(gap),
            'gap_N': float(gap * n),
            'E0': float(evals[0]),
        })
    dt = time.time() - t0
    gN = [d['gap_N'] for d in data]
    print(f"n={n}: {dt:.1f}s, gap×N=[{min(gN):.4f}, {max(gN):.4f}]", flush=True)
    gap_data_hy[n] = data

results['gap_data_hybrid'] = {str(k): v for k, v in gap_data_hy.items()}

# Find hybrid crossings
print("\n--- Hybrid crossings ---", flush=True)
crossings_hy = []
for i in range(len(n_sorted) - 1):
    n1, n2 = n_sorted[i], n_sorted[i + 1]
    g1 = np.array([d['g'] for d in gap_data_hy[n1]])
    v1 = np.array([d['gap_N'] for d in gap_data_hy[n1]])
    g2 = np.array([d['g'] for d in gap_data_hy[n2]])
    v2 = np.array([d['gap_N'] for d in gap_data_hy[n2]])
    g_common = np.linspace(max(g1[0], g2[0]), min(g1[-1], g2[-1]), 1000)
    iv1 = np.interp(g_common, g1, v1)
    iv2 = np.interp(g_common, g2, v2)
    diff = iv1 - iv2
    found = False
    for j in range(len(diff) - 1):
        if diff[j] * diff[j + 1] < 0:
            g_cross = g_common[j] - diff[j] * (g_common[j + 1] - g_common[j]) / (diff[j + 1] - diff[j])
            gN_at = float(np.interp(g_cross, g1, v1))
            print(f"  ({n1},{n2}): g_c = {g_cross:.5f}, gap×N = {gN_at:.4f}", flush=True)
            crossings_hy.append({
                'n_pair': [n1, n2], 'g_c': float(g_cross), 'gap_N_at': gN_at
            })
            found = True
            break
    if not found:
        print(f"  ({n1},{n2}): NO crossing", flush=True)

results['crossings_hybrid'] = crossings_hy

# Summary
print(f"\n{'='*60}", flush=True)
print("SUMMARY", flush=True)
print(f"{'='*60}", flush=True)

if crossings_sq:
    gc_sq = [c['g_c'] for c in crossings_sq if c.get('g_c') is not None]
    if gc_sq:
        print(f"S_q Potts g_c(q=5): {gc_sq}", flush=True)
    else:
        print("S_q Potts: NO gap×N crossing found — possible first-order signature!", flush=True)

if crossings_hy:
    gc_hy = [c['g_c'] for c in crossings_hy if c.get('g_c') is not None]
    if gc_hy:
        print(f"Hybrid g_c(q=5): {gc_hy}", flush=True)

print(f"\nS_q field (1-site): all-ones - I  (strength {q-1} per unit g)", flush=True)
print(f"Hybrid field (1-site): X + X†  (strength 2cos(2π/{q}) = {2*np.cos(2*np.pi/q):.4f} per unit g)", flush=True)
print(f"Field ratio (S_q/hybrid): {(q-1) / (2*np.cos(2*np.pi/q)):.3f}", flush=True)

# Save
with open("results/sprint_076a_sq_potts_gap_q5.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print("\nSaved to results/sprint_076a_sq_potts_gap_q5.json", flush=True)

from db_utils import record
for c in crossings_sq:
    if c.get('g_c') is not None:
        record(sprint=76, model='sq_potts', q=q, n=max(c['n_pair']),
               quantity='gc_gap', value=c['g_c'],
               method='gapN_crossing', notes=f"pair {c['n_pair']}")
for c in crossings_hy:
    if c.get('g_c') is not None:
        record(sprint=76, model='hybrid', q=q, n=max(c['n_pair']),
               quantity='gc_gap', value=c['g_c'],
               method='gapN_crossing', notes=f"pair {c['n_pair']}, verification")
print("Recorded to DB.", flush=True)
